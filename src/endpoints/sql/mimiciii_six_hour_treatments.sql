WITH assigned(icu_stay_id,trajectory_class) AS (VALUES %s), eligible AS MATERIALIZED (
        SELECT i.*,a.trajectory_class FROM assigned a JOIN mimiciii.icustays i
        ON i.icustay_id=a.icu_stay_id
    ), weight AS MATERIALIZED (
        SELECT DISTINCT ON(e.icustay_id) e.icustay_id,w.weight AS weight_raw_kg,
            CASE WHEN w.weight BETWEEN 20 AND 300 THEN w.weight END AS weight_kg
        FROM eligible e LEFT JOIN mimiciii.weight_durations w ON w.icustay_id=e.icustay_id
            AND w.starttime<e.intime+INTERVAL '24 hours' AND w.endtime>e.intime
        ORDER BY e.icustay_id,CASE WHEN w.starttime<=e.intime AND w.endtime>e.intime THEN 0 ELSE 1 END,w.starttime
    ), phenylephrine_cv_mcg_min AS MATERIALIZED (
    SELECT cv.icustay_id,cv.charttime,cv.rate
    FROM mimiciii.inputevents_cv AS cv INNER JOIN eligible AS e USING(icustay_id)
    WHERE cv.itemid=30127 AND cv.rate>0
      AND cv.charttime >= e.intime AND cv.charttime < e.intime + INTERVAL '24 hours'
), phenylephrine_segments AS MATERIALIZED (
    SELECT p.icustay_id,p.starttime,p.endtime,p.vaso_rate,
           BOOL_OR(cv.icustay_id IS NOT NULL) AS is_mcg_min
    FROM mimiciii.phenylephrine_dose AS p INNER JOIN eligible AS e USING(icustay_id)
    LEFT JOIN phenylephrine_cv_mcg_min AS cv
      ON cv.icustay_id=p.icustay_id AND cv.charttime>=p.starttime AND cv.charttime<p.endtime
     AND ABS(cv.rate-p.vaso_rate)<1e-6
    GROUP BY p.icustay_id,p.starttime,p.endtime,p.vaso_rate
), vasopressin_mv_stays AS MATERIALIZED (
    SELECT DISTINCT im.icustay_id FROM mimiciii.inputevents_mv AS im
    INNER JOIN eligible AS e USING(icustay_id)
    WHERE im.itemid=222315 AND im.starttime < e.intime + INTERVAL '24 hours'
      AND COALESCE(im.endtime, im.starttime) >= e.intime
), vasopressin_cv_u_hour AS MATERIALIZED (
    SELECT cv.icustay_id,cv.charttime,cv.rate
    FROM mimiciii.inputevents_cv AS cv INNER JOIN eligible AS e USING(icustay_id)
    WHERE cv.itemid=30051 AND LOWER(cv.rateuom)='uhr' AND cv.rate>0
      AND cv.charttime >= e.intime AND cv.charttime < e.intime + INTERVAL '24 hours'
), vasopressin_segments AS MATERIALIZED (
    SELECT p.icustay_id,p.starttime,p.endtime,p.vaso_rate,
           (mv.icustay_id IS NOT NULL OR BOOL_OR(cv.icustay_id IS NOT NULL)) AS is_u_hour
    FROM mimiciii.vasopressin_dose AS p INNER JOIN eligible AS e USING(icustay_id)
    LEFT JOIN vasopressin_mv_stays AS mv USING(icustay_id)
    LEFT JOIN vasopressin_cv_u_hour AS cv
      ON cv.icustay_id=p.icustay_id AND cv.charttime>=p.starttime AND cv.charttime<p.endtime
     AND ABS(cv.rate-p.vaso_rate)<1e-6
    GROUP BY p.icustay_id,p.starttime,p.endtime,p.vaso_rate,mv.icustay_id
), component_events AS MATERIALIZED (
    SELECT d.icustay_id, d.starttime, d.endtime, 'norepinephrine'::text AS drug, d.vaso_rate::numeric AS rate, 1.0::numeric AS factor
    FROM mimiciii.norepinephrine_dose d INNER JOIN eligible e USING(icustay_id)
    WHERE d.vaso_rate >= 0 AND d.starttime < e.intime + INTERVAL '24 hours' AND d.endtime > e.intime
    UNION ALL SELECT d.icustay_id,d.starttime,d.endtime,'epinephrine',d.vaso_rate::numeric,1.0 FROM mimiciii.epinephrine_dose d INNER JOIN eligible e USING(icustay_id)
    WHERE d.vaso_rate >= 0 AND d.starttime < e.intime + INTERVAL '24 hours' AND d.endtime > e.intime
    UNION ALL
    SELECT p.icustay_id,p.starttime,p.endtime,'phenylephrine',
           (CASE WHEN p.is_mcg_min THEN p.vaso_rate/NULLIF(w.weight_kg,0)
                 ELSE p.vaso_rate END)::numeric,0.1
    FROM phenylephrine_segments p
    LEFT JOIN weight w USING(icustay_id)
    WHERE p.vaso_rate >= 0
    UNION ALL SELECT d.icustay_id,d.starttime,d.endtime,'dopamine',d.vaso_rate::numeric,0.01 FROM mimiciii.dopamine_dose d INNER JOIN eligible e USING(icustay_id)
    WHERE d.vaso_rate >= 0 AND d.starttime < e.intime + INTERVAL '24 hours' AND d.endtime > e.intime
    UNION ALL SELECT icustay_id,starttime,endtime,'vasopressin',
      (CASE WHEN is_u_hour THEN vaso_rate/60.0 ELSE vaso_rate END)::numeric,2.5
      FROM vasopressin_segments WHERE vaso_rate >= 0
),
    bounds AS (
        SELECT icustay_id,intime AS t FROM eligible UNION
        SELECT icustay_id,intime+INTERVAL '24 hours' FROM eligible UNION
        SELECT n.icustay_id,GREATEST(n.starttime,e.intime) FROM component_events n JOIN eligible e USING(icustay_id)
        WHERE n.starttime<e.intime+INTERVAL '24 hours' AND n.endtime>e.intime UNION
        SELECT n.icustay_id,LEAST(n.endtime,e.intime+INTERVAL '24 hours') FROM component_events n JOIN eligible e USING(icustay_id)
        WHERE n.starttime<e.intime+INTERVAL '24 hours' AND n.endtime>e.intime
    ), segments AS (
        SELECT icustay_id,t AS starttime,LEAD(t) OVER(PARTITION BY icustay_id ORDER BY t) AS endtime FROM bounds
    ), per_drug AS (
        SELECT s.icustay_id,s.starttime,s.endtime,n.drug,MAX(n.rate*n.factor) AS dose,
            BOOL_OR(n.drug IS NOT NULL AND n.rate IS NULL) AS unquantified
        FROM segments s LEFT JOIN component_events n ON n.icustay_id=s.icustay_id
        AND n.starttime<=s.starttime AND n.endtime>s.starttime
        WHERE s.endtime>s.starttime GROUP BY s.icustay_id,s.starttime,s.endtime,n.drug
    ), nee AS (
        SELECT icustay_id,starttime,endtime,SUM(dose) AS dose,
            BOOL_OR(unquantified) AS unquantified FROM per_drug GROUP BY 1,2,3
    ), intervals AS (
        SELECT e.*,k AS interval_index,e.intime+(k-1)*INTERVAL '6 hours' AS t0,
            e.intime+k*INTERVAL '6 hours' AS t1,w.weight_raw_kg,w.weight_kg
        FROM eligible e CROSS JOIN generate_series(1,4) k
        LEFT JOIN weight w ON w.icustay_id=e.icustay_id
    ), input_events AS MATERIALIZED (
        SELECT im.icustay_id,im.starttime,COALESCE(im.endtime,im.starttime) AS endtime,im.amount,
            (LOWER(COALESCE(im.ordercategoryname,'')) SIMILAR TO '(02-%%|03-%%|04-%%)'
                AND (LOWER(COALESCE(di.label,'')) ~ 'albumin|hetastarch|dextran|saline|nacl|ringer|(^|[^a-z])lr([^a-z]|$)|lactated|dextrose|d5|sterile water|crystalloid'
                     OR LOWER(COALESCE(di.category,''))='fluids/intake')
                AND LOWER(COALESCE(di.label,'')) NOT IN ('piggyback','solution','multivitamins')
                AND LOWER(COALESCE(di.label,'')) !~ 'tube feed|replete|nepro|osmolite|jevity|promote|glucerna|nutren|packed red|platelet|fresh frozen|cryoprecip') AS audited,
            'MetaVision'::text AS source
        FROM eligible e JOIN mimiciii.inputevents_mv im USING(icustay_id)
        JOIN mimiciii.d_items di USING(itemid)
        WHERE im.starttime<e.intime+INTERVAL '24 hours' AND COALESCE(im.endtime,im.starttime)>=e.intime
            AND LOWER(im.amountuom)='ml' AND im.amount>0 AND im.cancelreason=0
            AND COALESCE(im.statusdescription,'') NOT IN ('Rewritten','Cancelled')
        UNION ALL
        SELECT cv.icustay_id,cv.charttime,cv.charttime,cv.amount,
            cv.itemid IN (30008,30009,30018,30020,30021,30143,30159,30160,30161,30176,30181,30352,30353,
                40548,42832,43237,43353,44203,44440,44952,45403,45532,46564),'CareVue'
        FROM eligible e JOIN mimiciii.inputevents_cv cv USING(icustay_id)
        WHERE cv.charttime>=e.intime AND cv.charttime<e.intime+INTERVAL '24 hours'
            AND LOWER(cv.amountuom)='ml' AND cv.amount>0
    ), input_interval AS (
        SELECT i.icustay_id,i.interval_index,COUNT(f.amount) AS input_events,
            SUM(CASE WHEN f.endtime>f.starttime THEN f.amount*EXTRACT(EPOCH FROM(LEAST(f.endtime,i.t1)-GREATEST(f.starttime,i.t0)))
                /EXTRACT(EPOCH FROM(f.endtime-f.starttime)) ELSE f.amount END) AS input_ml,
            COALESCE(SUM(CASE WHEN f.endtime>f.starttime THEN f.amount*EXTRACT(EPOCH FROM(LEAST(f.endtime,i.t1)-GREATEST(f.starttime,i.t0)))
                /EXTRACT(EPOCH FROM(f.endtime-f.starttime)) ELSE f.amount END) FILTER(WHERE f.audited),0) AS audited_iv_ml,
            COUNT(f.amount) FILTER(WHERE f.source='MetaVision') AS mv_input_events,
            COUNT(f.amount) FILTER(WHERE f.source='CareVue') AS cv_input_events
        FROM intervals i LEFT JOIN input_events f ON f.icustay_id=i.icustay_id AND f.starttime<i.t1
            AND ((f.endtime>f.starttime AND f.endtime>i.t0) OR (f.endtime=f.starttime AND f.starttime>=i.t0)) GROUP BY 1,2
    ), output_interval AS (
        SELECT i.icustay_id,i.interval_index,COUNT(o.value) AS output_events,SUM(o.value) AS output_ml
        FROM intervals i LEFT JOIN mimiciii.outputevents o ON o.icustay_id=i.icustay_id
            AND o.charttime>=i.t0 AND o.charttime<i.t1 AND LOWER(o.valueuom)='ml' AND o.value>0 AND COALESCE(o.iserror,0)=0 GROUP BY 1,2
    ), nee_interval AS (
        SELECT i.icustay_id,i.interval_index,CASE WHEN BOOL_OR(n.unquantified) THEN NULL ELSE COALESCE(MAX(n.dose),0) END AS nee_max,
            (COALESCE(MAX(n.dose),0)>0 OR COALESCE(BOOL_OR(n.unquantified),FALSE))::integer AS nee_any,
            COALESCE(BOOL_OR(n.unquantified),FALSE)::integer AS nee_unquantified
        FROM intervals i LEFT JOIN nee n ON n.icustay_id=i.icustay_id AND n.starttime<i.t1 AND n.endtime>i.t0 GROUP BY 1,2
    ), vitals AS (
        SELECT i.icustay_id,i.interval_index,COUNT(v.meanbp) FILTER(WHERE v.meanbp BETWEEN 20 AND 200) AS map_observed_n,
            percentile_cont(.5) WITHIN GROUP(ORDER BY v.meanbp) FILTER(WHERE v.meanbp BETWEEN 20 AND 200) AS map_median
        FROM intervals i LEFT JOIN mimiciii.pivoted_vital v ON v.icustay_id=i.icustay_id AND v.charttime>=i.t0 AND v.charttime<i.t1 GROUP BY 1,2
    ), lactate AS (
        SELECT i.icustay_id,i.interval_index,COUNT(l.valuenum) AS lactate_observed_n,
            percentile_cont(.5) WITHIN GROUP(ORDER BY l.valuenum) AS lactate_median
        FROM intervals i LEFT JOIN mimiciii.labevents l ON l.hadm_id=i.hadm_id AND l.charttime>=i.t0 AND l.charttime<i.t1
            AND l.itemid=50813 AND l.valuenum>0 GROUP BY 1,2
    ), urine AS (
        SELECT i.icustay_id,i.interval_index,COUNT(u.urineoutput) AS urine_observed_n,SUM(u.urineoutput) AS urine_output_ml
        FROM intervals i LEFT JOIN mimiciii.pivoted_uo u ON u.icustay_id=i.icustay_id AND u.charttime>=i.t0 AND u.charttime<i.t1
            AND u.urineoutput>=0 GROUP BY 1,2
    )
    SELECT i.icustay_id AS icu_stay_id,i.trajectory_class,i.interval_index,i.weight_raw_kg,i.weight_kg,
        1 AS at_risk_start,6 AS observed_window_hours,n.nee_max,n.nee_any,n.nee_unquantified,
        f.input_events,f.input_ml,f.audited_iv_ml,
        CASE WHEN f.audited_iv_ml=0 THEN 0 ELSE f.audited_iv_ml/i.weight_kg END AS audited_iv_ml_kg,
        f.mv_input_events,f.cv_input_events,o.output_events,o.output_ml,
        v.map_observed_n,v.map_median,l.lactate_observed_n,l.lactate_median,u.urine_observed_n,
        u.urine_output_ml,u.urine_output_ml/i.weight_kg/6.0 AS urine_output_ml_kg_h
    FROM intervals i LEFT JOIN nee_interval n USING(icustay_id,interval_index)
    LEFT JOIN input_interval f USING(icustay_id,interval_index) LEFT JOIN output_interval o USING(icustay_id,interval_index)
    LEFT JOIN vitals v USING(icustay_id,interval_index) LEFT JOIN lactate l USING(icustay_id,interval_index)
    LEFT JOIN urine u USING(icustay_id,interval_index) ORDER BY i.icustay_id,i.interval_index
