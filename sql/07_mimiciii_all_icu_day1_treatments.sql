-- MIMIC-III historical cross-version treatment extraction (ICU admission to 24 hours).
-- Read-only query; one row per eligible first ICU stay. Steroids are orders-only sensitivity data.

WITH icu_ranked AS (
    SELECT
        ie.subject_id, ie.hadm_id, ie.icustay_id, ie.intime, ie.outtime,
        adm.admittime, adm.dischtime, adm.has_chartevents_data,
        CASE WHEN DATE_PART('year', AGE(ie.intime, pat.dob))::integer > 89 THEN 90
             ELSE DATE_PART('year', AGE(ie.intime, pat.dob))::integer END AS age,
        ROW_NUMBER() OVER (PARTITION BY ie.subject_id ORDER BY ie.intime, ie.icustay_id) AS patient_icu_rank
    FROM mimiciii.icustays AS ie
    INNER JOIN mimiciii.admissions AS adm ON adm.hadm_id = ie.hadm_id
    INNER JOIN mimiciii.patients AS pat ON pat.subject_id = ie.subject_id
), eligible AS MATERIALIZED (
    SELECT * FROM icu_ranked
    WHERE patient_icu_rank = 1 AND age >= 18 AND has_chartevents_data = 1
      AND outtime - intime >= INTERVAL '24 hours'
      AND dischtime - admittime >= INTERVAL '24 hours'
), weight AS (
    SELECT e.icustay_id, wd.weight AS weight_kg
    FROM eligible AS e
    LEFT JOIN LATERAL (
        SELECT w.weight FROM mimiciii.weight_durations AS w
        WHERE w.icustay_id = e.icustay_id AND w.starttime < e.intime + INTERVAL '24 hours' AND w.endtime > e.intime
        ORDER BY CASE WHEN w.starttime <= e.intime AND w.endtime > e.intime THEN 0 ELSE 1 END, w.starttime
        LIMIT 1
    ) AS wd ON TRUE
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
), input_union AS MATERIALIZED (
    SELECT im.icustay_id, im.starttime, im.endtime, im.itemid, im.amount, im.amountuom
    FROM mimiciii.inputevents_mv AS im
    INNER JOIN eligible AS e USING(icustay_id)
    WHERE im.cancelreason = 0
      AND im.starttime < e.intime + INTERVAL '24 hours'
      AND COALESCE(im.endtime, im.starttime) >= e.intime
      AND im.itemid IN (30008,30009,30018,30020,30021,30143,30159,30160,30161,30176,30181,30352,30353,
                        40548,42832,43237,43353,44203,44440,44952,45403,45532,46564)
    UNION ALL
    SELECT cv.icustay_id, cv.charttime, cv.charttime, cv.itemid, cv.amount, cv.amountuom
    FROM mimiciii.inputevents_cv AS cv
    INNER JOIN eligible AS e USING(icustay_id)
    WHERE COALESCE(cv.stopped, '') NOT ILIKE '%stop%'
      AND cv.charttime >= e.intime AND cv.charttime < e.intime + INTERVAL '24 hours'
      AND cv.itemid IN (30008,30009,30018,30020,30021,30143,30159,30160,30161,30176,30181,30352,30353,
                        40548,42832,43237,43353,44203,44440,44952,45403,45532,46564)
), fluid_input AS (
    SELECT
        e.icustay_id,
        COUNT(iu.*) FILTER (WHERE LOWER(iu.amountuom) = 'ml' AND iu.amount > 0) AS n_input_events,
        SUM(CASE
            WHEN LOWER(iu.amountuom) <> 'ml' OR iu.amount <= 0 THEN NULL
            WHEN iu.endtime > iu.starttime THEN iu.amount *
                EXTRACT(EPOCH FROM (LEAST(iu.endtime, e.intime + INTERVAL '24 hours') - GREATEST(iu.starttime, e.intime)))
                / NULLIF(EXTRACT(EPOCH FROM (iu.endtime - iu.starttime)), 0)
            WHEN iu.starttime >= e.intime AND iu.starttime < e.intime + INTERVAL '24 hours' THEN iu.amount
        END) AS fluid_input_ml_24h,
        SUM(CASE WHEN iu.itemid IN (30008,30009,30018,30020,30021,30143,30159,30160,30161,30176,30181,30352,30353,
                                    40548,42832,43237,43353,44203,44440,44952,45403,45532,46564)
             AND LOWER(iu.amountuom) = 'ml' AND iu.amount > 0 THEN
             CASE WHEN iu.endtime > iu.starttime THEN iu.amount *
               EXTRACT(EPOCH FROM (LEAST(iu.endtime,e.intime+INTERVAL '24 hours')-GREATEST(iu.starttime,e.intime))) /
               NULLIF(EXTRACT(EPOCH FROM(iu.endtime-iu.starttime)),0) ELSE iu.amount END END) AS resuscitation_fluid_ml_24h,
        SUM(CASE WHEN iu.itemid IN (30018,30020,30021,30143,30159,30160,30161,30176,30352,30353,44440,45532)
             AND LOWER(iu.amountuom) = 'ml' AND iu.amount > 0 THEN
             CASE WHEN iu.endtime > iu.starttime THEN iu.amount *
               EXTRACT(EPOCH FROM (LEAST(iu.endtime,e.intime+INTERVAL '24 hours')-GREATEST(iu.starttime,e.intime))) /
               NULLIF(EXTRACT(EPOCH FROM(iu.endtime-iu.starttime)),0) ELSE iu.amount END END) AS crystalloid_fluid_ml_24h,
        SUM(CASE WHEN iu.itemid IN (30008,30009,30181,40548,42832,43237,43353,44203,44952,45403,46564)
             AND LOWER(iu.amountuom) = 'ml' AND iu.amount > 0 THEN
             CASE WHEN iu.endtime > iu.starttime THEN iu.amount *
               EXTRACT(EPOCH FROM (LEAST(iu.endtime,e.intime+INTERVAL '24 hours')-GREATEST(iu.starttime,e.intime))) /
               NULLIF(EXTRACT(EPOCH FROM(iu.endtime-iu.starttime)),0) ELSE iu.amount END END) AS colloid_fluid_ml_24h
    FROM eligible AS e
    LEFT JOIN input_union AS iu
        ON iu.icustay_id = e.icustay_id
       AND iu.starttime < e.intime + INTERVAL '24 hours'
       AND iu.endtime >= e.intime
    GROUP BY e.icustay_id
), fluid_output AS MATERIALIZED (
    SELECT
        e.icustay_id,
        COUNT(oe.*) FILTER (WHERE LOWER(oe.valueuom) = 'ml' AND oe.value > 0 AND COALESCE(oe.iserror, 0) = 0) AS n_output_events,
        SUM(oe.value) FILTER (WHERE LOWER(oe.valueuom) = 'ml' AND oe.value > 0 AND COALESCE(oe.iserror, 0) = 0) AS fluid_output_ml_24h
    FROM eligible AS e
    LEFT JOIN mimiciii.outputevents AS oe
        ON oe.icustay_id = e.icustay_id AND oe.charttime >= e.intime AND oe.charttime < e.intime + INTERVAL '24 hours'
       AND LOWER(oe.valueuom) = 'ml' AND oe.value > 0 AND COALESCE(oe.iserror, 0) = 0
    GROUP BY e.icustay_id
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
), nee_boundaries AS (
    SELECT icustay_id, intime AS boundary FROM eligible
    UNION SELECT icustay_id, intime + INTERVAL '24 hours' FROM eligible
    UNION SELECT ce.icustay_id, GREATEST(ce.starttime, e.intime) FROM component_events ce JOIN eligible e USING(icustay_id)
          WHERE ce.starttime < e.intime + INTERVAL '24 hours' AND ce.endtime > e.intime
    UNION SELECT ce.icustay_id, LEAST(ce.endtime, e.intime + INTERVAL '24 hours') FROM component_events ce JOIN eligible e USING(icustay_id)
          WHERE ce.starttime < e.intime + INTERVAL '24 hours' AND ce.endtime > e.intime
), nee_segments AS (
    SELECT icustay_id, boundary AS segment_start,
           LEAD(boundary) OVER (PARTITION BY icustay_id ORDER BY boundary) AS segment_end
    FROM nee_boundaries
), per_drug AS (
    SELECT s.icustay_id,s.segment_start,s.segment_end,ce.drug,MAX(ce.rate*ce.factor) AS component_nee
    FROM nee_segments s
    LEFT JOIN component_events ce ON ce.icustay_id=s.icustay_id AND ce.starttime<=s.segment_start AND ce.endtime>s.segment_start
    WHERE s.segment_end>s.segment_start
    GROUP BY s.icustay_id,s.segment_start,s.segment_end,ce.drug
), nee_rates AS (
    SELECT icustay_id,segment_start,segment_end,COALESCE(SUM(component_nee),0) AS nee
    FROM per_drug GROUP BY icustay_id,segment_start,segment_end
), nee_summary AS (
    SELECT icustay_id,(MAX(nee)>0)::integer AS nee_any_24h,MAX(nee) AS nee_max_24h,
           SUM(nee*EXTRACT(EPOCH FROM(segment_end-segment_start))/3600.0)/24.0 AS nee_twa_24h,
           SUM(nee*EXTRACT(EPOCH FROM(segment_end-segment_start))/3600.0) AS nee_auc_24h,
           SUM(EXTRACT(EPOCH FROM(segment_end-segment_start))/3600.0) FILTER(WHERE nee>0) AS nee_duration_hours_24h
    FROM nee_rates GROUP BY icustay_id
), vaso_components AS (
    SELECT e.icustay_id,
           MAX(ce.rate) FILTER(WHERE ce.drug='norepinephrine') AS norepinephrine_max,
           MAX(ce.rate) FILTER(WHERE ce.drug='epinephrine') AS epinephrine_max,
           MAX(ce.rate) FILTER(WHERE ce.drug='phenylephrine') AS phenylephrine_max,
           MAX(ce.rate) FILTER(WHERE ce.drug='dopamine') AS dopamine_max,
           MAX(ce.rate) FILTER(WHERE ce.drug='vasopressin') AS vasopressin_u_min_max
    FROM eligible e LEFT JOIN component_events ce ON ce.icustay_id=e.icustay_id
      AND ce.starttime<e.intime+INTERVAL '24 hours' AND ce.endtime>e.intime
    GROUP BY e.icustay_id
), steroid_rows AS MATERIALIZED (
    SELECT e.icustay_id,LOWER(p.drug) AS medication,LOWER(COALESCE(p.route,'')) AS route,
           CASE WHEN COALESCE(p.dose_val_rx,'') ~ '[0-9]' THEN SUBSTRING(p.dose_val_rx FROM '[0-9]+[.]?[0-9]*')::numeric END
           * CASE WHEN LOWER(p.dose_unit_rx) IN ('mcg','ug') THEN 0.001 WHEN LOWER(p.dose_unit_rx)='g' THEN 1000.0
                  WHEN LOWER(p.dose_unit_rx)='mg' THEN 1.0 END AS dose_mg
    FROM eligible e INNER JOIN mimiciii.prescriptions p ON p.hadm_id=e.hadm_id
      AND p.startdate < e.intime + INTERVAL '24 hours' AND COALESCE(p.enddate,p.startdate)>=e.intime
    WHERE LOWER(p.drug) ~ 'hydrocortisone|methylprednisolone|prednisolone|prednisone|dexamethasone'
      AND LOWER(COALESCE(p.route,'')) ~ '(^|[^a-z])(iv|intravenous|po|oral|ng|gt|im)([^a-z]|$)'
      AND LOWER(p.drug) !~ 'ophth|eye|drop|topical|cream|ointment|inhal|nasal|rectal'
), steroid_summary AS (
    SELECT icustay_id,(COUNT(*) FILTER(WHERE dose_mg>0)>0)::integer AS systemic_steroid_any_24h,
      SUM(dose_mg) FILTER(WHERE medication~'hydrocortisone') AS hydrocortisone_mg_24h,
      SUM(dose_mg) FILTER(WHERE medication~'methylprednisolone') AS methylprednisolone_mg_24h,
      SUM(dose_mg) FILTER(WHERE medication~'prednisolone' AND medication!~'methyl') AS prednisolone_mg_24h,
      SUM(dose_mg) FILTER(WHERE medication~'prednisone' AND medication!~'methyl') AS prednisone_mg_24h,
      SUM(dose_mg) FILTER(WHERE medication~'dexamethasone') AS dexamethasone_mg_24h,
      COUNT(*) FILTER(WHERE dose_mg IS NULL) AS steroid_unparsed_rows
    FROM steroid_rows GROUP BY icustay_id
)
SELECT 'mimiciii'::text AS dataset,e.subject_id::text AS subject_id,e.hadm_id::text AS encounter_id,
  e.icustay_id::text AS icu_stay_id,e.intime AS exposure_start,e.intime+INTERVAL '24 hours' AS exposure_end,
  w.weight_kg,fi.n_input_events,fo.n_output_events,(fi.n_input_events>0)::integer AS input_documented,
  (fo.n_output_events>0)::integer AS output_documented,fi.fluid_input_ml_24h,fo.fluid_output_ml_24h,
  CASE WHEN fi.fluid_input_ml_24h IS NOT NULL AND fo.fluid_output_ml_24h IS NOT NULL THEN fi.fluid_input_ml_24h-fo.fluid_output_ml_24h END AS fluid_balance_ml_24h,
  CASE WHEN fi.n_input_events>0 THEN COALESCE(fi.resuscitation_fluid_ml_24h,0) END AS resuscitation_fluid_ml_24h,
  CASE WHEN fi.n_input_events>0 THEN COALESCE(fi.crystalloid_fluid_ml_24h,0) END AS crystalloid_fluid_ml_24h,
  CASE WHEN fi.n_input_events>0 THEN COALESCE(fi.colloid_fluid_ml_24h,0) END AS colloid_fluid_ml_24h,
  fi.fluid_input_ml_24h/NULLIF(w.weight_kg,0) AS fluid_input_ml_kg_24h,
  CASE WHEN fi.fluid_input_ml_24h IS NOT NULL AND fo.fluid_output_ml_24h IS NOT NULL THEN (fi.fluid_input_ml_24h-fo.fluid_output_ml_24h)/NULLIF(w.weight_kg,0) END AS fluid_balance_ml_kg_24h,
  COALESCE(ns.nee_any_24h,0) AS nee_any_24h,COALESCE(ns.nee_max_24h,0) AS nee_max_24h,
  COALESCE(ns.nee_twa_24h,0) AS nee_twa_24h,COALESCE(ns.nee_auc_24h,0) AS nee_auc_24h,
  COALESCE(ns.nee_duration_hours_24h,0) AS nee_duration_hours_24h,
  vc.norepinephrine_max,vc.epinephrine_max,vc.phenylephrine_max,vc.dopamine_max,
  0::numeric AS metaraminol_max,vc.vasopressin_u_min_max,0::numeric AS angiotensin_ii_max,
  1::integer AS metaraminol_source_unavailable,1::integer AS angiotensin_ii_source_unavailable,
  COALESCE(ss.systemic_steroid_any_24h,0) AS systemic_steroid_any_24h,
  COALESCE(ss.hydrocortisone_mg_24h,0) AS hydrocortisone_mg_24h,
  COALESCE(ss.methylprednisolone_mg_24h,0) AS methylprednisolone_mg_24h,
  COALESCE(ss.prednisolone_mg_24h,0) AS prednisolone_mg_24h,
  COALESCE(ss.prednisone_mg_24h,0) AS prednisone_mg_24h,
  COALESCE(ss.dexamethasone_mg_24h,0) AS dexamethasone_mg_24h,
  COALESCE(ss.hydrocortisone_mg_24h,0)+COALESCE(ss.methylprednisolone_mg_24h,0)*5.0+COALESCE(ss.prednisolone_mg_24h,0)*4.0+COALESCE(ss.dexamethasone_mg_24h,0)*26.67 AS hce_mg_24h,
  COALESCE(ss.hydrocortisone_mg_24h,0)+COALESCE(ss.methylprednisolone_mg_24h,0)*5.0+(COALESCE(ss.prednisolone_mg_24h,0)+COALESCE(ss.prednisone_mg_24h,0))*4.0+COALESCE(ss.dexamethasone_mg_24h,0)*26.67 AS hce_with_prednisone_mg_24h,
  COALESCE(ss.steroid_unparsed_rows,0) AS steroid_unparsed_rows,'orders_only_prescriptions'::text AS steroid_source,
  'NE+E+phenylephrine/10+dopamine/100+metaraminol/8+vasopressin(U/min)*2.5+angiotensinII*10'::text AS nee_formula,
  'hydrocortisone+methylprednisolone*5+prednisolone*4+dexamethasone*26.67'::text AS hce_formula
FROM eligible e LEFT JOIN weight w USING(icustay_id) LEFT JOIN fluid_input fi USING(icustay_id)
LEFT JOIN fluid_output fo USING(icustay_id) LEFT JOIN nee_summary ns USING(icustay_id)
LEFT JOIN vaso_components vc USING(icustay_id) LEFT JOIN steroid_summary ss USING(icustay_id);
