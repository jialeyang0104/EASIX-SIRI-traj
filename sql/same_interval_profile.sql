/*
   Same-interval descriptive profile extraction for the frozen MIMIC-IV K=3
   trajectory cohort.  The caller must create a temporary table named
   tmp_trajectory_assignment (icu_stay_id bigint, trajectory_class integer).

   The query deliberately uses the actual event timestamp relative to ICU
   intime.  Intervals are [0,6), [6,12), ... [162,168) hours; an event exactly
   on a boundary is assigned to the interval on its right.  Rows are retained
   for all assigned stays and intervals; support states remain NULL when the
   source has no observation rather than being recoded to absence.  The
   audited-IV action flag is complete for at-risk rows; no audited event is a
   structural zero for the aggregate exposure, whereas an observed event with
   no valid weight remains missing only for mL/kg.
*/

WITH ranked AS (
    SELECT
        ie.subject_id, ie.hadm_id, ie.stay_id, ie.intime, ie.outtime,
        adm.admittime, adm.dischtime, adm.deathtime, adm.race,
        pat.gender, pat.dod,
        LEAST(
            pat.anchor_age + DATE_PART(
                'year', AGE(adm.admittime,
                    MAKE_TIMESTAMP(pat.anchor_year, 1, 1, 0, 0, 0))
            )::integer,
            90
        ) AS age,
        ROW_NUMBER() OVER (
            PARTITION BY ie.subject_id ORDER BY ie.intime, ie.stay_id
        ) AS patient_icu_rank
    FROM mimiciv_icu.icustays AS ie
    INNER JOIN mimiciv_hosp.admissions AS adm ON adm.hadm_id = ie.hadm_id
    INNER JOIN mimiciv_hosp.patients AS pat ON pat.subject_id = ie.subject_id
), eligible AS (
    SELECT
        r.*,
        CASE
            WHEN r.deathtime IS NOT NULL THEN r.deathtime
            WHEN r.dod IS NOT NULL THEN GREATEST(r.dod::timestamp, r.dischtime)
        END AS analysis_death_ts
    FROM ranked AS r
    WHERE r.patient_icu_rank = 1
      AND r.age >= 18
      AND r.outtime - r.intime >= INTERVAL '24 hours'
      AND r.dischtime - r.admittime >= INTERVAL '24 hours'
), cohort AS (
    SELECT
        e.*, a.trajectory_class,
        COALESCE(w.weight_admit, w.weight)::double precision AS weight_kg
    FROM eligible AS e
    INNER JOIN tmp_trajectory_assignment AS a ON a.icu_stay_id = e.stay_id
    LEFT JOIN mimiciv_derived.first_day_weight AS w ON w.stay_id = e.stay_id
), intervals AS (
    SELECT
        c.stay_id, c.subject_id, c.hadm_id, c.intime, c.outtime, c.dischtime,
        c.age, c.gender, c.race, c.analysis_death_ts, c.weight_kg,
        c.trajectory_class, g.interval_index,
        c.intime + (g.interval_index - 1) * INTERVAL '6 hours' AS interval_start,
        c.intime + g.interval_index * INTERVAL '6 hours' AS interval_end,
        GREATEST(
            EXTRACT(EPOCH FROM (
                LEAST(
                    c.intime + g.interval_index * INTERVAL '6 hours',
                    c.outtime,
                    COALESCE(c.analysis_death_ts,
                        c.intime + g.interval_index * INTERVAL '6 hours')
                ) - (c.intime + (g.interval_index - 1) * INTERVAL '6 hours')
            )) / 3600.0,
            0.0
        ) AS observed_window_hours,
        CASE
            WHEN c.outtime > c.intime + (g.interval_index - 1) * INTERVAL '6 hours'
             AND (c.analysis_death_ts IS NULL OR c.analysis_death_ts >
                    c.intime + (g.interval_index - 1) * INTERVAL '6 hours')
            THEN 1 ELSE 0
        END AS at_risk_start
    FROM cohort AS c
    CROSS JOIN LATERAL generate_series(1, 4) AS g(interval_index)
), interval_bounds AS (
    SELECT
        i.*,
        LEAST(i.interval_end, i.outtime,
              COALESCE(i.analysis_death_ts, i.interval_end)) AS effective_end,
        CASE WHEN i.observed_window_hours >= 6.0 THEN 1 ELSE 0 END AS complete_interval,
        CASE WHEN i.observed_window_hours > 0 AND i.observed_window_hours < 6.0
             THEN 1 ELSE 0 END AS partial_interval
    FROM intervals AS i
),
/* NEE is the frozen project's derived norepinephrine-equivalent dose. */
nee_events AS (
    SELECT
        i.stay_id, i.interval_index, i.interval_start, i.interval_end,
        n.starttime, n.endtime,
        n.norepinephrine_equivalent_dose::double precision AS nee
    FROM interval_bounds AS i
    INNER JOIN mimiciv_derived.norepinephrine_equivalent_dose AS n
        ON n.stay_id = i.stay_id
       AND n.starttime < i.effective_end
       AND n.endtime > i.interval_start
       AND n.norepinephrine_equivalent_dose > 0
    WHERE i.at_risk_start = 1
), nee_boundaries AS (
    SELECT stay_id, interval_index, interval_start AS boundary FROM interval_bounds
    WHERE at_risk_start = 1
    UNION
    SELECT stay_id, interval_index, interval_end FROM interval_bounds
    WHERE at_risk_start = 1
    UNION
    SELECT stay_id, interval_index, GREATEST(starttime, interval_start)
    FROM nee_events
    UNION
    SELECT stay_id, interval_index, LEAST(endtime, interval_end)
    FROM nee_events
), nee_segments AS (
    SELECT
        stay_id, interval_index, boundary AS segment_start,
        LEAD(boundary) OVER (PARTITION BY stay_id, interval_index ORDER BY boundary)
            AS segment_end
    FROM nee_boundaries
), nee_segment_rates AS (
    SELECT
        s.stay_id, s.interval_index, s.segment_start, s.segment_end,
        COALESCE(MAX(n.nee), 0.0) AS nee
    FROM nee_segments AS s
    LEFT JOIN nee_events AS n
      ON n.stay_id = s.stay_id
     AND n.interval_index = s.interval_index
     AND n.starttime <= s.segment_start
     AND n.endtime > s.segment_start
    WHERE s.segment_end > s.segment_start
    GROUP BY s.stay_id, s.interval_index, s.segment_start, s.segment_end
), nee_profile AS (
    SELECT
        i.stay_id, i.interval_index,
        CASE WHEN i.at_risk_start = 1 THEN COALESCE(MAX(r.nee), 0.0) END AS nee_max,
        CASE WHEN i.at_risk_start = 1 THEN COALESCE(SUM(
            r.nee * EXTRACT(EPOCH FROM (r.segment_end-r.segment_start))/3600.0
        ), 0.0) END AS nee_auc,
        CASE WHEN i.at_risk_start = 1 THEN COALESCE(SUM(
            EXTRACT(EPOCH FROM (r.segment_end-r.segment_start))/3600.0
        ) FILTER (WHERE r.nee > 0), 0.0) END AS nee_positive_duration_hours
    FROM interval_bounds AS i
    LEFT JOIN nee_segment_rates AS r
      ON r.stay_id = i.stay_id AND r.interval_index = i.interval_index
    GROUP BY i.stay_id, i.interval_index, i.at_risk_start
),
/* Audited high-certainty IV-fluid mapping exactly reproduces the frozen
   longitudinal contract; low-certainty enteral, nutrition, carrier, blood and
   unclassified input is not silently promoted to a treatment action. */
fluid_events AS (
    SELECT
        i.stay_id, i.interval_index,
        ie.starttime, COALESCE(ie.endtime, ie.starttime) AS endtime,
        ie.amount::double precision AS amount,
        CASE
            WHEN (
                CASE
                    WHEN LOWER(COALESCE(ie.ordercategoryname,'')) LIKE '%oral%'
                      OR LOWER(COALESCE(ie.ordercategoryname,'')) LIKE '%gastric%'
                      OR LOWER(COALESCE(di.label,'')) ~ '(^|[^a-z])(po|gt|gastric|free water)([^a-z]|$)'
                        THEN 'enteral_oral'
                    WHEN LOWER(COALESCE(di.category,'')) LIKE '%nutrition%'
                      OR LOWER(COALESCE(di.label,'')) ~ 'tube feed|replete|nepro|osmolite|jevity|promote|glucerna|nutren'
                        THEN 'nutrition'
                    WHEN LOWER(COALESCE(di.category,'')) LIKE '%blood%'
                      AND LOWER(COALESCE(di.label,'')) !~ 'albumin|hetastarch|dextran'
                        THEN 'blood_product'
                    WHEN LOWER(COALESCE(ie.ordercategoryname,'')) LIKE '%antibiotic%'
                      OR LOWER(COALESCE(ie.ordercategoryname,'')) LIKE '%prophylaxis%'
                        THEN 'medication_carrier'
                    WHEN LOWER(COALESCE(di.label,'')) ~ 'albumin|hetastarch|dextran|colloid'
                        THEN 'colloid'
                    WHEN LOWER(COALESCE(di.category,''))='fluids/intake'
                      AND (LOWER(COALESCE(ie.ordercategoryname,'')) LIKE '02-%fluids (crystalloids)%'
                        OR LOWER(COALESCE(ie.ordercategoryname,'')) LIKE '03-%iv fluid bolus%')
                      AND LOWER(COALESCE(di.label,'')) NOT IN ('piggyback','solution','multivitamins')
                        THEN 'crystalloid'
                    WHEN LOWER(COALESCE(di.label,'')) ~ 'saline|nacl|ringer|(^|[^a-z])lr([^a-z]|$)|lactated|dextrose|d5|sterile water|crystalloid'
                        THEN 'crystalloid'
                    WHEN LOWER(COALESCE(ie.ordercategoryname,'')) LIKE '%drip%'
                        THEN 'maintenance_or_drip'
                    ELSE 'other_unclassified'
                END
            ) IN ('crystalloid','colloid')
            AND (
                LOWER(COALESCE(ie.ordercategoryname,'')) LIKE '02-%fluids (crystalloids)%'
                OR LOWER(COALESCE(ie.ordercategoryname,'')) LIKE '03-%iv fluid bolus%'
                OR LOWER(COALESCE(ie.ordercategoryname,'')) LIKE '04-%fluids (colloids)%'
            )
            AND LOWER(COALESCE(di.label,'')) NOT IN ('piggyback','solution','multivitamins')
            THEN 1 ELSE 0 END AS audited_item,
        i.interval_start, i.interval_end, i.effective_end
    FROM interval_bounds AS i
    INNER JOIN mimiciv_icu.inputevents AS ie
      ON ie.stay_id = i.stay_id
    INNER JOIN mimiciv_icu.d_items AS di ON di.itemid = ie.itemid
     AND ie.starttime < i.effective_end
     AND COALESCE(ie.endtime, ie.starttime) > i.interval_start
     AND ie.amount IS NOT NULL
     AND ie.amount > 0
     AND LOWER(COALESCE(ie.amountuom, '')) = 'ml'
     AND COALESCE(ie.statusdescription, '') NOT IN ('Rewritten', 'Cancelled')
    WHERE i.at_risk_start = 1
), fluid_profile AS (
    SELECT
        i.stay_id, i.interval_index,
        COUNT(f.*) FILTER (WHERE f.audited_item = 1) AS audited_event_n,
        COALESCE(SUM(
            CASE
                WHEN f.audited_item = 1 AND f.endtime > f.starttime
                THEN f.amount * EXTRACT(EPOCH FROM (
                    LEAST(f.endtime, i.effective_end) -
                    GREATEST(f.starttime, i.interval_start)
                )) / NULLIF(EXTRACT(EPOCH FROM (f.endtime-f.starttime)), 0)
                WHEN f.audited_item = 1
                 AND f.starttime >= i.interval_start AND f.starttime < i.effective_end
                THEN f.amount
            END
        ), 0.0) AS audited_iv_ml,
        COUNT(f.*) AS all_ml_input_event_n
    FROM interval_bounds AS i
    LEFT JOIN fluid_events AS f
      ON f.stay_id = i.stay_id AND f.interval_index = i.interval_index
    GROUP BY i.stay_id, i.interval_index
),
vital_events AS (
    SELECT
        i.stay_id, i.interval_index, v.charttime,
        CASE
            WHEN v.mbp BETWEEN 20 AND 200 THEN v.mbp::double precision
            WHEN v.mbp_ni BETWEEN 20 AND 200 THEN v.mbp_ni::double precision
        END AS map_value,
        CASE WHEN v.heart_rate BETWEEN 20 AND 250
             THEN v.heart_rate::double precision END AS hr_value
    FROM interval_bounds AS i
    INNER JOIN mimiciv_derived.vitalsign AS v
      ON v.stay_id = i.stay_id
     AND v.charttime >= i.interval_start
     AND v.charttime < i.effective_end
    WHERE i.at_risk_start = 1
), vital_profile AS (
    SELECT
        i.stay_id, i.interval_index,
        COUNT(v.map_value) AS map_observed_n,
        percentile_cont(0.50) WITHIN GROUP (ORDER BY v.map_value) AS map_median,
        MIN(v.map_value) AS map_min,
        (array_agg(v.map_value ORDER BY v.charttime DESC)
            FILTER (WHERE v.map_value IS NOT NULL))[1] AS map_last,
        MAX(v.charttime) FILTER (WHERE v.map_value IS NOT NULL) AS map_last_time,
        COUNT(v.hr_value) AS hr_observed_n,
        percentile_cont(0.50) WITHIN GROUP (ORDER BY v.hr_value) AS hr_median,
        (array_agg(v.hr_value ORDER BY v.charttime DESC)
            FILTER (WHERE v.hr_value IS NOT NULL))[1] AS hr_last,
        MAX(v.charttime) FILTER (WHERE v.hr_value IS NOT NULL) AS hr_last_time
    FROM interval_bounds AS i
    LEFT JOIN vital_events AS v
      ON v.stay_id = i.stay_id AND v.interval_index = i.interval_index
    GROUP BY i.stay_id, i.interval_index
),
urine_events AS (
    SELECT i.stay_id, i.interval_index, u.charttime,
           u.urineoutput::double precision AS urine_output_ml
    FROM interval_bounds AS i
    INNER JOIN mimiciv_derived.urine_output AS u
      ON u.stay_id = i.stay_id
     AND u.charttime >= i.interval_start
     AND u.charttime < i.effective_end
     AND u.urineoutput IS NOT NULL
     AND u.urineoutput >= 0
    WHERE i.at_risk_start = 1
), urine_profile AS (
    SELECT
        i.stay_id, i.interval_index,
        COUNT(u.*) AS urine_observed_n,
        SUM(u.urine_output_ml) AS urine_output_ml,
        MAX(u.charttime) AS urine_last_time
    FROM interval_bounds AS i
    LEFT JOIN urine_events AS u
      ON u.stay_id = i.stay_id AND u.interval_index = i.interval_index
    GROUP BY i.stay_id, i.interval_index
),
lactate_events AS (
    SELECT
        i.stay_id, i.interval_index, le.charttime,
        le.valuenum::double precision AS lactate_value
    FROM interval_bounds AS i
    INNER JOIN mimiciv_hosp.labevents AS le
      ON le.hadm_id = i.hadm_id
     AND le.charttime >= i.interval_start
     AND le.charttime < i.effective_end
     AND le.itemid IN (50813, 52442, 53154)
     AND le.valuenum IS NOT NULL
     AND le.valuenum > 0
    WHERE i.at_risk_start = 1
), lactate_profile AS (
    SELECT
        i.stay_id, i.interval_index,
        COUNT(l.*) AS lactate_observed_n,
        percentile_cont(0.50) WITHIN GROUP (ORDER BY l.lactate_value) AS lactate_median,
        MAX(l.lactate_value) AS lactate_max,
        (array_agg(l.lactate_value ORDER BY l.charttime DESC)
            FILTER (WHERE l.lactate_value IS NOT NULL))[1] AS lactate_last,
        MAX(l.charttime) FILTER (WHERE l.lactate_value IS NOT NULL) AS lactate_last_time
    FROM interval_bounds AS i
    LEFT JOIN lactate_events AS l
      ON l.stay_id = i.stay_id AND l.interval_index = i.interval_index
    GROUP BY i.stay_id, i.interval_index
),
ventilation_profile AS (
    SELECT
        i.stay_id, i.interval_index,
        COUNT(v.*) AS ventilation_observed_n,
        MAX(CASE WHEN v.ventilation_status = 'InvasiveVent' THEN 1 ELSE 0 END)
            AS ventilation_active
    FROM interval_bounds AS i
    LEFT JOIN mimiciv_derived.ventilation AS v
      ON v.stay_id = i.stay_id
     AND v.starttime < i.effective_end
     AND v.endtime > i.interval_start
    WHERE i.at_risk_start = 1
    GROUP BY i.stay_id, i.interval_index
), crrt_profile AS (
    SELECT
        i.stay_id, i.interval_index,
        COUNT(c.*) AS crrt_observed_n,
        CASE WHEN COUNT(c.system_active) > 0 THEN MAX(c.system_active) END
            AS crrt_active
    FROM interval_bounds AS i
    LEFT JOIN mimiciv_derived.crrt AS c
      ON c.stay_id = i.stay_id
     AND c.charttime >= i.interval_start
     AND c.charttime < i.effective_end
    WHERE i.at_risk_start = 1
    GROUP BY i.stay_id, i.interval_index
), sofa_profile AS (
    SELECT
        i.stay_id, i.interval_index,
        COUNT(s.cardiovascular) AS sofa_cv_observed_n,
        percentile_cont(0.50) WITHIN GROUP (ORDER BY s.cardiovascular)
            FILTER (WHERE s.cardiovascular IS NOT NULL) AS sofa_cv_median,
        MAX(s.cardiovascular) FILTER (WHERE s.cardiovascular IS NOT NULL)
            AS sofa_cv_max,
        COUNT(s.respiration) AS sofa_resp_observed_n,
        percentile_cont(0.50) WITHIN GROUP (ORDER BY s.respiration)
            FILTER (WHERE s.respiration IS NOT NULL) AS sofa_resp_median,
        MAX(s.respiration) FILTER (WHERE s.respiration IS NOT NULL)
            AS sofa_resp_max
    FROM interval_bounds AS i
    LEFT JOIN mimiciv_derived.sofa AS s
      ON s.stay_id = i.stay_id
     AND s.starttime < i.effective_end
     AND s.endtime > i.interval_start
    WHERE i.at_risk_start = 1
    GROUP BY i.stay_id, i.interval_index
)
SELECT
    'mimiciv'::text AS dataset,
    i.stay_id::text AS icu_stay_id,
    i.subject_id::text AS subject_id,
    i.hadm_id::text AS encounter_id,
    i.trajectory_class,
    i.interval_index,
    i.intime AS icu_intime,
    i.interval_start,
    i.interval_end,
    i.effective_end,
    i.observed_window_hours,
    i.at_risk_start,
    i.complete_interval,
    i.partial_interval,
    CASE WHEN i.analysis_death_ts > i.interval_start
              AND i.analysis_death_ts <= i.interval_end THEN 1 ELSE 0 END
        AS death_in_interval,
    CASE WHEN i.outtime > i.interval_start
              AND i.outtime <= i.interval_end THEN 1 ELSE 0 END
        AS icu_discharge_in_interval,
    CASE WHEN i.dischtime > i.interval_start
              AND i.dischtime <= i.interval_end THEN 1 ELSE 0 END
        AS hospital_discharge_in_interval,
    i.age,
    i.gender::text AS sex,
    i.race::text AS race_ethnicity,
    i.weight_kg,
    n.nee_max,
    n.nee_auc,
    n.nee_positive_duration_hours,
    CASE WHEN i.at_risk_start = 1 AND n.nee_max > 0 THEN 1
         WHEN i.at_risk_start = 1 THEN 0 END AS nee_any,
    f.audited_event_n,
    CASE WHEN i.at_risk_start = 1 AND f.audited_event_n > 0 THEN f.audited_iv_ml
         WHEN i.at_risk_start = 1 THEN 0.0 END AS audited_iv_ml,
    CASE WHEN i.at_risk_start = 1 AND f.audited_event_n > 0
              AND i.weight_kg BETWEEN 20 AND 300
         THEN f.audited_iv_ml / NULLIF(i.weight_kg, 0)
         WHEN i.at_risk_start = 1 AND COALESCE(f.audited_event_n, 0) = 0 THEN 0.0 END AS audited_iv_ml_kg,
    CASE WHEN f.audited_event_n > 0 THEN 1
         WHEN i.at_risk_start = 1 THEN 0 END AS audited_iv_any,
    f.all_ml_input_event_n,
    v.map_observed_n, v.map_median, v.map_min, v.map_last, v.map_last_time,
    v.hr_observed_n, v.hr_median, v.hr_last, v.hr_last_time,
    u.urine_observed_n, u.urine_output_ml, u.urine_last_time,
    CASE WHEN u.urine_output_ml IS NOT NULL
              AND i.weight_kg BETWEEN 20 AND 300
              AND i.observed_window_hours > 0
         THEN u.urine_output_ml / i.weight_kg / i.observed_window_hours END
         AS urine_output_ml_kg_h,
    l.lactate_observed_n, l.lactate_median, l.lactate_max, l.lactate_last,
    l.lactate_last_time,
    so.sofa_cv_observed_n, so.sofa_cv_median, so.sofa_cv_max,
    so.sofa_resp_observed_n, so.sofa_resp_median, so.sofa_resp_max,
    ve.ventilation_observed_n,
    CASE WHEN ve.ventilation_observed_n > 0 THEN ve.ventilation_active END
        AS ventilation_active,
    cr.crrt_observed_n,
    cr.crrt_active
FROM interval_bounds AS i
LEFT JOIN nee_profile AS n USING (stay_id, interval_index)
LEFT JOIN fluid_profile AS f USING (stay_id, interval_index)
LEFT JOIN vital_profile AS v USING (stay_id, interval_index)
LEFT JOIN urine_profile AS u USING (stay_id, interval_index)
LEFT JOIN lactate_profile AS l USING (stay_id, interval_index)
LEFT JOIN ventilation_profile AS ve USING (stay_id, interval_index)
LEFT JOIN crrt_profile AS cr USING (stay_id, interval_index)
LEFT JOIN sofa_profile AS so USING (stay_id, interval_index)
ORDER BY i.stay_id, i.interval_index;
