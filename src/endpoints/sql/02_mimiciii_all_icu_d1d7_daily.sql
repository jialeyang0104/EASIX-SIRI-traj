-- MIMIC-III historical cross-version transport all-ICU D1-D7 daily EASIX/SIRI dataset.
-- Inclusion: adult, first ICU stay per patient, ICU LOS >=24h, hospital LOS >=24h.

WITH icu_ranked AS (
    SELECT
        ie.subject_id,
        ie.hadm_id,
        ie.icustay_id,
        ie.intime,
        ie.outtime,
        adm.admittime,
        adm.dischtime,
        adm.deathtime,
        adm.hospital_expire_flag,
        pat.dod,
        pat.gender,
        adm.ethnicity,
        CASE
            WHEN DATE_PART('year', AGE(ie.intime, pat.dob))::integer > 89 THEN 90
            ELSE DATE_PART('year', AGE(ie.intime, pat.dob))::integer
        END AS age,
        ROW_NUMBER() OVER (
            PARTITION BY ie.subject_id
            ORDER BY ie.intime, ie.icustay_id
        ) AS patient_icu_rank
    FROM mimiciii.icustays AS ie
    INNER JOIN mimiciii.admissions AS adm
        ON adm.hadm_id = ie.hadm_id
    INNER JOIN mimiciii.patients AS pat
        ON pat.subject_id = ie.subject_id
),
eligible AS (
    SELECT
        *,
        EXTRACT(EPOCH FROM (outtime - intime)) / 86400.0 AS icu_los_days,
        EXTRACT(EPOCH FROM (dischtime - admittime)) / 86400.0 AS hospital_los_days,
        CASE
            WHEN deathtime >= admittime THEN deathtime
            WHEN dod IS NOT NULL THEN GREATEST(dod::timestamp, dischtime)
        END AS analysis_deathtime
    FROM icu_ranked
    WHERE patient_icu_rank = 1
      AND age >= 18
      AND outtime - intime >= INTERVAL '24 hours'
      AND dischtime - admittime >= INTERVAL '24 hours'
),
sepsis_flag AS (
    SELECT DISTINCT s.icustay_id, 1::integer AS sepsis
    FROM mimiciii.sepsis3 AS s
    INNER JOIN eligible AS e
        ON e.icustay_id = s.icustay_id
    WHERE s.sepsis3 = TRUE
      AND COALESCE(s.sepsis_onset, s.sofa_time, s.suspected_infection_time) <= e.intime + INTERVAL '24 hours'
),
day1_sofa AS (
    SELECT e.icustay_id, MAX(s.sofa_24hours)::integer AS sofa_day1
    FROM eligible AS e
    INNER JOIN mimiciii.pivoted_sofa AS s
        ON s.icustay_id = e.icustay_id
       AND s.starttime < e.intime + INTERVAL '24 hours'
       AND s.endtime > e.intime
    GROUP BY e.icustay_id
),
day1_vasopressor AS (
    SELECT e.icustay_id, 1::integer AS vasopressor_day1
    FROM eligible AS e
    INNER JOIN mimiciii.vasopressor_durations AS v
        ON v.icustay_id = e.icustay_id
       AND v.starttime < e.intime + INTERVAL '24 hours'
       AND v.endtime > e.intime
    GROUP BY e.icustay_id
),
day1_crrt AS (
    SELECT e.icustay_id, 1::integer AS crrt_day1
    FROM eligible AS e
    INNER JOIN mimiciii.crrt_durations AS c
        ON c.icustay_id = e.icustay_id
       AND c.starttime < e.intime + INTERVAL '24 hours'
       AND c.endtime > e.intime
    GROUP BY e.icustay_id
),
day1_aki AS (
    SELECT e.icustay_id, MAX(k.aki_stage)::integer AS aki_stage_day1
    FROM eligible AS e
    INNER JOIN mimiciii.kdigo_stages AS k
        ON k.icustay_id = e.icustay_id
       AND k.charttime >= e.intime
       AND k.charttime < e.intime + INTERVAL '24 hours'
    GROUP BY e.icustay_id
),
daily_labs AS (
    SELECT
        e.icustay_id,
        FLOOR(EXTRACT(EPOCH FROM (le.charttime - e.intime)) / 86400)::integer + 1 AS day,
        MAX(le.valuenum) FILTER (WHERE le.itemid = 50954 AND le.valuenum > 0) AS ldh_u_l,
        MAX(le.valuenum) FILTER (WHERE le.itemid = 50912 AND le.valuenum > 0) AS creatinine_mg_dl,
        MIN(le.valuenum) FILTER (WHERE le.itemid = 51265 AND le.valuenum > 0) AS platelets_109_l,
        MAX(le.valuenum) FILTER (WHERE le.itemid = 51301 AND le.valuenum > 0) AS wbc_109_l,
        MAX(le.valuenum) FILTER (WHERE le.itemid = 51256 AND le.valuenum > 0) AS neutrophils_pct,
        MAX(le.valuenum) FILTER (WHERE le.itemid = 51254 AND le.valuenum > 0) AS monocytes_pct,
        MIN(le.valuenum) FILTER (WHERE le.itemid = 51244 AND le.valuenum > 0) AS lymphocytes_pct
    FROM eligible AS e
    INNER JOIN mimiciii.labevents AS le
        ON le.hadm_id = e.hadm_id
       AND le.charttime >= e.intime
       AND le.charttime < LEAST(e.outtime, e.intime + INTERVAL '7 days')
    WHERE le.itemid IN (50954, 50912, 51265, 51301, 51256, 51254, 51244)
    GROUP BY
        e.icustay_id,
        FLOOR(EXTRACT(EPOCH FROM (le.charttime - e.intime)) / 86400)::integer + 1
),
day_grid AS (
    SELECT e.*, gs.day
    FROM eligible AS e
    CROSS JOIN GENERATE_SERIES(1, 7) AS gs(day)
)
SELECT
    'mimiciii'::text AS dataset,
    g.subject_id::text AS subject_id,
    g.hadm_id::text AS encounter_id,
    g.icustay_id::text AS icu_stay_id,
    g.day,
    g.age,
    g.gender::text AS sex,
    g.ethnicity::text AS race_ethnicity,
    ROUND(g.icu_los_days::numeric, 4) AS icu_los_days,
    ROUND(g.hospital_los_days::numeric, 4) AS hospital_los_days,
    COALESCE(sf.sepsis, 0) AS sepsis,
    s.sirs AS sirs_day1,
    sofa.sofa_day1,
    oasis.oasis AS oasis_day1,
    aps.apsiii AS apsiii_day1,
    COALESCE(vent.vent, 0) AS mechanical_ventilation_day1,
    COALESCE(dv.vasopressor_day1, 0) AS vasopressor_day1,
    COALESCE(dc.crrt_day1, 0) AS crrt_day1,
    COALESCE(da.aki_stage_day1, 0) AS aki_stage_day1,
    l.ldh_u_l,
    l.creatinine_mg_dl,
    l.platelets_109_l,
    l.wbc_109_l,
    l.neutrophils_pct,
    l.monocytes_pct,
    l.lymphocytes_pct,
    NULL::numeric AS neutrophils_abs,
    NULL::numeric AS monocytes_abs,
    NULL::numeric AS lymphocytes_abs,
    CASE
        WHEN l.ldh_u_l > 0 AND l.creatinine_mg_dl > 0 AND l.platelets_109_l > 0
            THEN l.ldh_u_l * l.creatinine_mg_dl / l.platelets_109_l
    END AS easix_day,
    CASE
        WHEN l.wbc_109_l > 0
         AND l.neutrophils_pct > 0
         AND l.monocytes_pct > 0
         AND l.lymphocytes_pct > 0
            THEN (l.wbc_109_l * l.neutrophils_pct / 100.0)
               * (l.wbc_109_l * l.monocytes_pct / 100.0)
               / NULLIF(l.wbc_109_l * l.lymphocytes_pct / 100.0, 0)
    END AS siri_day,
    CASE
        WHEN g.analysis_deathtime IS NOT NULL
            THEN GREATEST(EXTRACT(EPOCH FROM (g.analysis_deathtime - g.intime)) / 86400.0, 0.01)
    END AS death_time_days,
    CASE
        WHEN g.analysis_deathtime IS NOT NULL AND g.analysis_deathtime <= g.intime + INTERVAL '28 days'
            THEN GREATEST(EXTRACT(EPOCH FROM (g.analysis_deathtime - g.intime)) / 86400.0, 0.01)
        ELSE 28.0
    END AS followup_28d_days,
    CASE WHEN g.analysis_deathtime IS NOT NULL AND g.analysis_deathtime <= g.intime + INTERVAL '28 days' THEN 1 ELSE 0 END AS event_28d,
    CASE
        WHEN g.analysis_deathtime IS NOT NULL AND g.analysis_deathtime <= g.intime + INTERVAL '90 days'
            THEN GREATEST(EXTRACT(EPOCH FROM (g.analysis_deathtime - g.intime)) / 86400.0, 0.01)
        ELSE 90.0
    END AS followup_90d_days,
    CASE WHEN g.analysis_deathtime IS NOT NULL AND g.analysis_deathtime <= g.intime + INTERVAL '90 days' THEN 1 ELSE 0 END AS event_90d,
    CASE WHEN g.analysis_deathtime IS NOT NULL AND g.analysis_deathtime <= g.intime + INTERVAL '7 days' THEN 1 ELSE 0 END AS mortality_7d,
    CASE WHEN g.analysis_deathtime IS NOT NULL AND g.analysis_deathtime <= g.intime + INTERVAL '28 days' THEN 1 ELSE 0 END AS mortality_28d,
    CASE WHEN g.analysis_deathtime IS NOT NULL AND g.analysis_deathtime <= g.intime + INTERVAL '90 days' THEN 1 ELSE 0 END AS mortality_90d
FROM day_grid AS g
LEFT JOIN daily_labs AS l
    ON l.icustay_id = g.icustay_id
   AND l.day = g.day
LEFT JOIN sepsis_flag AS sf
    ON sf.icustay_id = g.icustay_id
LEFT JOIN mimiciii.sirs AS s
    ON s.icustay_id = g.icustay_id
LEFT JOIN day1_sofa AS sofa
    ON sofa.icustay_id = g.icustay_id
LEFT JOIN mimiciii.oasis AS oasis
    ON oasis.icustay_id = g.icustay_id
LEFT JOIN mimiciii.apsiii AS aps
    ON aps.icustay_id = g.icustay_id
LEFT JOIN mimiciii.ventilation_first_day AS vent
    ON vent.icustay_id = g.icustay_id
LEFT JOIN day1_vasopressor AS dv
    ON dv.icustay_id = g.icustay_id
LEFT JOIN day1_crrt AS dc
    ON dc.icustay_id = g.icustay_id
LEFT JOIN day1_aki AS da
    ON da.icustay_id = g.icustay_id
ORDER BY g.icustay_id, g.day;
