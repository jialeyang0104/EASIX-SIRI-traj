-- MIMIC-IV all-ICU D1-D7 daily EASIX/SIRI dataset.
-- Inclusion: adult, first ICU stay per patient, ICU LOS >=24h, hospital LOS >=24h.

WITH icu_ranked AS (
    SELECT
        ie.subject_id,
        ie.hadm_id,
        ie.stay_id,
        ie.intime,
        ie.outtime,
        adm.admittime,
        adm.dischtime,
        adm.deathtime,
        adm.hospital_expire_flag,
        pat.dod,
        pat.gender,
        adm.race,
        LEAST(
            pat.anchor_age
                + DATE_PART(
                    'year',
                    AGE(adm.admittime, MAKE_TIMESTAMP(pat.anchor_year, 1, 1, 0, 0, 0))
                )::integer,
            90
        ) AS age,
        ROW_NUMBER() OVER (
            PARTITION BY ie.subject_id
            ORDER BY ie.intime, ie.stay_id
        ) AS patient_icu_rank
    FROM mimiciv_icu.icustays AS ie
    INNER JOIN mimiciv_hosp.admissions AS adm
        ON adm.hadm_id = ie.hadm_id
    INNER JOIN mimiciv_hosp.patients AS pat
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
    SELECT DISTINCT s.stay_id, 1::integer AS sepsis
    FROM mimiciv_derived.sepsis3 AS s
    INNER JOIN eligible AS e
        ON e.stay_id = s.stay_id
    WHERE s.sepsis3 = TRUE
      AND COALESCE(s.sofa_time, s.suspected_infection_time) <= e.intime + INTERVAL '24 hours'
),
day1_support AS (
    SELECT
        e.stay_id,
        MAX(CASE WHEN v.ventilation_status = 'InvasiveVent' THEN 1 ELSE 0 END) AS mechanical_ventilation_day1
    FROM eligible AS e
    LEFT JOIN mimiciv_derived.ventilation AS v
        ON v.stay_id = e.stay_id
       AND v.starttime < e.intime + INTERVAL '24 hours'
       AND v.endtime > e.intime
    GROUP BY e.stay_id
),
day1_vasopressor AS (
    SELECT e.stay_id, 1::integer AS vasopressor_day1
    FROM eligible AS e
    INNER JOIN mimiciv_derived.vasoactive_agent AS v
        ON v.stay_id = e.stay_id
       AND v.starttime < e.intime + INTERVAL '24 hours'
       AND v.endtime > e.intime
    GROUP BY e.stay_id
),
day1_crrt AS (
    SELECT e.stay_id, 1::integer AS crrt_day1
    FROM eligible AS e
    INNER JOIN mimiciv_derived.crrt AS c
        ON c.stay_id = e.stay_id
       AND c.charttime >= e.intime
       AND c.charttime < e.intime + INTERVAL '24 hours'
    GROUP BY e.stay_id
),
day1_aki AS (
    SELECT e.stay_id, MAX(k.aki_stage)::integer AS aki_stage_day1
    FROM eligible AS e
    INNER JOIN mimiciv_derived.kdigo_stages AS k
        ON k.stay_id = e.stay_id
       AND k.charttime >= e.intime
       AND k.charttime < e.intime + INTERVAL '24 hours'
    GROUP BY e.stay_id
),
daily_labs AS (
    SELECT
        e.stay_id,
        FLOOR(EXTRACT(EPOCH FROM (le.charttime - e.intime)) / 86400)::integer + 1 AS day,
        MAX(le.valuenum) FILTER (WHERE le.itemid = 50954 AND le.valuenum > 0) AS ldh_u_l,
        MAX(le.valuenum) FILTER (WHERE le.itemid = 50912 AND le.valuenum > 0) AS creatinine_mg_dl,
        MIN(le.valuenum) FILTER (WHERE le.itemid = 51265 AND le.valuenum > 0) AS platelets_109_l,
        MAX(le.valuenum) FILTER (WHERE le.itemid = 51301 AND le.valuenum > 0) AS wbc_109_l,
        MAX(le.valuenum) FILTER (WHERE le.itemid = 51256 AND le.valuenum > 0) AS neutrophils_pct,
        MAX(le.valuenum) FILTER (WHERE le.itemid = 51254 AND le.valuenum > 0) AS monocytes_pct,
        MIN(le.valuenum) FILTER (WHERE le.itemid = 51244 AND le.valuenum > 0) AS lymphocytes_pct,
        MAX(le.valuenum) FILTER (WHERE le.itemid = 52075 AND le.valuenum > 0) AS neutrophils_abs,
        MAX(
            CASE
                WHEN le.itemid = 52074 AND le.valuenum > 0 THEN le.valuenum
                WHEN le.itemid = 51253 AND le.valuenum > 0 THEN le.valuenum / 1000.0
            END
        ) AS monocytes_abs,
        MIN(
            CASE
                WHEN le.itemid = 51133 AND le.valuenum > 0 THEN le.valuenum
                WHEN le.itemid = 52769 AND le.valuenum > 0 THEN le.valuenum / 1000.0
            END
        ) AS lymphocytes_abs
    FROM eligible AS e
    INNER JOIN mimiciv_hosp.labevents AS le
        ON le.hadm_id = e.hadm_id
       AND le.charttime >= e.intime
       AND le.charttime < LEAST(e.outtime, e.intime + INTERVAL '7 days')
    WHERE le.itemid IN (
        50954, 50912, 51265, 51301, 51256, 51254, 51244,
        52075, 52074, 51253, 51133, 52769
    )
    GROUP BY
        e.stay_id,
        FLOOR(EXTRACT(EPOCH FROM (le.charttime - e.intime)) / 86400)::integer + 1
),
day_grid AS (
    SELECT e.*, gs.day
    FROM eligible AS e
    CROSS JOIN GENERATE_SERIES(1, 7) AS gs(day)
),
assembled AS (
    SELECT
        g.*,
        COALESCE(sf.sepsis, 0) AS sepsis,
        s.sirs AS sirs_day1,
        sofa.sofa AS sofa_day1,
        oasis.oasis AS oasis_day1,
        aps.apsiii AS apsiii_day1,
        COALESCE(ds.mechanical_ventilation_day1, 0) AS mechanical_ventilation_day1,
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
        l.neutrophils_abs,
        l.monocytes_abs,
        l.lymphocytes_abs,
        CASE
            WHEN l.ldh_u_l > 0 AND l.creatinine_mg_dl > 0 AND l.platelets_109_l > 0
                THEN l.ldh_u_l * l.creatinine_mg_dl / l.platelets_109_l
        END AS easix_day,
        CASE
            WHEN l.neutrophils_abs > 0 AND l.monocytes_abs > 0 AND l.lymphocytes_abs > 0
                THEN l.neutrophils_abs * l.monocytes_abs / l.lymphocytes_abs
            WHEN l.wbc_109_l > 0
             AND l.neutrophils_pct > 0
             AND l.monocytes_pct > 0
             AND l.lymphocytes_pct > 0
                THEN (l.wbc_109_l * l.neutrophils_pct / 100.0)
                   * (l.wbc_109_l * l.monocytes_pct / 100.0)
                   / NULLIF(l.wbc_109_l * l.lymphocytes_pct / 100.0, 0)
        END AS siri_day
    FROM day_grid AS g
    LEFT JOIN daily_labs AS l
        ON l.stay_id = g.stay_id
       AND l.day = g.day
    LEFT JOIN sepsis_flag AS sf
        ON sf.stay_id = g.stay_id
    LEFT JOIN mimiciv_derived.sirs AS s
        ON s.stay_id = g.stay_id
    LEFT JOIN mimiciv_derived.first_day_sofa AS sofa
        ON sofa.stay_id = g.stay_id
    LEFT JOIN mimiciv_derived.oasis AS oasis
        ON oasis.stay_id = g.stay_id
    LEFT JOIN mimiciv_derived.apsiii AS aps
        ON aps.stay_id = g.stay_id
    LEFT JOIN day1_support AS ds
        ON ds.stay_id = g.stay_id
    LEFT JOIN day1_vasopressor AS dv
        ON dv.stay_id = g.stay_id
    LEFT JOIN day1_crrt AS dc
        ON dc.stay_id = g.stay_id
    LEFT JOIN day1_aki AS da
        ON da.stay_id = g.stay_id
)
SELECT
    'mimiciv'::text AS dataset,
    subject_id::text AS subject_id,
    hadm_id::text AS encounter_id,
    stay_id::text AS icu_stay_id,
    day,
    age,
    gender::text AS sex,
    race::text AS race_ethnicity,
    -- Preserve timestamp-derived precision for the D1–D7 denominator. The
    -- day-start contract compares these values (in days) with d−1; rounding
    -- to four decimals can move an ICU discharge a few seconds across a day
    -- boundary (for example, an actual LOS of 4.000046 days becoming 4.0000).
    icu_los_days,
    hospital_los_days,
    sepsis,
    sirs_day1,
    sofa_day1,
    oasis_day1,
    apsiii_day1,
    mechanical_ventilation_day1,
    vasopressor_day1,
    crrt_day1,
    aki_stage_day1,
    ldh_u_l,
    creatinine_mg_dl,
    platelets_109_l,
    wbc_109_l,
    neutrophils_pct,
    monocytes_pct,
    lymphocytes_pct,
    neutrophils_abs,
    monocytes_abs,
    lymphocytes_abs,
    easix_day,
    siri_day,
    CASE
        WHEN analysis_deathtime IS NOT NULL
            THEN GREATEST(EXTRACT(EPOCH FROM (analysis_deathtime - intime)) / 86400.0, 0.01)
    END AS death_time_days,
    CASE
        WHEN analysis_deathtime IS NOT NULL AND analysis_deathtime <= intime + INTERVAL '28 days'
            THEN GREATEST(EXTRACT(EPOCH FROM (analysis_deathtime - intime)) / 86400.0, 0.01)
        ELSE 28.0
    END AS followup_28d_days,
    CASE WHEN analysis_deathtime IS NOT NULL AND analysis_deathtime <= intime + INTERVAL '28 days' THEN 1 ELSE 0 END AS event_28d,
    CASE
        WHEN analysis_deathtime IS NOT NULL AND analysis_deathtime <= intime + INTERVAL '90 days'
            THEN GREATEST(EXTRACT(EPOCH FROM (analysis_deathtime - intime)) / 86400.0, 0.01)
        ELSE 90.0
    END AS followup_90d_days,
    CASE WHEN analysis_deathtime IS NOT NULL AND analysis_deathtime <= intime + INTERVAL '90 days' THEN 1 ELSE 0 END AS event_90d,
    CASE WHEN analysis_deathtime IS NOT NULL AND analysis_deathtime <= intime + INTERVAL '7 days' THEN 1 ELSE 0 END AS mortality_7d,
    CASE WHEN analysis_deathtime IS NOT NULL AND analysis_deathtime <= intime + INTERVAL '28 days' THEN 1 ELSE 0 END AS mortality_28d,
    CASE WHEN analysis_deathtime IS NOT NULL AND analysis_deathtime <= intime + INTERVAL '90 days' THEN 1 ELSE 0 END AS mortality_90d
FROM assembled
ORDER BY stay_id, day;
