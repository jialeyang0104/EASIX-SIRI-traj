BEGIN READ ONLY;
COPY (-- MIMIC-IV all-ICU D1-D7 daily EASIX/SIRI dataset.
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
) SELECT stay_id AS icu_stay_id, subject_id, hadm_id AS encounter_id,
        intime,outtime,admittime,dischtime,deathtime,dod,analysis_deathtime AS corrected_death FROM eligible) TO STDOUT WITH CSV HEADER;
COMMIT;
