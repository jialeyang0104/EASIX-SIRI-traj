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
), labs AS (
        SELECT e.stay_id, AVG(l.valuenum) FILTER (WHERE l.itemid=50861 AND l.valuenum >0) AS alt_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=50878 AND l.valuenum >0) AS ast_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=50863 AND l.valuenum >0) AS alp_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=50885 AND l.valuenum >=0) AS bilirubin_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=50862 AND l.valuenum >0) AS albumin_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=50912 AND l.valuenum >0) AS creatinine_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=51006 AND l.valuenum >0) AS bun_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=50954 AND l.valuenum >0) AS ldh_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=51301 AND l.valuenum >0) AS wbc_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=51265 AND l.valuenum >0) AS platelet_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=51222 AND l.valuenum >0) AS hemoglobin_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=51221 AND l.valuenum >0) AS hematocrit_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=51279 AND l.valuenum >0) AS rbc_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=51250 AND l.valuenum >0) AS mcv_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=51237 AND l.valuenum >0) AS inr_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=51274 AND l.valuenum >0) AS pt_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=51275 AND l.valuenum >0) AS ptt_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=50983 AND l.valuenum >0) AS sodium_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=50971 AND l.valuenum >0) AS potassium_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=50902 AND l.valuenum >0) AS chloride_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=50882 AND l.valuenum >0) AS bicarbonate_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=50868 AND l.valuenum >=0) AS aniongap_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=50893 AND l.valuenum >0) AS calcium_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=50960 AND l.valuenum >0) AS magnesium_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=50970 AND l.valuenum >0) AS phosphate_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=50931 AND l.valuenum >0) AS glucose_mean,
AVG(l.valuenum) FILTER (WHERE l.itemid=50813 AND l.valuenum >0) AS lactate_mean FROM eligible e JOIN mimiciv_hosp.labevents l
        ON l.hadm_id=e.hadm_id AND l.charttime>=e.intime AND l.charttime<e.intime+INTERVAL '24 hours'
        WHERE l.itemid IN (50861,50878,50863,50885,50862,50912,51006,50954,51301,51265,51222,51221,51279,51250,51237,51274,51275,50983,50971,50902,50882,50868,50893,50960,50970,50931,50813) GROUP BY e.stay_id
    ), vitals AS (
        SELECT e.stay_id, AVG(v.heart_rate) AS heart_rate_mean,
AVG(v.sbp) AS sbp_mean,
AVG(v.dbp) AS dbp_mean,
AVG(v.mbp) AS map_mean,
AVG(v.resp_rate) AS resp_rate_mean,
AVG(v.temperature) AS temperature_mean,
AVG(v.spo2) AS spo2_mean FROM eligible e JOIN mimiciv_derived.vitalsign v
        ON v.stay_id=e.stay_id AND v.charttime>=e.intime AND v.charttime<e.intime+INTERVAL '24 hours' GROUP BY e.stay_id
    ), neurologic AS (
        SELECT e.stay_id, AVG(g.gcs) AS gcs_mean FROM eligible e JOIN mimiciv_derived.gcs g
        ON g.stay_id=e.stay_id AND g.charttime>=e.intime AND g.charttime<e.intime+INTERVAL '24 hours' GROUP BY e.stay_id
    ), bloodgas AS (
        SELECT e.stay_id, AVG(b.ph) AS ph_mean,
AVG(b.po2) AS pao2_mean,
AVG(b.pco2) AS paco2_mean,
AVG(b.pao2fio2ratio) AS pf_ratio_mean,
AVG(b.baseexcess) AS baseexcess_mean FROM eligible e JOIN mimiciv_derived.bg b
        ON b.hadm_id=e.hadm_id
        AND b.charttime>=e.intime AND b.charttime<e.intime+INTERVAL '24 hours' AND b.specimen='ART.' GROUP BY e.stay_id
    ), urine AS (
        SELECT e.stay_id, SUM(u.urineoutput)/24.0 AS urine_hourly_mean FROM eligible e JOIN mimiciv_derived.urine_output u
        ON u.stay_id=e.stay_id AND u.charttime>=e.intime AND u.charttime<e.intime+INTERVAL '24 hours'
        WHERE u.urineoutput>=0 GROUP BY e.stay_id
    ) SELECT e.stay_id AS icu_stay_id,
        l.alt_mean,l.ast_mean,l.alp_mean,l.bilirubin_mean,l.albumin_mean,l.creatinine_mean,l.bun_mean,l.ldh_mean,l.wbc_mean,l.platelet_mean,l.hemoglobin_mean,l.hematocrit_mean,l.rbc_mean,l.mcv_mean,l.inr_mean,l.pt_mean,l.ptt_mean,l.sodium_mean,l.potassium_mean,l.chloride_mean,l.bicarbonate_mean,l.aniongap_mean,l.calcium_mean,l.magnesium_mean,l.phosphate_mean,l.glucose_mean,l.lactate_mean, v.heart_rate_mean,v.sbp_mean,v.dbp_mean,v.map_mean,v.resp_rate_mean,v.temperature_mean,v.spo2_mean,
        g.gcs_mean, b.ph_mean,b.pao2_mean,b.paco2_mean,b.pf_ratio_mean,b.baseexcess_mean, u.urine_hourly_mean
        FROM eligible e LEFT JOIN labs l USING (stay_id) LEFT JOIN vitals v USING (stay_id)
        LEFT JOIN neurologic g USING (stay_id) LEFT JOIN bloodgas b USING (stay_id) LEFT JOIN urine u USING (stay_id) ORDER BY e.stay_id) TO STDOUT WITH CSV HEADER;
COMMIT;
