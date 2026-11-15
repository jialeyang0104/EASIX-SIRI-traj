-- MIMIC-IV all-ICU treatment exposure extraction (ICU admission to 24 hours).
-- Read-only query; one row per eligible first ICU stay.
-- NEE is evaluated on concurrent rate intervals. Vasopressin is converted from U/h to U/min.

WITH icu_ranked AS (
    SELECT
        ie.subject_id, ie.hadm_id, ie.stay_id, ie.intime, ie.outtime,
        adm.admittime, adm.dischtime,
        LEAST(
            pat.anchor_age + DATE_PART(
                'year', AGE(adm.admittime, MAKE_TIMESTAMP(pat.anchor_year, 1, 1, 0, 0, 0))
            )::integer,
            90
        ) AS age,
        ROW_NUMBER() OVER (PARTITION BY ie.subject_id ORDER BY ie.intime, ie.stay_id) AS patient_icu_rank
    FROM mimiciv_icu.icustays AS ie
    INNER JOIN mimiciv_hosp.admissions AS adm ON adm.hadm_id = ie.hadm_id
    INNER JOIN mimiciv_hosp.patients AS pat ON pat.subject_id = ie.subject_id
), eligible AS (
    SELECT *
    FROM icu_ranked
    WHERE patient_icu_rank = 1
      AND age >= 18
      AND outtime - intime >= INTERVAL '24 hours'
      AND dischtime - admittime >= INTERVAL '24 hours'
), weight AS (
    SELECT e.stay_id, COALESCE(fd.weight_admit, fd.weight) AS weight_kg
    FROM eligible AS e
    LEFT JOIN mimiciv_derived.first_day_weight AS fd ON fd.stay_id = e.stay_id
), fluid_input AS (
    SELECT
        e.stay_id,
        COUNT(ie.*) FILTER (WHERE LOWER(ie.amountuom) = 'ml' AND ie.amount > 0) AS n_input_events,
        SUM(
            CASE
                WHEN LOWER(ie.amountuom) <> 'ml' OR ie.amount <= 0 THEN NULL
                WHEN ie.endtime > ie.starttime THEN ie.amount *
                    EXTRACT(EPOCH FROM (
                        LEAST(ie.endtime, e.intime + INTERVAL '24 hours')
                        - GREATEST(ie.starttime, e.intime)
                    )) / NULLIF(EXTRACT(EPOCH FROM (ie.endtime - ie.starttime)), 0)
                WHEN ie.starttime >= e.intime AND ie.starttime < e.intime + INTERVAL '24 hours' THEN ie.amount
            END
        ) AS fluid_input_ml_24h,
        SUM(
            CASE
                WHEN ie.itemid IN (220861, 220862, 220863, 220864, 220953, 220954, 220955, 220956, 220967, 220968,
                                   220979, 221001, 221003, 221195, 226364, 226365, 226375, 226376)
                 AND LOWER(ie.amountuom) = 'ml' AND ie.amount > 0
                THEN CASE
                    WHEN ie.endtime > ie.starttime THEN ie.amount *
                        EXTRACT(EPOCH FROM (
                            LEAST(ie.endtime, e.intime + INTERVAL '24 hours')
                            - GREATEST(ie.starttime, e.intime)
                        )) / NULLIF(EXTRACT(EPOCH FROM (ie.endtime - ie.starttime)), 0)
                    ELSE ie.amount
                END
            END
        ) AS resuscitation_fluid_ml_24h,
        SUM(
            CASE
                WHEN ie.itemid IN (220953, 220954, 220955, 220956, 220967, 220968, 226364, 226375)
                 AND LOWER(ie.amountuom) = 'ml' AND ie.amount > 0
                THEN CASE WHEN ie.endtime > ie.starttime THEN ie.amount *
                    EXTRACT(EPOCH FROM (LEAST(ie.endtime, e.intime + INTERVAL '24 hours') - GREATEST(ie.starttime, e.intime))) /
                    NULLIF(EXTRACT(EPOCH FROM (ie.endtime - ie.starttime)), 0) ELSE ie.amount END
            END
        ) AS crystalloid_fluid_ml_24h,
        SUM(
            CASE
                WHEN ie.itemid IN (220861, 220862, 220863, 220864, 220979, 221001, 221003, 221195, 226365, 226376)
                 AND LOWER(ie.amountuom) = 'ml' AND ie.amount > 0
                THEN CASE WHEN ie.endtime > ie.starttime THEN ie.amount *
                    EXTRACT(EPOCH FROM (LEAST(ie.endtime, e.intime + INTERVAL '24 hours') - GREATEST(ie.starttime, e.intime))) /
                    NULLIF(EXTRACT(EPOCH FROM (ie.endtime - ie.starttime)), 0) ELSE ie.amount END
            END
        ) AS colloid_fluid_ml_24h
    FROM eligible AS e
    LEFT JOIN mimiciv_icu.inputevents AS ie
        ON ie.stay_id = e.stay_id
       AND ie.starttime < e.intime + INTERVAL '24 hours'
       AND COALESCE(ie.endtime, ie.starttime) >= e.intime
       AND COALESCE(ie.statusdescription, '') NOT IN ('Rewritten', 'Cancelled')
    GROUP BY e.stay_id
), fluid_output AS (
    SELECT
        e.stay_id,
        COUNT(oe.*) FILTER (WHERE LOWER(oe.valueuom) = 'ml' AND oe.value > 0) AS n_output_events,
        SUM(oe.value) FILTER (WHERE LOWER(oe.valueuom) = 'ml' AND oe.value > 0) AS fluid_output_ml_24h
    FROM eligible AS e
    LEFT JOIN mimiciv_icu.outputevents AS oe
        ON oe.stay_id = e.stay_id
       AND oe.charttime >= e.intime
       AND oe.charttime < e.intime + INTERVAL '24 hours'
    GROUP BY e.stay_id
), component_events AS MATERIALIZED (
    SELECT
        ie.stay_id, ie.starttime, ie.endtime,
        CASE ie.itemid
            WHEN 221906 THEN 'norepinephrine'
            WHEN 221289 THEN 'epinephrine'
            WHEN 221749 THEN 'phenylephrine'
            WHEN 221662 THEN 'dopamine'
            WHEN 222315 THEN 'vasopressin'
            WHEN 229709 THEN 'angiotensin_ii'
            WHEN 229764 THEN 'angiotensin_ii'
        END::text AS drug,
        CASE
            WHEN ie.itemid = 222315 AND LOWER(ie.rateuom) IN ('units/hour','unit/hour','u/hour') THEN ie.rate / 60.0
            WHEN ie.itemid = 222315 AND LOWER(ie.rateuom) IN ('units/min','unit/min','u/min') THEN ie.rate
            WHEN ie.itemid IN (229709,229764) AND LOWER(ie.rateuom) = 'ng/kg/min' THEN ie.rate / 1000.0
            WHEN ie.itemid IN (229709,229764) AND LOWER(ie.rateuom) IN ('mcg/kg/min','ug/kg/min') THEN ie.rate
            WHEN LOWER(ie.rateuom) IN ('mcg/kg/min','ug/kg/min') THEN ie.rate
            WHEN LOWER(ie.rateuom) = 'mg/kg/min' THEN ie.rate * 1000.0
            WHEN LOWER(ie.rateuom) IN ('mcg/min','ug/min') THEN ie.rate / NULLIF(COALESCE(ie.patientweight,w.weight_kg),0)
        END::numeric AS rate,
        CASE ie.itemid WHEN 221749 THEN 0.1 WHEN 221662 THEN 0.01
             WHEN 222315 THEN 2.5 WHEN 229709 THEN 10.0 WHEN 229764 THEN 10.0 ELSE 1.0 END::numeric AS factor
    FROM mimiciv_icu.inputevents AS ie
    INNER JOIN eligible AS e ON e.stay_id = ie.stay_id
    LEFT JOIN weight AS w ON w.stay_id = ie.stay_id
    WHERE ie.itemid IN (221906,221289,221749,221662,222315,229709,229764)
      AND ie.rate > 0
      AND ie.starttime < e.intime + INTERVAL '24 hours'
      AND ie.endtime > e.intime
      AND COALESCE(ie.statusdescription,'') NOT IN ('Rewritten','Cancelled')
), nee_boundaries AS (
    SELECT stay_id, intime AS boundary FROM eligible
    UNION
    SELECT stay_id, intime + INTERVAL '24 hours' FROM eligible
    UNION
    SELECT n.stay_id, GREATEST(n.starttime, e.intime)
    FROM component_events AS n INNER JOIN eligible AS e ON e.stay_id = n.stay_id
    WHERE n.rate IS NOT NULL
    UNION
    SELECT n.stay_id, LEAST(n.endtime, e.intime + INTERVAL '24 hours')
    FROM component_events AS n INNER JOIN eligible AS e ON e.stay_id = n.stay_id
    WHERE n.rate IS NOT NULL
), nee_segments AS (
    SELECT stay_id, boundary AS segment_start,
           LEAD(boundary) OVER (PARTITION BY stay_id ORDER BY boundary) AS segment_end
    FROM nee_boundaries
), per_drug AS (
    SELECT s.stay_id,s.segment_start,s.segment_end,n.drug,MAX(n.rate*n.factor) AS component_nee
    FROM nee_segments AS s
    LEFT JOIN component_events AS n
      ON n.stay_id=s.stay_id AND n.starttime<=s.segment_start AND n.endtime>s.segment_start AND n.rate IS NOT NULL
    WHERE s.segment_end>s.segment_start
    GROUP BY s.stay_id,s.segment_start,s.segment_end,n.drug
), nee_rates AS (
    SELECT stay_id,segment_start,segment_end,COALESCE(SUM(component_nee),0) AS nee
    FROM per_drug GROUP BY stay_id,segment_start,segment_end
), nee_summary AS (
    SELECT
        stay_id,
        (MAX(nee) > 0)::integer AS nee_any_24h,
        MAX(nee) AS nee_max_24h,
        SUM(nee * EXTRACT(EPOCH FROM (segment_end - segment_start)) / 3600.0) / 24.0 AS nee_twa_24h,
        SUM(nee * EXTRACT(EPOCH FROM (segment_end - segment_start)) / 3600.0) AS nee_auc_24h,
        SUM(EXTRACT(EPOCH FROM (segment_end - segment_start)) / 3600.0) FILTER (WHERE nee > 0) AS nee_duration_hours_24h
    FROM nee_rates
    GROUP BY stay_id
), vaso_components AS (
    SELECT
        e.stay_id,
        MAX(v.rate) FILTER(WHERE v.drug='norepinephrine') AS norepinephrine_max,
        MAX(v.rate) FILTER(WHERE v.drug='epinephrine') AS epinephrine_max,
        MAX(v.rate) FILTER(WHERE v.drug='phenylephrine') AS phenylephrine_max,
        MAX(v.rate) FILTER(WHERE v.drug='dopamine') AS dopamine_max,
        MAX(v.rate) FILTER(WHERE v.drug='vasopressin') AS vasopressin_u_min_max,
        MAX(v.rate) FILTER(WHERE v.drug='angiotensin_ii') AS angiotensin_ii_max
    FROM eligible AS e
    LEFT JOIN component_events AS v ON v.stay_id=e.stay_id
    GROUP BY e.stay_id
), steroid_emar AS MATERIALIZED (
    SELECT
        e.stay_id, em.subject_id, em.emar_id, em.emar_seq,
        LOWER(em.medication) AS medication
    FROM eligible AS e
    INNER JOIN mimiciv_hosp.emar AS em
        ON em.subject_id = e.subject_id AND em.hadm_id = e.hadm_id
       AND em.charttime >= e.intime AND em.charttime < e.intime + INTERVAL '24 hours'
    WHERE LOWER(em.medication) ~ 'hydrocortisone|methylprednisolone|prednisolone|prednisone|dexamethasone'
      AND LOWER(em.medication) !~ 'ophth|eye|drop|topical|cream|ointment|inhal|nasal|rectal'
      AND LOWER(COALESCE(em.event_txt, '')) = 'administered'
), emar_detail_agg AS (
    SELECT
        se.stay_id, se.medication, se.subject_id, se.emar_id, se.emar_seq,
        MAX(NULLIF(TRIM(ed.dose_given), '')) AS dose_given,
        MAX(NULLIF(TRIM(ed.dose_given_unit), '')) AS dose_given_unit,
        MAX(NULLIF(TRIM(ed.route), '')) AS route
    FROM steroid_emar AS se
    INNER JOIN mimiciv_hosp.emar_detail AS ed
        ON ed.subject_id=se.subject_id AND ed.emar_id=se.emar_id AND ed.emar_seq=se.emar_seq
    GROUP BY se.stay_id, se.medication, se.subject_id, se.emar_id, se.emar_seq
), steroid_rows AS (
    SELECT
        ed.stay_id,
        ed.medication,
        LOWER(COALESCE(ed.route, '')) AS route,
        CASE
            WHEN COALESCE(ed.dose_given, '') ~ '^[[:space:]]*[0-9]+([.][0-9]+)?[[:space:]]*$'
            THEN TRIM(ed.dose_given)::numeric
        END * CASE
            WHEN LOWER(ed.dose_given_unit) IN ('mcg', 'ug') THEN 0.001
            WHEN LOWER(ed.dose_given_unit) = 'g' THEN 1000.0
            WHEN LOWER(ed.dose_given_unit) = 'mg' THEN 1.0
        END AS dose_mg
    FROM emar_detail_agg AS ed
    WHERE LOWER(COALESCE(ed.route, '')) ~ '(^|[^a-z])(iv|intravenous|po|oral|ng|gt|g[- ]?tube|im)([^a-z]|$)'
), steroid_summary AS (
    SELECT
        stay_id,
        (COUNT(*) FILTER (WHERE dose_mg > 0) > 0)::integer AS systemic_steroid_any_24h,
        SUM(dose_mg) FILTER (WHERE medication ~ 'hydrocortisone') AS hydrocortisone_mg_24h,
        SUM(dose_mg) FILTER (WHERE medication ~ 'methylprednisolone') AS methylprednisolone_mg_24h,
        SUM(dose_mg) FILTER (WHERE medication ~ 'prednisolone' AND medication !~ 'methyl') AS prednisolone_mg_24h,
        SUM(dose_mg) FILTER (WHERE medication ~ 'prednisone' AND medication !~ 'methyl') AS prednisone_mg_24h,
        SUM(dose_mg) FILTER (WHERE medication ~ 'dexamethasone') AS dexamethasone_mg_24h,
        COUNT(*) FILTER (WHERE dose_mg IS NULL) AS steroid_unparsed_rows
    FROM steroid_rows GROUP BY stay_id
)
SELECT
    'mimiciv'::text AS dataset,
    e.subject_id::text AS subject_id,
    e.hadm_id::text AS encounter_id,
    e.stay_id::text AS icu_stay_id,
    e.intime AS exposure_start,
    e.intime + INTERVAL '24 hours' AS exposure_end,
    w.weight_kg,
    fi.n_input_events,
    fo.n_output_events,
    (fi.n_input_events > 0)::integer AS input_documented,
    (fo.n_output_events > 0)::integer AS output_documented,
    fi.fluid_input_ml_24h,
    fo.fluid_output_ml_24h,
    CASE WHEN fi.fluid_input_ml_24h IS NOT NULL AND fo.fluid_output_ml_24h IS NOT NULL
         THEN fi.fluid_input_ml_24h - fo.fluid_output_ml_24h END AS fluid_balance_ml_24h,
    CASE WHEN fi.n_input_events > 0 THEN COALESCE(fi.resuscitation_fluid_ml_24h, 0) END AS resuscitation_fluid_ml_24h,
    CASE WHEN fi.n_input_events > 0 THEN COALESCE(fi.crystalloid_fluid_ml_24h, 0) END AS crystalloid_fluid_ml_24h,
    CASE WHEN fi.n_input_events > 0 THEN COALESCE(fi.colloid_fluid_ml_24h, 0) END AS colloid_fluid_ml_24h,
    fi.fluid_input_ml_24h / NULLIF(w.weight_kg, 0) AS fluid_input_ml_kg_24h,
    CASE WHEN fi.fluid_input_ml_24h IS NOT NULL AND fo.fluid_output_ml_24h IS NOT NULL
         THEN (fi.fluid_input_ml_24h - fo.fluid_output_ml_24h) / NULLIF(w.weight_kg, 0) END AS fluid_balance_ml_kg_24h,
    COALESCE(ns.nee_any_24h, 0) AS nee_any_24h,
    COALESCE(ns.nee_max_24h, 0) AS nee_max_24h,
    COALESCE(ns.nee_twa_24h, 0) AS nee_twa_24h,
    COALESCE(ns.nee_auc_24h, 0) AS nee_auc_24h,
    COALESCE(ns.nee_duration_hours_24h, 0) AS nee_duration_hours_24h,
    vc.norepinephrine_max, vc.epinephrine_max, vc.phenylephrine_max, vc.dopamine_max,
    0::numeric AS metaraminol_max, vc.vasopressin_u_min_max, vc.angiotensin_ii_max,
    1::integer AS metaraminol_source_unavailable,
    0::integer AS angiotensin_ii_source_unavailable,
    COALESCE(ss.systemic_steroid_any_24h, 0) AS systemic_steroid_any_24h,
    COALESCE(ss.hydrocortisone_mg_24h, 0) AS hydrocortisone_mg_24h,
    COALESCE(ss.methylprednisolone_mg_24h, 0) AS methylprednisolone_mg_24h,
    COALESCE(ss.prednisolone_mg_24h, 0) AS prednisolone_mg_24h,
    COALESCE(ss.prednisone_mg_24h, 0) AS prednisone_mg_24h,
    COALESCE(ss.dexamethasone_mg_24h, 0) AS dexamethasone_mg_24h,
    COALESCE(ss.hydrocortisone_mg_24h, 0)
      + COALESCE(ss.methylprednisolone_mg_24h, 0) * 5.0
      + COALESCE(ss.prednisolone_mg_24h, 0) * 4.0
      + COALESCE(ss.dexamethasone_mg_24h, 0) * 26.67 AS hce_mg_24h,
    COALESCE(ss.hydrocortisone_mg_24h, 0)
      + COALESCE(ss.methylprednisolone_mg_24h, 0) * 5.0
      + (COALESCE(ss.prednisolone_mg_24h, 0) + COALESCE(ss.prednisone_mg_24h, 0)) * 4.0
      + COALESCE(ss.dexamethasone_mg_24h, 0) * 26.67 AS hce_with_prednisone_mg_24h,
    COALESCE(ss.steroid_unparsed_rows, 0) AS steroid_unparsed_rows,
    'administered_emar'::text AS steroid_source,
    'NE+E+phenylephrine/10+dopamine/100+metaraminol/8+vasopressin(U/min)*2.5+angiotensinII*10'::text AS nee_formula,
    'hydrocortisone+methylprednisolone*5+prednisolone*4+dexamethasone*26.67'::text AS hce_formula
FROM eligible AS e
LEFT JOIN weight AS w ON w.stay_id = e.stay_id
LEFT JOIN fluid_input AS fi ON fi.stay_id = e.stay_id
LEFT JOIN fluid_output AS fo ON fo.stay_id = e.stay_id
LEFT JOIN nee_summary AS ns ON ns.stay_id = e.stay_id
LEFT JOIN vaso_components AS vc ON vc.stay_id = e.stay_id
LEFT JOIN steroid_summary AS ss ON ss.stay_id = e.stay_id
ORDER BY e.stay_id;
