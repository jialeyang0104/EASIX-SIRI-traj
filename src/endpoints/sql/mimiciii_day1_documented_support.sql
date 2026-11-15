WITH assigned(icu_stay_id,trajectory_class) AS (VALUES %s), cohort AS MATERIALIZED (SELECT i.*,a.trajectory_class FROM assigned a JOIN mimiciii.icustays i
        ON i.icustay_id=a.icu_stay_id), mv AS (
            SELECT e.icustay_id,COUNT(v.icustay_id) AS source_rows,
                COALESCE(MAX((TRUE)::integer) FILTER(WHERE v.icustay_id IS NOT NULL),0) AS active
            FROM cohort e LEFT JOIN mimiciii.ventilation_durations v ON v.icustay_id=e.icustay_id
                AND v.starttime<e.intime+INTERVAL '24 hours' AND v.endtime>e.intime GROUP BY e.icustay_id
        ), cr AS (
            SELECT e.icustay_id,COUNT(c.icustay_id) AS source_rows,
                COALESCE(MAX((TRUE)::integer) FILTER(WHERE c.icustay_id IS NOT NULL),0) AS active
            FROM cohort e LEFT JOIN mimiciii.crrt_durations c ON c.icustay_id=e.icustay_id AND c.starttime<e.intime+INTERVAL '24 hours' AND c.endtime>e.intime GROUP BY e.icustay_id
        )
        SELECT e.icustay_id AS icu_stay_id, o.mechvent_score AS oasis_mv_points,
            o.oasis-o.mechvent_score AS oasis_without_mv,
            mv.source_rows AS mv_source_rows,mv.active AS documented_mv_24h,
            cr.source_rows AS crrt_source_rows,cr.active AS documented_crrt_24h
        FROM cohort e LEFT JOIN mimiciii.oasis o ON o.icustay_id=e.icustay_id
        LEFT JOIN mv ON mv.icustay_id=e.icustay_id LEFT JOIN cr ON cr.icustay_id=e.icustay_id ORDER BY e.icustay_id
