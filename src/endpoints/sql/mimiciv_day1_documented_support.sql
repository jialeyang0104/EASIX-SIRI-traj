WITH assigned(icu_stay_id,trajectory_class) AS (VALUES %s), cohort AS MATERIALIZED (SELECT i.*,a.trajectory_class FROM assigned a JOIN mimiciv_icu.icustays i
        ON i.stay_id=a.icu_stay_id), mv AS (
            SELECT e.stay_id,COUNT(v.stay_id) AS source_rows,
                COALESCE(MAX((v.ventilation_status='InvasiveVent')::integer) FILTER(WHERE v.stay_id IS NOT NULL),0) AS active
            FROM cohort e LEFT JOIN mimiciv_derived.ventilation v ON v.stay_id=e.stay_id
                AND v.starttime<e.intime+INTERVAL '24 hours' AND v.endtime>e.intime GROUP BY e.stay_id
        ), cr AS (
            SELECT e.stay_id,COUNT(c.stay_id) AS source_rows,
                COALESCE(MAX(((c.system_active=1 OR (c.system_active IS NULL AND c.blood_flow>0)))::integer) FILTER(WHERE c.stay_id IS NOT NULL),0) AS active
            FROM cohort e LEFT JOIN mimiciv_derived.crrt c ON c.stay_id=e.stay_id AND c.charttime>=e.intime AND c.charttime<e.intime+INTERVAL '24 hours' GROUP BY e.stay_id
        )
        SELECT e.stay_id AS icu_stay_id, o.mechvent_score AS oasis_mv_points,
            o.oasis-o.mechvent_score AS oasis_without_mv,
            mv.source_rows AS mv_source_rows,mv.active AS documented_mv_24h,
            cr.source_rows AS crrt_source_rows,cr.active AS documented_crrt_24h
        FROM cohort e LEFT JOIN mimiciv_derived.oasis o ON o.stay_id=e.stay_id
        LEFT JOIN mv ON mv.stay_id=e.stay_id LEFT JOIN cr ON cr.stay_id=e.stay_id ORDER BY e.stay_id
