-- Read-only equality audit of the local MIMIC-IV NEE concept.
-- Formula matches the local mimic-code PostgreSQL concept; vasopressin is units/hour.
BEGIN READ ONLY;
COPY (
WITH reconstructed AS (
  SELECT stay_id, starttime, endtime,
    ROUND(CAST(
      COALESCE(norepinephrine, 0) + COALESCE(epinephrine, 0)
      + COALESCE(phenylephrine::double precision / 10, 0)
      + COALESCE(dopamine::double precision / 100, 0)
      + COALESCE((vasopressin * 2.5)::double precision / 60, 0)
      AS DECIMAL), 4) AS norepinephrine_equivalent_dose
  FROM mimiciv_derived.vasoactive_agent
  WHERE norepinephrine IS NOT NULL OR epinephrine IS NOT NULL
     OR phenylephrine IS NOT NULL OR dopamine IS NOT NULL OR vasopressin IS NOT NULL
), excess_actual AS (
  SELECT * FROM mimiciv_derived.norepinephrine_equivalent_dose
  EXCEPT ALL SELECT * FROM reconstructed
), excess_reconstructed AS (
  SELECT * FROM reconstructed
  EXCEPT ALL SELECT * FROM mimiciv_derived.norepinephrine_equivalent_dose
)
SELECT (SELECT count(*) FROM mimiciv_derived.norepinephrine_equivalent_dose) AS actual_rows,
       (SELECT count(*) FROM reconstructed) AS reconstructed_rows,
       (SELECT count(*) FROM excess_actual) AS actual_only_rows,
       (SELECT count(*) FROM excess_reconstructed) AS reconstructed_only_rows
) TO STDOUT WITH CSV HEADER;
COMMIT;
