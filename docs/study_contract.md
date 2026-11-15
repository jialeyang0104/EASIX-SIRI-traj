# Scientific contract

Adults, first ICU stay, ICU and hospital length of stay at least 24 hours; neither database requires admission-level `has_chartevents_data=1`. Daily index definitions and raw measurement units are encoded in `src/endpoints/sql/*daily.sql`.

Primary endpoint: all-cause death by ICU-admission day 28. Static models enter at 24 hours; trajectory models enter at day 7 among patients alive and still in ICU. Both use the same admission-based day-28 cutoff (27 and 21 remaining days). The retained scripts also generate historical 90-day sensitivity outputs.

Secondary endpoint: first-ICU death; live discharge is a competing event. Fine–Gray association models follow the entire first ICU stay, without 28/90-day administrative truncation. CIF display uses admission days 7–60. Unresolved death-date ordering remains unknown and is excluded from time-to-event fitting, not recoded as live discharge.

Trajectory training uses MIMIC-IV; two or more observed positive days per index, log2 transform, IV 1st/99th percentile winsorization, linear interpolation with boundary carry, fixed IV scaling. Candidate K=1–6 models use quadratic degree and pruning; the selected model is projected into III without refitting. Degree 1–4 with/without pruning and five-seed stability are separate diagnostics. No fitted weights or individual labels are distributed. Fresh fitting may permute class labels or change convergence; compare fitted profiles and aggregate candidate metrics before claiming exact replication.

The final Day-1 classifier uses patient/class-stratified IV 70/30 training/internal validation and fixed III external evaluation, not the historical crossfit bridge. It considers 55 candidates; training-only consensus selection originally retained 22. Preprocessing, selection and tuning stay inside the training data. The clinical profile display contains 42 features.

Cox adjusts age, sex, race/ethnicity and Day-1 OASIS. ICU models combine Other/Unknown race after support checks. Seven subgroup interactions are BH-adjusted within each database and endpoint/horizon. Treatment analyses retain retrospective future-trajectory stratification; the original fixed cutoffs and missing-source semantics remain in the code.

## Source provenance

The original extraction record dated 2026-07-30 used PostgreSQL 16.13. The local databases did not preserve a machine-readable PhysioNet release identifier; exact MIMIC release and derived-concept commit must be recovered from the data custodian's download/build records. The finalized manuscript reports MIMIC-IV 3.1; the local download metadata and derived-concept commit have not independently confirmed that release attribution. The required derived schema is discoverable from SQL references; build it before extraction. This missing release provenance limits exact numerical replication on a new database installation.

`source_provenance.tsv` maps copied study code to its original project-relative source. Public edits change paths, remove archive restoration and private skill dependencies, and add explicit extraction/checking entrypoints; scientific formulas are retained.
