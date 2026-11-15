# Ordered dataflow audit

The README sequence supplies every consumed private cache:

| Producer | Private products | Consumers |
|---|---|---|
| `extract_sources.py` | Both daily caches, raw endpoints, IV treatment extract | trajectory fit, systems/Cox/ICU/classifier, treatment cleaning |
| `clean_treatments.py` | Plausibility-clean IV treatment extract | supplementary treatments, parent management |
| `fit_gbmtm.R` | Local K1–6/selected models, preprocessing, assignments, prepared training | sensitivity, stability, fixed projection |
| `reproject_gbmtm.R` | Current assignments and dual-index eligibility | Cox, treatment, ICU, classifier, profiles |
| `extract_day1_systems.py` | Both first-day system means | classifier, clinical displays |
| `05_58_v5_landmark.R` and `check_results.py` | Scale and landmark aggregate contracts | ICU preparation, cohort displays |
| `extract_supplementary_treatments.py` | Current assigned Day-1 treatments and III six-hour source | treatment models, ICU preparation |
| `icu/prepare_analysis.py` | ICU endpoint analysis table | Fine–Gray, CIF, clinical flow preparation |
| `train_day1_classifier.py` | Trained classifier, individual probabilities, internal/external SHAP arrays | Figure 7 (internal SHAP only) |
| `extract_interval_profiles.py` | IV same-interval patient profiles | parent/six-hour management preparation |
| `prepare_management_inputs.py` | Current parent treatments and both current six-hour tables | clinical display preparation |

Display metadata (`display_variables.csv` and related label/system definitions) is distributed in `data/aggregate/display/`. Preparation overwrites parent/class baseline rows and availability directly from the current caches; it no longer appends clinical rows or availability entries from an older result. `phenotypes()` creates `subtype_radar.csv` before `system_radar()` reads it. `prepare_display_data.py` produces display source tables before flow/outcome binding. The final refresh copies model CSVs and the classifier manifest to released display paths; rendering then uses the refreshed aggregates plus locally generated internal SHAP.

The released aggregates support display-only use without executing the private-data sequence. No licensed-data model run was performed during this dataflow audit.
