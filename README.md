# EASIX–SIRI trajectories in critically ill adults

Reproducibility code for joint seven-day EASIX/SIRI trajectories, 28-day all-cause mortality, first-ICU mortality, clinical profiles, and a Day-1 trajectory classifier in MIMIC-IV and MIMIC-III.

## Start here

- [Scientific definitions and source provenance](docs/study_contract.md)
- [Ordered producer/consumer audit](docs/dataflow.md)
- [Final article figure/table mapping](docs/figure_mapping.md)
- [Aggregate display sources](data/aggregate/) — no patient data required to inspect these tables
- [Runtime setup](environment/README.md) and [controlled-data prerequisites](data/README.md)

| Path | Purpose | Inputs and outputs | Entry point |
|---|---|---|---|
| `src/extract/`, `sql/` | Licensed-data extraction | PostgreSQL → ignored private caches | `extract_sources.py` |
| `src/trajectory/` | K selection, GBMTM and diagnostics | Daily inputs → local models and aggregate metrics | `fit_gbmtm.R` |
| `src/endpoints/` | Current daily SQL, Cox, subgroup/treatment models | Private inputs → local aggregates | commands below |
| `src/icu/` | Competing-risk endpoint and models | Private caches → Fine–Gray/CIF aggregates | `prepare_analysis.py` |
| `src/treatment/` | Exploratory first-day treatment models for Figures S6 and S7 | Local ICU and classifier caches → local estimates and probability support | `prepare_day1_probabilities.py`, `fit_associations.R` |
| `src/classifier/` | Final 70/30 classifier | Day-1 inputs and trajectory labels → metrics, private model/SHAP arrays | `train_day1_classifier.py` |
| `src/display/` | Original statistical figure composition, local export helpers | Released aggregates and local SHAP → `results/figures/` | `render_figures.py` |
| `data/aggregate/` | Released reference aggregate sources | Frozen study estimates → figure display | inspect CSVs |
| `config/` | ICU endpoint/display contract | JSON → ICU analysis/display | `icu_analysis_contract.json` |
| `environment/`, `tests/` | Environment provenance and public checks | Code/aggregates → checks | commands below |

## Reproduce displays from released aggregates

After environment setup, run from the repository root:

```sh
export PROJECT_ROOT="$PWD"
environment/.venv/bin/python src/display/render_figures.py --figure forest
```

Replace `forest` with a selection in [the mapping](docs/figure_mapping.md). Use `--figure treatment` for the current fixed-class Figure S6 or `--figure treatment-probability` for the four-panel Figure S7. Both render from the released aggregates without patient data. See the [figure legends](docs/figure_mapping.md#exploratory-treatment-figure-legends) for the contrasts and tests. `--figure all` regenerates every group, including Figure 7, which requires locally regenerated `src/classifier/data_private/shap_internal.npz`. The patient-level SHAP arrays are deliberately not included. Output files are written to `results/figures/`; stable analysis-era filenames map to final submission labels through the mapping document. The local support module retains the article palette and PDF/SVG/300-dpi TIFF export settings without requiring private Codex skills.

## Recompute from licensed databases

Install the runtime and configure PostgreSQL authentication first. These commands create patient-level local files; run them only in an approved MIMIC workspace. They do not run manuscript/Word production.

```sh
export PROJECT_ROOT="$PWD"
environment/.venv/bin/python src/extract/extract_sources.py
environment/.venv/bin/python src/extract/clean_treatments.py
Rscript src/trajectory/fit_gbmtm.R
Rscript src/trajectory/audit_multistart.R
Rscript src/trajectory/fit_sensitivity.R
Rscript src/endpoints/scripts/reproject_gbmtm.R
environment/.venv/bin/python src/endpoints/scripts/extract_day1_systems.py
Rscript src/endpoints/scripts/05_58_v5_landmark.R
environment/.venv/bin/python src/endpoints/scripts/check_results.py
environment/.venv/bin/python src/endpoints/scripts/extract_supplementary_treatments.py
Rscript src/endpoints/scripts/fit_supplementary_treatments.R
environment/.venv/bin/python src/icu/prepare_analysis.py
Rscript src/icu/fit_competing_risks.R
Rscript src/icu/prepare_cif_display.R
environment/.venv/bin/python src/classifier/train_day1_classifier.py
environment/.venv/bin/python src/treatment/prepare_day1_probabilities.py
Rscript src/treatment/fit_associations.R
environment/.venv/bin/python src/extract/extract_interval_profiles.py
environment/.venv/bin/python src/endpoints/scripts/prepare_management_inputs.py
environment/.venv/bin/python src/display/prepare_display_data.py
environment/.venv/bin/python src/display/prepare_flows.py
environment/.venv/bin/python src/display/prepare_outcome_displays.py
environment/.venv/bin/python src/display/refresh_model_aggregates.py
environment/.venv/bin/python src/display/render_figures.py --figure all --treatment-source refit
```

The source fit can take substantial time, particularly degree/pruning sensitivity and classifier selection/bootstrap. Existing generated models are reused by the original fitting scripts. `data/aggregate/` is the published reference set; the explicit refresh step replaces the earlier endpoint/classifier/trajectory reference tables with your new run. S6 and S7 refit outputs stay under `src/treatment/aggregate/`; the renderer uses them only with `--treatment-source refit`. Review regenerated profiles and class mapping before comparing classes numerically.

## Verification and current status

```sh
environment/.venv/bin/python tests/check_public.py
Rscript tests/parse_r.R
```

The current S6 and S7 figures were rendered from the released aggregates. Their treatment models were also reproduced against the original controlled caches: 48 fixed-class rows (46 estimable), 16 probability-model tests, 3,232 curve points and 1,280 support-bin counts matched the released sources. The release checks cover aggregate schemas, synthetic ICU endpoints and source parsing. The source release and derived-concept build remain prerequisites for exact replication from a fresh database.

## Data and license

[MIT](LICENSE) covers code and documentation. Controlled MIMIC source data remain subject to PhysioNet agreements. Patient data, fitted weights, individual probabilities, SHAP arrays, credentials, manuscript files and proprietary PDFs are not distributed. Generated private directories and `results/` are ignored by Git; keep those exclusions when extending the repository.
