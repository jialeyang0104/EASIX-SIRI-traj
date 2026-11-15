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

Replace `forest` with a selection in [the mapping](docs/figure_mapping.md). `--figure all` regenerates every group, including Figure 7, which requires locally regenerated `src/classifier/data_private/shap_internal.npz`. The patient-level SHAP arrays are deliberately not included. Output files are written to `results/figures/`; stable analysis-era filenames map to final submission labels through the mapping document. The local support module retains the article palette and PDF/SVG/300-dpi TIFF export settings without requiring private Codex skills.

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
environment/.venv/bin/python src/extract/extract_interval_profiles.py
environment/.venv/bin/python src/endpoints/scripts/prepare_management_inputs.py
environment/.venv/bin/python src/display/prepare_display_data.py
environment/.venv/bin/python src/display/prepare_flows.py
environment/.venv/bin/python src/display/prepare_outcome_displays.py
environment/.venv/bin/python src/display/refresh_model_aggregates.py
environment/.venv/bin/python src/display/render_figures.py --figure all
```

The source fit can take substantial time, particularly degree/pruning sensitivity and classifier selection/bootstrap. Existing generated models are reused by the original fitting scripts. `data/aggregate/` is the published reference set; the explicit refresh step replaces model aggregate sources with your new run. Review regenerated profiles and class mapping before comparing classes numerically.

## Verification and current status

```sh
environment/.venv/bin/python tests/check_public.py
Rscript tests/parse_r.R
```

This checkout has passed Python compilation, R parsing, controlled-output path scanning and synthetic ICU endpoint checks. All 11 aggregate-only renderer selections have exported successfully (PDF/SVG/TIFF/PNG); no image QA was performed. Licensed-data extraction, model refitting and the classifier/SHAP figure have not been run here. The source release/derived-concept version remains a data-custodian prerequisite for exact numerical replication.

## Data and license

[MIT](LICENSE) covers code and documentation. Controlled MIMIC source data remain subject to PhysioNet agreements. Patient data, fitted weights, individual probabilities, SHAP arrays, credentials, manuscript files and proprietary PDFs are not distributed. Generated private directories and `results/` are ignored by Git; keep those exclusions when extending the repository.
