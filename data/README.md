# Controlled-access data

Obtain MIMIC-IV and MIMIC-III through PhysioNet credentialing and their data-use agreements. Load PostgreSQL databases `mimiciv` (port 5432) and `mimiciii` (5433); build the derived concepts referenced by the extraction SQL. Configure authentication using normal PostgreSQL mechanisms, outside this repository. Queries run read-only.

Source schemas: `mimiciv_hosp`, `mimiciv_icu`, `mimiciv_derived` and `mimiciii`. The SQL files are the field/item/unit contracts; [required relation inventory](../docs/database_tables.md). Database release and concept-build versions affect exact reproduction; see the source-version notes in `docs/study_contract.md`.

Extraction creates `data/private/`, source caches in `src/endpoints/data_private/`, and private model/intermediate outputs under `results/`. Classifier probabilities, ICU timestamps, identifiers, SHAP arrays and fitted objects stay local and are ignored by Git.
