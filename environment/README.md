# Runtime

Python 3.12 and R with PostgreSQL client `psql` are required. The original `uv.lock` and `renv.lock` preserve analysis environment provenance. The public Python manifest omits the unused PDF-reader dependency absent from the original lock; its direct requirements now match the lock. `DESCRIPTION` also declares the ICU-model dependencies. Install the additional database driver and explicit SciPy dependency after syncing:

```sh
uv sync --project environment --frozen
uv pip install --python environment/.venv/bin/python -r environment/requirements-extra.txt
Rscript -e 'if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv", repos="https://cloud.r-project.org"); renv::restore(lockfile="environment/renv.lock"); renv::install(c("data.table@1.18.6.1", "survival@3.8-12", "cmprsk@2.2-12", "jsonlite@2.0.0"))'
```

The R lock records the original trajectory-fitting package versions. The explicit installs apply the later endpoint/ICU runtime versions recorded in [analysis_runtime.json](analysis_runtime.json). A fresh machine also needs the platform system libraries required by these packages. Use `environment/.venv/bin/python` for the Python commands below. No patient-level computation has been rerun in this public checkout.

## Verification status

Python syntax, aggregate-source checks, synthetic ICU endpoint tests and the aggregate-only figure exports passed using the existing study Python environment recorded in `analysis_runtime.json`. R scripts passed parsing. A fresh `uv sync --frozen` installation was attempted but the dependency download did not complete; finish setup before using the commands above. Database extraction and model fitting were not rerun for this release.
