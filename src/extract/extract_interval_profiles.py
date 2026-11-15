#!/usr/bin/env python3
"""Extract timestamp-correct MIMIC-IV clinical profiles for the frozen K=3 cohort.

The script uses a temporary database table containing only the frozen trajectory
assignment key and writes the row-level extract to the local ``private`` folder.
That folder is deliberately excluded from the public repository.  Public
outputs are produced by the aggregation script after this extraction succeeds.
"""

from __future__ import annotations

import io
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import psycopg2


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ASSIGNMENT_FILE = ROOT/'src/endpoints/data_private/all_icu_gbmtm_assignments_mimiciv_mimiciii.csv'
SQL_FILE = ROOT/'sql/same_interval_profile.sql'
PRIVATE_ROOT = ROOT/'data/private'
PRIVATE_ROOT.mkdir(parents=True,exist_ok=True)
OUTPUT_FILE = PRIVATE_ROOT/'mimiciv_same_interval_profiles_patient_level.csv'
META_FILE = PRIVATE_ROOT/'same_interval_extraction_metadata.json'

def fail(message: str) -> None:
    raise RuntimeError(message)


def load_assignments() -> pd.DataFrame:
    if not ASSIGNMENT_FILE.exists():
        fail(f"Frozen assignment file is missing: {ASSIGNMENT_FILE}")
    assignments = pd.read_csv(ASSIGNMENT_FILE, low_memory=False)
    required = {"dataset", "icu_stay_id", "gbmtm_class", "selected_k"}
    missing = required.difference(assignments.columns)
    if missing:
        fail(f"Assignment file missing columns: {sorted(missing)}")
    x = assignments.loc[
        (assignments["dataset"].astype(str) == "mimiciv")
        & (pd.to_numeric(assignments["selected_k"], errors="coerce") == 3),
        ["icu_stay_id", "gbmtm_class"],
    ].copy()
    x["icu_stay_id"] = pd.to_numeric(x["icu_stay_id"], errors="coerce")
    x["trajectory_class"] = pd.to_numeric(x["gbmtm_class"], errors="coerce")
    x = x.drop(columns="gbmtm_class")
    if x["icu_stay_id"].isna().any() or x["trajectory_class"].isna().any():
        fail("Frozen assignments contain non-numeric stay or class values")
    x["icu_stay_id"] = x["icu_stay_id"].astype("int64")
    x["trajectory_class"] = x["trajectory_class"].astype("int64")
    if x["icu_stay_id"].duplicated().any():
        fail("Frozen MIMIC-IV K=3 assignment has duplicate ICU stays")
    if not set(x["trajectory_class"]).issubset({1, 2, 3}):
        fail("Unexpected trajectory class outside 1, 2, 3")
    if len(x) != 2698:
        fail(f"Expected 2,698 frozen MIMIC-IV assignments, found {len(x)}")
    return x.sort_values("icu_stay_id").reset_index(drop=True)


def connect() -> psycopg2.extensions.connection:
    return psycopg2.connect(
        host=os.environ.get("PGHOST", "localhost"),
        port=int(os.environ.get("PGPORT", "5432")),
        dbname=os.environ.get("PGDATABASE", "mimiciv"),
        user=os.environ.get("PGUSER"),
        password=os.environ.get("PGPASSWORD") or None,
    )


def run_query(assignments: pd.DataFrame, n_intervals: int) -> pd.DataFrame:
    sql = SQL_FILE.read_text(encoding="utf-8")
    if "generate_series(1, 4)" not in sql:
        fail("Extraction SQL no longer contains its interval template")
    sql = sql.replace("generate_series(1, 4)", f"generate_series(1, {n_intervals})")
    connection = connect()
    try:
        with connection.cursor() as cursor:
            cursor.execute("SET statement_timeout = 0")
            cursor.execute("SET work_mem = '512MB'")
            cursor.execute(
                "CREATE TEMP TABLE tmp_trajectory_assignment "
                "(icu_stay_id bigint PRIMARY KEY, trajectory_class integer NOT NULL) "
                "ON COMMIT DROP"
            )
            buffer = io.StringIO()
            assignments[["icu_stay_id", "trajectory_class"]].to_csv(
                buffer, index=False, header=False
            )
            buffer.seek(0)
            cursor.copy_from(
                buffer,
                "tmp_trajectory_assignment",
                columns=("icu_stay_id", "trajectory_class"),
                sep=",",
            )
            cursor.execute("ANALYZE tmp_trajectory_assignment")
            cursor.execute(sql)
            rows = cursor.fetchall()
            columns = [item.name for item in cursor.description]
    finally:
        connection.close()
    result = pd.DataFrame(rows, columns=columns)
    expected_rows = len(assignments) * n_intervals
    if len(result) != expected_rows:
        fail(f"Expected {expected_rows} assigned interval rows, found {len(result)}")
    if result[["icu_stay_id", "interval_index"]].duplicated().any():
        fail("Extracted profile has duplicate stay/interval keys")
    class_counts = result.groupby("trajectory_class")["icu_stay_id"].nunique().to_dict()
    expected_counts = assignments.groupby("trajectory_class")["icu_stay_id"].nunique().to_dict()
    if class_counts != expected_counts:
        fail(f"Class counts changed during extraction: {class_counts} vs {expected_counts}")
    return result


def write_outputs(assignments: pd.DataFrame, result: pd.DataFrame, n_intervals: int) -> None:
    # PostgreSQL timestamps are kept as ISO strings for a lossless local audit.
    result.to_csv(OUTPUT_FILE, index=False, float_format="%.10g")
    metadata = {
        "analysis": "trajectory_clinical_management_profiles",
        "extraction_time_utc": datetime.now(timezone.utc).isoformat(),
        "database": os.environ.get("PGDATABASE", "mimiciv"),
        "database_host": os.environ.get("PGHOST", "localhost"),
        "database_port": int(os.environ.get("PGPORT", "5432")),
        "n_intervals": n_intervals,
        "interval_hours": 6,
        "interval_rule": "[intime + (k-1)*6 h, intime + k*6 h), k=1..28; boundary events enter the right interval",
        "trajectory_assignment_source": str(ASSIGNMENT_FILE),
        "trajectory_assignment_n": int(len(assignments)),
        "trajectory_assignment_class_counts": {
            str(k): int(v)
            for k, v in assignments["trajectory_class"].value_counts().sort_index().items()
        },
        "row_level_output": str(OUTPUT_FILE),
        "row_level_output_is_restricted": True,
        "row_count": int(len(result)),
        "source_tables": [
            "mimiciv_icu.icustays",
            "mimiciv_hosp.admissions",
            "mimiciv_hosp.patients",
            "mimiciv_derived.first_day_weight",
            "mimiciv_derived.norepinephrine_equivalent_dose",
            "mimiciv_icu.inputevents",
            "mimiciv_derived.vitalsign",
            "mimiciv_derived.urine_output",
            "mimiciv_hosp.labevents",
            "mimiciv_derived.ventilation",
            "mimiciv_derived.crrt",
            "mimiciv_derived.sofa",
        ],
    }
    META_FILE.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")


def main() -> None:
    assignments = load_assignments()
    n_intervals = int(os.environ.get("PROFILE_INTERVALS", "28"))
    if n_intervals not in {4, 28}:
        fail("PROFILE_INTERVALS must be 4 (early window) or 28 (D1-D7 extraction)")
    result = run_query(assignments, n_intervals)
    write_outputs(assignments, result, n_intervals)
    print(json.dumps({"rows": len(result), "output": str(OUTPUT_FILE)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
