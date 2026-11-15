#!/usr/bin/env python3
"""Apply post-extraction plausibility rules to Day 1 treatment tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "data/private"
REPORT = ROOT / "results/treatment_cleaning.json"


def clean(path: Path) -> dict[str, int | str]:
    data = pd.read_csv(path, low_memory=False)
    if data["icu_stay_id"].duplicated().any():
        raise RuntimeError(f"Duplicate ICU stays in {path.name}")

    weight = pd.to_numeric(data["weight_kg"], errors="coerce")
    if "weight_raw_kg" not in data:
        data.insert(data.columns.get_loc("weight_kg"), "weight_raw_kg", weight)
    valid_weight = weight.between(20, 300, inclusive="both")
    data["weight_kg"] = weight.where(valid_weight)

    fluid_input = pd.to_numeric(data["fluid_input_ml_24h"], errors="coerce")
    fluid_balance = pd.to_numeric(data["fluid_balance_ml_24h"], errors="coerce")
    data["fluid_input_ml_kg_24h"] = fluid_input / data["weight_kg"]
    data["fluid_balance_ml_kg_24h"] = fluid_balance / data["weight_kg"]
    data.replace([np.inf, -np.inf], np.nan, inplace=True)

    temp = path.with_suffix(".cleaning.csv")
    data.to_csv(temp, index=False)
    temp.replace(path)
    return {
        "file": path.name,
        "rows": len(data),
        "unique_icu_stays": int(data["icu_stay_id"].nunique()),
        "invalid_weight_rows_set_missing": int((weight.notna() & ~valid_weight).sum()),
        "valid_weight_rows": int(data["weight_kg"].notna().sum()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", type=Path)
    args = parser.parse_args()
    paths = args.paths or [
        TABLES / "all_icu_day1_treatments_mimiciv.csv",
    ]
    missing = [str(path) for path in paths if not path.exists() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing or empty treatment extracts: " + ", ".join(missing))
    results = [clean(path.resolve()) for path in paths]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(results, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
