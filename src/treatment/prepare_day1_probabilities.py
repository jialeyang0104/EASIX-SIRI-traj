"""Generate local, out-of-sample probabilities for the treatment analysis."""
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.model_selection import StratifiedKFold

SRC = Path(__file__).resolve().parents[1]
PRIVATE = SRC / "treatment/data_private"
CLASSIFIER = SRC / "classifier"
ENDPOINTS = SRC / "endpoints/data_private"
KEY = ["dataset", "icu_stay_id", "subject_id"]
LABEL = "gbmtm_class"
PROBS = [f"prob_class_{k}" for k in range(1, 4)]

def load_day1(database, labels):
    d = pd.read_csv(ENDPOINTS / f"all_icu_d1d7_daily_{database}.csv", low_memory=False).query("day == 1").copy()
    cells = ["neutrophils", "monocytes", "lymphocytes"]
    absolute = d[[c + "_abs" for c in cells]].notna().all(axis=1) & (database == "mimiciv")
    for cell in cells:
        d[cell + "_harmonised"] = d[cell + "_abs"].where(absolute, d.wbc_109_l * d[cell + "_pct"] / 100)
    d["male"] = d.sex.eq("M").astype(int)
    d = d.merge(labels[labels.dataset.eq(database)], on=KEY[:2], validate="one_to_one")
    systems = pd.read_csv(ENDPOINTS / f"{database}_day1_system_means.csv")
    d = d.merge(systems, on="icu_stay_id", how="left", validate="one_to_one", indicator=True)
    assert d._merge.eq("both").all() and d.subject_id.is_unique
    return d.drop(columns="_merge")

def main():
    spec = json.loads((CLASSIFIER / "aggregate/classifier_manifest.json").read_text())
    saved = joblib.load(CLASSIFIER / "data_private/day1_classifier.joblib")
    features, template = saved["features"], saved["pipeline"]
    assert features == spec["selected_features"]
    split = pd.read_csv(CLASSIFIER / "data_private/iv_patient_split.csv").assign(dataset="mimiciv")
    labels = pd.read_csv(ENDPOINTS / "all_icu_gbmtm_assignments_mimiciv_mimiciii.csv")
    labels = labels.loc[labels.selected_k.eq(3), KEY[:2] + [LABEL]]
    stored = pd.read_csv(CLASSIFIER / "data_private/day1_fixed_predictions.csv")
    validation = {"mimiciv": "MIMIC-IV internal validation", "mimiciii": "MIMIC-III external validation"}
    results = []
    for db, cohort in validation.items():
        d = load_day1(db, labels)
        for c in PROBS:
            d[c] = np.nan
        d["fold"] = pd.Series(pd.NA, index=d.index, dtype="Int64")
        d["prediction_source"] = "frozen_model_external"
        if db == "mimiciv":
            d = d.merge(split[KEY + [LABEL, "partition"]], on=KEY + [LABEL], validate="one_to_one")
            train = d.partition.eq("train")
            hold = d.partition.eq("internal validation")
            assert (train | hold).all() and len(d) == len(split)
            assert set(d.loc[train, "subject_id"]).isdisjoint(d.loc[hold, "subject_id"])
            training = d.loc[train].sort_values("subject_id")
            cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=spec["seed"])
            for fold, (i, j) in enumerate(cv.split(training[features], training[LABEL]), 1):
                fit, held = training.iloc[i], training.iloc[j]
                assert set(fit.subject_id).isdisjoint(held.subject_id)
                estimator = clone(template).fit(fit[features], fit[LABEL] - 1)
                assert np.array_equal(estimator.classes_, np.arange(3))
                d.loc[held.index, PROBS] = estimator.predict_proba(held[features])
                d.loc[held.index, "fold"] = fold
            assert d.loc[train, "fold"].notna().all()
            d.loc[train, "prediction_source"] = "fixed_recipe_patient_OOF"
            d.loc[hold, "prediction_source"] = "frozen_model_internal"
        else:
            hold = pd.Series(True, index=d.index)
        held = d.loc[hold]
        reference = held[KEY + [LABEL]].merge(
            stored.loc[stored.cohort.eq(cohort), KEY[1:] + [LABEL] + PROBS],
            on=KEY[1:] + [LABEL], validate="one_to_one", sort=False)
        assert np.array_equal(reference.icu_stay_id, held.icu_stay_id)
        pp = template.predict_proba(held[features])
        assert np.allclose(pp, reference[PROBS], rtol=0, atol=2e-7)
        d.loc[hold, PROBS] = pp
        p = d[PROBS].to_numpy()
        assert np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all()
        assert np.allclose(p.sum(axis=1), 1, rtol=0, atol=1e-6)
        results.append(d[KEY + ["prediction_source", "fold"] + PROBS])
    PRIVATE.mkdir(exist_ok=True)
    probabilities = pd.concat(results, ignore_index=True)
    probabilities.to_csv(PRIVATE / "day1_probabilities.csv", index=False)
    print(probabilities.groupby(["dataset", "prediction_source"]).size().to_string())

if __name__ == "__main__":
    main()
