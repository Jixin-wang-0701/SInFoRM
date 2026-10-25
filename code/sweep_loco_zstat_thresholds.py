"""Evaluate fixed LOCO scores at stronger ZSTAT-positive definitions."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
THRESHOLDS = (1.96, 2.5, 3.0, 4.0, 5.0)
SOURCES = {
    "full_feature_fixed_alpha": (
        HERE / "loco_full_feature_fixed_alpha" / "oof_predictions.csv",
        "pfam_new_component", "rebuilt_score",
    ),
    "fold_selected_min30": (
        HERE / "loco_rebuilt_pfam" / "oof_predictions.csv",
        "rebuilt_pfam_component", "rebuilt_feature_only_score",
    ),
    "half_training_chromosomes": (
        HERE / "loco_half_chromosome_screen" / "oof_predictions.csv",
        "rebuilt_pfam_component", "rebuilt_feature_only_score",
    ),
    "threequarter_training_chromosomes": (
        HERE / "loco_threequarter_chromosome_screen" / "oof_predictions.csv",
        "rebuilt_pfam_component", "rebuilt_feature_only_score",
    ),
}


def metrics(positive: np.ndarray, score: np.ndarray) -> dict:
    n = len(score)
    k = math.ceil(n * 0.1)
    order = np.argsort(-score, kind="stable")[:k]
    boundary = np.partition(score, -k)[-k]
    above = score > boundary
    tied = score == boundary
    tie_slots = k - int(above.sum())
    expected_hits = float(positive[above].sum()) + tie_slots * float(positive[tied].mean())
    prevalence = float(positive.mean())
    ap = float(average_precision_score(positive, score))
    return {
        "n": n,
        "positive_n": int(positive.sum()),
        "positive_prevalence": prevalence,
        "auprc": ap,
        "auprc_over_prevalence": ap / prevalence,
        "auroc": float(roc_auc_score(positive, score)),
        "top10_n": k,
        "top10_positive_hits_stable_ties": int(positive[order].sum()),
        "top10_positive_hits_expected_random_ties": expected_hits,
        "sensitivity_at_10_expected_random_ties": expected_hits / int(positive.sum()),
        "top10_precision_expected_random_ties": expected_hits / k,
        "top10_boundary_tied_n": int(tied.sum()),
        "zero_score_fraction": float(np.mean(score == 0)),
    }


def main() -> None:
    flags = pd.read_csv(ROOT / "sinform_loco_full" / "data" / "rows.csv",
                        usecols=["bed_id", "ZSTAT", "has_expression", "has_go"])
    common_ids = set(flags.loc[flags.has_expression & flags.has_go, "bed_id"])
    output = []
    for model, (path, pfam_column, combined_column) in SOURCES.items():
        frame = pd.read_csv(path)
        assert len(frame) == len(flags) == 95750 and frame.bed_id.is_unique
        frame = frame.merge(flags, on="bed_id", validate="one_to_one",
                            suffixes=("", "_source"))
        assert np.allclose(frame.ZSTAT, frame.ZSTAT_source)
        frame = frame.drop(columns="ZSTAT_source")
        for subset, part in (("all_95750", frame),
                             ("common_expression_go_pfam", frame[frame.bed_id.isin(common_ids)])):
            score_columns = {
                "expression": np.abs(part.expression_component.to_numpy(float)),
                "go": np.abs(part.go_component.to_numpy(float)),
                "pfam_rebuilt": np.abs(part[pfam_column].to_numpy(float)),
                "combined": part[combined_column].to_numpy(float),
            }
            for threshold in THRESHOLDS:
                positive = (part.ZSTAT.to_numpy(float) > threshold)
                for block, score in score_columns.items():
                    output.append({"model": model, "subset": subset,
                                   "endpoint": "signed_ZSTAT", "threshold": threshold,
                                   "block": block, **metrics(positive, score)})
    pd.DataFrame(output).to_csv(HERE / "loco_threshold_sweep.csv", index=False)

    egf = pd.read_csv(HERE / "egf_ca_single_feature_loco_predictions.csv")
    abs_egf = pd.read_csv(HERE / "egf_ca_absolute_zstat_loco_predictions.csv",
                          usecols=["bed_id", "oof_feature_only_abs_target_score"])
    egf = egf.merge(abs_egf, on="bed_id", validate="one_to_one")
    assert len(egf) == len(flags) and egf.bed_id.is_unique
    output = []
    for subset, part in (("all_95750", egf),
                         ("common_expression_go_pfam", egf[egf.bed_id.isin(common_ids)])):
        for threshold in THRESHOLDS:
            z = part.ZSTAT.to_numpy(float)
            for endpoint, positive, score_column in (
                ("signed_ZSTAT", z > threshold, "oof_feature_only_abs_score"),
                ("absolute_ZSTAT", np.abs(z) > threshold,
                 "oof_feature_only_abs_target_score"),
            ):
                score = part[score_column].to_numpy(float)
                output.append({"model": "EGF_CA_single_feature_LOCO",
                               "subset": subset, "endpoint": endpoint,
                               "threshold": threshold, "block": "EGF_CA",
                               **metrics(positive, score)})
    pd.DataFrame(output).to_csv(HERE / "egf_ca_threshold_sweep.csv", index=False)
    print("Wrote", len(pd.read_csv(HERE / "loco_threshold_sweep.csv")),
          "model rows and", len(output), "EGF_CA rows")


if __name__ == "__main__":
    main()
