"""Compare historical coefficients and full-junction block LASSO fits after Pfam rebuilding."""

from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import Lasso, LassoCV
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DATA = ROOT / "sinform_loco_full" / "data"
TARGET = "chr4:105898001:105927336"
COEF_PATH = {
    "expression": Path("D:/Phd_project/model_df/exp_lasso_feature_coefficients.csv"),
    "go": Path("D:/Phd_project/model_df/go_lasso_feature_coefficients.csv"),
    "pfam": Path("D:/Phd_project/pfam_go/pfam_coef.csv"),
}


def rank(scores: np.ndarray, ix: int) -> int:
    order = np.argsort(-scores, kind="stable")
    ranks = np.empty(len(scores), dtype=np.int32)
    ranks[order] = np.arange(1, len(scores) + 1)
    return int(ranks[ix])


def original_names(block: str) -> list[str]:
    if block == "pfam":
        header = pd.read_csv(HERE / "pfam_feature_map.csv")
        return header.legacy_raw_name.dropna().tolist()
    prefix = {"expression": "EXP__", "go": "GO__"}[block]
    return pd.read_csv(DATA / f"features_{block}.csv").feature.str.removeprefix(prefix).tolist()


def legacy_weights(block: str, names: list[str]) -> np.ndarray:
    coefficients = pd.read_csv(COEF_PATH[block]).set_index("feature").coefficient
    return coefficients.reindex(names).fillna(0).to_numpy(float)


def fit_block(label: str, X: sparse.csr_matrix, y: np.ndarray) -> tuple[np.ndarray, dict]:
    scaler = StandardScaler(with_mean=False)
    Xscaled = scaler.fit_transform(X.astype(np.float64))
    fixed_alphas = {
        "expression": 0.0018863663299020059,  # full-junction 5-fold CV
        "go": 0.000829523098237332,           # archived GO LASSO model
        "pfam_old": 0.0019986224674819718,     # archived Pfam LASSO model
        "pfam_new": 0.0019986224674819718,
    }
    use_fixed_alpha = "--historical-alpha" in sys.argv
    if use_fixed_alpha:
        model = Lasso(alpha=fixed_alphas[label], max_iter=10000,
                      tol=1e-4, selection="cyclic")
    else:
        cv = KFold(n_splits=5, shuffle=True, random_state=123)
        model = LassoCV(alphas=100, cv=cv, n_jobs=4, random_state=123,
                        max_iter=10000, tol=1e-4, selection="cyclic")
    print(f"Fitting {label}: {X.shape}, nnz={X.nnz}", flush=True)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always", ConvergenceWarning)
        model.fit(Xscaled, y)
    coef = np.asarray(model.coef_, dtype=np.float64)
    # Reproduce the archived SInFoRM score convention: standardized-fit
    # coefficients applied directly to raw input columns, without intercept.
    component = np.asarray(X @ coef).ravel()
    details = {
        "alpha": float(model.alpha if use_fixed_alpha else model.alpha_),
        "nonzero_coefficients": int(np.count_nonzero(coef)),
        "convergence_warnings": sum(issubclass(w.category, ConvergenceWarning) for w in captured),
        "cv_mse_min": float(np.min(model.mse_path_.mean(axis=1))) if not use_fixed_alpha else None,
    }
    np.save(HERE / f"coef_{label}.npy", coef)
    print(label, details, flush=True)
    return component, details


def main() -> None:
    rows = pd.read_csv(DATA / "rows.csv", usecols=["bed_id", "ZSTAT"])
    assert len(rows) == 95750 and rows.bed_id.is_unique
    ix = int(np.flatnonzero(rows.bed_id.to_numpy() == TARGET)[0])
    y = rows.ZSTAT.to_numpy(np.float64)
    matrices = {
        "expression": sparse.load_npz(DATA / "X_expression.npz").tocsr(),
        "go": sparse.load_npz(DATA / "X_go.npz").tocsr(),
        "pfam_old": sparse.load_npz(DATA / "X_pfam.npz").tocsr(),
        "pfam_new": sparse.csr_matrix(mmread(HERE / "X_pfam_rebuilt.mtx")),
    }
    assert all(X.shape[0] == len(rows) for X in matrices.values())
    assert matrices["pfam_new"].shape[1] == len(pd.read_csv(HERE / "pfam_feature_map.csv"))

    fixed = {}
    for block in ("expression", "go"):
        names = original_names(block)
        fixed[block] = np.asarray(matrices[block] @ legacy_weights(block, names)).ravel()
    old_names = original_names("pfam")
    old_weights = legacy_weights("pfam", old_names)
    fixed["pfam_old"] = np.asarray(matrices["pfam_old"] @ old_weights).ravel()
    fixed["pfam_new"] = np.asarray(matrices["pfam_new"][:, :len(old_names)] @ old_weights).ravel()
    old_fixed_score = np.abs(fixed["expression"] + fixed["go"] + fixed["pfam_old"])
    new_fixed_score = np.abs(fixed["expression"] + fixed["go"] + fixed["pfam_new"])
    historical = pd.read_excel("D:/Phd_project/model_df/common_df.xlsx").set_index("region").reindex(rows.bed_id)
    historical_diff = np.max(np.abs(old_fixed_score - historical["abs"].to_numpy(float)))
    assert historical_diff < 1e-7, historical_diff
    if "--fixed-only" in sys.argv:
        threshold = math.ceil(len(rows) * .1)
        for label, score in (("old", old_fixed_score), ("rebuilt", new_fixed_score)):
            print(label, "NPNT rank", rank(score, ix), "of", len(rows),
                  "top10", rank(score, ix) <= threshold,
                  "score", float(score[ix]), flush=True)
        return

    refit = {}
    fit_details = {}
    for label in ("expression", "go", "pfam_old", "pfam_new"):
        refit[label], fit_details[label] = fit_block(label, matrices[label], y)
    old_refit_score = np.abs(refit["expression"] + refit["go"] + refit["pfam_old"])
    new_refit_score = np.abs(refit["expression"] + refit["go"] + refit["pfam_new"])
    threshold = math.ceil(len(rows) * .1)
    score_sets = {
        "historical_coefficients_old_pfam": old_fixed_score,
        "historical_coefficients_new_pfam": new_fixed_score,
        "all_junction_refit_old_pfam": old_refit_score,
        "all_junction_refit_new_pfam": new_refit_score,
    }
    comparison = {}
    for label, score in score_sets.items():
        position = rank(score, ix)
        comparison[label] = {
            "npnt_rank": position,
            "npnt_percentile_top": float(100 * position / len(rows)),
            "npnt_score_absolute": float(score[ix]),
            "top_10_percent": bool(position <= threshold),
            "top_10_cutoff_score": float(np.partition(score, -threshold)[-threshold]),
        }
    frame = pd.DataFrame({"bed_id": rows.bed_id, **score_sets})
    frame.to_csv(HERE / "scores_comparison.csv", index=False)
    result = {
        "target": TARGET, "rows": len(rows), "top_10_rank_cutoff": threshold,
        "historical_score_reproduction_max_abs_difference": float(historical_diff),
        "training": ("Separate LASSO per block on all 95,750 junctions; StandardScaler variance scaling; raw-column dot products and absolute sum for ranking; no intercept in score. GO and Pfam penalties were selected in the archived original models; expression penalty was selected by five-fold CV on this full universe."
                     if "--historical-alpha" in sys.argv else
                     "Separate LassoCV per block on all 95,750 junctions; five shuffled folds, 100 alphas, StandardScaler variance scaling; raw-column dot products and absolute sum for ranking; no intercept in score."),
        "comparisons": comparison,
        "npnt_components": {
            "historical_new_pfam": {k: float(fixed[k][ix]) for k in ("expression", "go", "pfam_new")},
            "refit_new_pfam": {k: float(refit[k][ix]) for k in ("expression", "go", "pfam_new")},
        },
        "fit_details": fit_details,
        "interpretation_limit": "All coefficients were fit on the same full junction universe being ranked; these are descriptive in-sample ranks, not held-out performance.",
    }
    (HERE / "result.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
