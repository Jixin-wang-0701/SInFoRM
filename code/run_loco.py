"""Chromosome-held-out SInFoRM analysis on the original 95,750-junction universe.

All supervised screening, scaling, baseline fitting and LASSO tuning are fitted
within each outer training set. Inner tuning holds out entire chromosomes.
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.metrics import average_precision_score, r2_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler


BLOCKS = ("expression", "go", "pfam")
ALPHA_RATIOS = (0.05, 0.2, 0.8)
MIN_NONZERO = {"expression": 5, "go": 30, "pfam": 30}
TOP_K = {"expression": 17, "go": 100, "pfam": 100}
BASELINE_RIDGE_ALPHA = 100.0
SEED = 20260615


def baseline_matrix(rows: pd.DataFrame) -> np.ndarray:
    length = np.maximum(1.0, rows.j_end.to_numpy() - rows.j_start.to_numpy())
    midpoint = (rows.j_start.to_numpy() + rows.j_end.to_numpy()) / 2
    columns = [np.log1p(length), np.log1p(midpoint)]
    names = ["log_junction_length", "log_chromosome_midpoint"]
    for name in (
        "magma_nsnps", "magma_nparam", "overlap_gene_count",
        "overlap_protein_coding_count", "overlap_gene_size_mean",
        "overlap_gene_size_max", "local_gene_density_1mb",
        "local_protein_coding_density_1mb", "nearest_gene_distance",
    ):
        columns.append(np.log1p(np.maximum(0, rows[name].to_numpy(dtype=float))))
        names.append("log1p_" + name)
    raw = np.column_stack(columns)
    assert np.isfinite(raw).all()
    return np.column_stack((raw, raw**2))


def fit_baseline(B: np.ndarray, y: np.ndarray, train: np.ndarray, test: np.ndarray):
    scaler = StandardScaler().fit(B[train])
    model = Ridge(alpha=BASELINE_RIDGE_ALPHA).fit(scaler.transform(B[train]), y[train])
    return model.predict(scaler.transform(B[train])), model.predict(scaler.transform(B[test]))


def screen_block(X: sparse.csr_matrix, train: np.ndarray, target: np.ndarray,
                 block: str) -> tuple[np.ndarray, np.ndarray]:
    Xtr = X[train]
    n = len(train)
    sx = np.asarray(Xtr.sum(axis=0)).ravel().astype(np.float64)
    sx2 = np.asarray(Xtr.power(2).sum(axis=0)).ravel().astype(np.float64)
    sy = float(target.sum())
    vary = float(target @ target - sy * sy / n)
    sxy = np.asarray(Xtr.T @ target).ravel().astype(np.float64)
    varx = sx2 - sx * sx / n
    correlation = (sxy - sx * sy / n) / np.sqrt(np.maximum(varx * vary, 1e-20))
    nonzero = np.asarray(Xtr.getnnz(axis=0)).ravel()
    eligible = np.isfinite(correlation) & (varx > 1e-12) & (
        nonzero >= MIN_NONZERO[block]) & (nonzero <= n - MIN_NONZERO[block])
    indices = np.flatnonzero(eligible)
    indices = indices[np.argsort(-np.abs(correlation[indices]), kind="stable")]
    indices = indices[: TOP_K[block]]
    return indices, correlation[indices]


def screened_matrix(matrices: dict[str, sparse.csr_matrix], blocks: tuple[str, ...],
                    train: np.ndarray, target: np.ndarray):
    parts = []
    detail = []
    for block in blocks:
        selected, correlations = screen_block(matrices[block], train, target, block)
        if len(selected):
            parts.append(matrices[block][:, selected])
            detail.extend((block, int(j), float(corr))
                          for j, corr in zip(selected, correlations))
    if not parts:
        return sparse.csr_matrix((next(iter(matrices.values())).shape[0], 0)), detail
    return sparse.hstack(parts, format="csr"), detail


def relative_alpha(X: sparse.csr_matrix, target: np.ndarray, l1_ratio: float) -> float:
    centered = target - target.mean()
    cross = np.asarray(X.T @ centered).ravel()
    return max(float(np.max(np.abs(cross))) / (len(target) * l1_ratio), 1e-8)


def fit_regularized(X: sparse.csr_matrix, target: np.ndarray,
                    alpha_ratio: float, l1_ratio: float):
    scaler = StandardScaler(with_mean=False).fit(X)
    Xs = scaler.transform(X)
    alpha = relative_alpha(Xs, target, l1_ratio) * alpha_ratio
    model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, max_iter=3000,
                       tol=1e-3, random_state=SEED, selection="cyclic")
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always", ConvergenceWarning)
        model.fit(Xs, target)
    nonconverged = any(issubclass(w.category, ConvergenceWarning) for w in captured)
    return scaler, model, alpha, nonconverged


def tune_model(matrices: dict[str, sparse.csr_matrix], B: np.ndarray,
               y: np.ndarray, chromosome: np.ndarray, outer_train: np.ndarray,
               blocks: tuple[str, ...], target_kind: str, l1_ratio: float):
    splitter = GroupKFold(n_splits=3)
    scores = np.zeros((len(ALPHA_RATIOS), 3), dtype=float)
    for inner_fold, (train_rel, valid_rel) in enumerate(
        splitter.split(outer_train, groups=chromosome[outer_train])
    ):
        inner_train = outer_train[train_rel]
        inner_valid = outer_train[valid_rel]
        if target_kind == "residual":
            base_train, base_valid = fit_baseline(B, y, inner_train, inner_valid)
            train_target = y[inner_train] - base_train
            valid_target = y[inner_valid] - base_valid
            positive = valid_target >= np.quantile(train_target, 0.9)
        else:
            train_target = y[inner_train]
            positive = y[inner_valid] > 1.96
        Xall, _ = screened_matrix(matrices, blocks, inner_train, train_target)
        if Xall.shape[1] == 0:
            scores[:, inner_fold] = positive.mean()
            continue
        Xtrain = Xall[inner_train]
        Xvalid = Xall[inner_valid]
        for j, ratio in enumerate(ALPHA_RATIOS):
            scaler, model, _, _ = fit_regularized(Xtrain, train_target, ratio, l1_ratio)
            prediction = model.predict(scaler.transform(Xvalid))
            scores[j, inner_fold] = average_precision_score(positive, prediction)
    mean_scores = scores.mean(axis=1)
    chosen = int(np.argmax(mean_scores))
    return ALPHA_RATIOS[chosen], mean_scores, scores


def fit_candidate(matrices: dict[str, sparse.csr_matrix], B: np.ndarray,
                  y: np.ndarray, chromosome: np.ndarray,
                  train: np.ndarray, test: np.ndarray,
                  blocks: tuple[str, ...], target_kind: str, l1_ratio: float,
                  feature_names: dict[str, np.ndarray]):
    ratio, mean_cv, fold_cv = tune_model(matrices, B, y, chromosome, train,
                                         blocks, target_kind, l1_ratio)
    if target_kind == "residual":
        baseline_train, baseline_test = fit_baseline(B, y, train, test)
        target = y[train] - baseline_train
    else:
        baseline_test = None
        target = y[train]
    Xall, detail = screened_matrix(matrices, blocks, train, target)
    if Xall.shape[1] == 0:
        prediction = np.full(len(test), target.mean())
        component = np.zeros(len(test))
        coefficients = np.array([])
        alpha = np.nan
        nonconverged = False
    else:
        scaler, model, alpha, nonconverged = fit_regularized(
            Xall[train], target, ratio, l1_ratio)
        prediction = model.predict(scaler.transform(Xall[test]))
        # This reproduces the legacy score convention: standardized-model
        # coefficients applied to raw annotation columns, with no intercept.
        component = np.asarray(Xall[test] @ model.coef_).ravel()
        coefficients = model.coef_
    feature_rows = [
        {"block": block, "feature": feature_names[block][index],
         "screen_correlation": corr, "coefficient_scaled": float(coef),
         "nonzero": bool(abs(coef) > 1e-10)}
        for (block, index, corr), coef in zip(detail, coefficients)
    ]
    tuning = {
        "blocks": "+".join(blocks), "target_kind": target_kind,
        "l1_ratio": l1_ratio, "alpha_ratio": ratio, "alpha": alpha,
        "inner_mean_ap": float(mean_cv[np.where(np.array(ALPHA_RATIOS) == ratio)[0][0]]),
        "inner_ap_by_alpha": mean_cv.tolist(),
        "inner_ap_by_alpha_and_fold": fold_cv.tolist(),
        "selected_features": int(len(detail)),
        "nonzero_coefficients": int(np.count_nonzero(np.abs(coefficients) > 1e-10)),
        "convergence_warning": nonconverged,
    }
    return prediction, component, baseline_test, feature_rows, tuning


def top10_metrics(y: np.ndarray, score: np.ndarray) -> dict[str, float]:
    k = int(np.ceil(len(y) * 0.1))
    chosen = np.argsort(-score, kind="stable")[:k]
    precision = float(y[chosen].mean())
    recall = float(y[chosen].sum() / y.sum()) if y.sum() else np.nan
    return {"top10_n": k, "top10_precision": precision,
            "top10_sensitivity": recall,
            "top10_enrichment": precision / y.mean() if y.mean() else np.nan}


def summarize(predictions: pd.DataFrame, output_dir: Path) -> None:
    positive = predictions.positive.to_numpy(dtype=bool)
    residual_positive = predictions.residual_positive.to_numpy(dtype=bool)
    score_columns = {
        "genomic_baseline": ("baseline_score", positive),
        "legacy_architecture_feature_only": ("legacy_architecture_feature_only_score", positive),
        "baseline_plus_residual_annotation": ("full_score", positive),
        "feature_only_residual": ("feature_only_residual_score", residual_positive),
    }
    for block in BLOCKS:
        score_columns[f"{block}_only"] = (f"{block}_component", positive)
    summary_rows = []
    chromosome_rows = []
    for name, (column, label) in score_columns.items():
        score = predictions[column].to_numpy(dtype=float)
        base = float(label.mean())
        overall = {"model": name, "n": len(label), "positives": int(label.sum()),
                   "prevalence": base, "auprc": average_precision_score(label, score),
                   **top10_metrics(label, score)}
        overall["auprc_delta_vs_prevalence"] = overall["auprc"] - base
        summary_rows.append(overall)
        for chrom, group in predictions.groupby("Chromosome", sort=False):
            idx = group.index.to_numpy()
            this_label = label[idx]
            this_score = score[idx]
            row = {"model": name, "chromosome": chrom, "n": len(idx),
                   "positives": int(this_label.sum()),
                   "prevalence": float(this_label.mean()),
                   "auprc": average_precision_score(this_label, this_score)
                   if this_label.any() else np.nan,
                   **top10_metrics(this_label, this_score)}
            row["auprc_delta_vs_prevalence"] = row["auprc"] - row["prevalence"]
            chromosome_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    by_chromosome = pd.DataFrame(chromosome_rows)
    baseline_ap = float(summary.loc[summary.model == "genomic_baseline", "auprc"].iloc[0])
    summary["auprc_delta_vs_genomic_baseline"] = summary.auprc - baseline_ap
    base_chr = by_chromosome.loc[by_chromosome.model == "genomic_baseline"].set_index("chromosome").auprc
    by_chromosome["auprc_delta_vs_genomic_baseline"] = [
        row.auprc - base_chr[row.chromosome] if row.model != "feature_only_residual" else np.nan
        for row in by_chromosome.itertuples()
    ]
    summary.to_csv(output_dir / "metrics_summary.csv", index=False)
    by_chromosome.to_csv(output_dir / "metrics_by_chromosome.csv", index=False)
    print(summary.to_string(index=False), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).parent / "data")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent / "results_primary")
    parser.add_argument("--max-folds", type=int, default=None)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = pd.read_csv(args.data_dir / "rows.csv")
    matrices = {block: sparse.load_npz(args.data_dir / f"X_{block}.npz").tocsr()
                for block in BLOCKS}
    feature_names = {block: pd.read_csv(args.data_dir / f"features_{block}.csv").feature.to_numpy()
                     for block in BLOCKS}
    assert len(rows) == 95750 and rows.bed_id.is_unique
    assert all(X.shape[0] == len(rows) for X in matrices.values())
    B = baseline_matrix(rows)
    y = rows.ZSTAT.to_numpy(dtype=float)
    chromosome = rows.Chromosome.to_numpy()
    chromosomes = sorted(set(chromosome), key=lambda name: int(name[3:]))
    if args.max_folds is not None:
        chromosomes = chromosomes[: args.max_folds]
    for fold_number, held_out in enumerate(chromosomes, start=1):
        fold_file = args.output_dir / f"fold_{held_out}.csv"
        if fold_file.exists():
            print(f"Skipping completed {held_out}", flush=True)
            continue
        train = np.flatnonzero(chromosome != held_out)
        test = np.flatnonzero(chromosome == held_out)
        print(f"Fold {fold_number}/{len(chromosomes)} {held_out}: train={len(train)}, test={len(test)}", flush=True)
        baseline_train, baseline_test = fit_baseline(B, y, train, test)
        frame = rows.loc[test, ["bed_id", "Chromosome", "ZSTAT"]].copy()
        frame["positive"] = y[test] > 1.96
        frame["baseline_score"] = baseline_test
        components = []
        features = []
        tuning_rows = []
        for block in BLOCKS:
            prediction, component, _, selected, tuning = fit_candidate(
                matrices, B, y, chromosome, train, test, (block,), "raw", 1.0, feature_names)
            frame[f"{block}_prediction"] = prediction
            frame[f"{block}_component"] = component
            components.append(component)
            features.extend({**item, "model": block, "holdout_chromosome": held_out}
                            for item in selected)
            tuning_rows.append({**tuning, "model": block, "holdout_chromosome": held_out})
            print(f"  {block}: CV AP={tuning['inner_mean_ap']:.4f}, nonzero={tuning['nonzero_coefficients']}", flush=True)
        frame["legacy_architecture_feature_only_score"] = np.abs(np.sum(components, axis=0))
        residual_prediction, _, _, selected, tuning = fit_candidate(
            matrices, B, y, chromosome, train, test, BLOCKS, "residual", 1.0, feature_names)
        frame["feature_only_residual_score"] = residual_prediction
        frame["full_score"] = baseline_test + residual_prediction
        heldout_residual = y[test] - baseline_test
        residual_threshold = float(np.quantile(y[train] - baseline_train, 0.9))
        frame["residual_zstat"] = heldout_residual
        frame["residual_positive"] = heldout_residual >= residual_threshold
        frame["residual_training_q90"] = residual_threshold
        features.extend({**item, "model": "combined_residual", "holdout_chromosome": held_out}
                        for item in selected)
        tuning_rows.append({**tuning, "model": "combined_residual", "holdout_chromosome": held_out})
        print(f"  combined residual: CV AP={tuning['inner_mean_ap']:.4f}, nonzero={tuning['nonzero_coefficients']}", flush=True)
        assert frame.bed_id.is_unique and frame.notna().all().all()
        frame.to_csv(fold_file, index=False)
        pd.DataFrame(features).to_csv(args.output_dir / f"features_{held_out}.csv", index=False)
        pd.DataFrame(tuning_rows).to_json(args.output_dir / f"tuning_{held_out}.json",
                                          orient="records", indent=2)
    fold_files = [args.output_dir / f"fold_{chrom}.csv" for chrom in chromosomes]
    if all(path.exists() for path in fold_files):
        all_predictions = pd.concat([pd.read_csv(path) for path in fold_files], ignore_index=True)
        if args.max_folds is None:
            assert len(all_predictions) == len(rows) and all_predictions.bed_id.is_unique
            assert set(all_predictions.bed_id) == set(rows.bed_id)
        all_predictions.to_csv(args.output_dir / "oof_predictions.csv", index=False)
        if args.max_folds is None:
            summarize(all_predictions, args.output_dir)
            all_features = pd.concat(
                [pd.read_csv(args.output_dir / f"features_{chrom}.csv") for chrom in chromosomes],
                ignore_index=True)
            all_features.to_csv(args.output_dir / "selected_features_by_fold.csv", index=False)
            stability = all_features.groupby(["model", "block", "feature"], as_index=False).agg(
                selected_folds=("holdout_chromosome", "nunique"),
                nonzero_folds=("nonzero", "sum"),
                median_abs_coefficient=("coefficient_scaled", lambda x: float(np.median(np.abs(x)))),
                median_screen_correlation=("screen_correlation", "median"),
            ).sort_values(["model", "selected_folds", "nonzero_folds"], ascending=[True, False, False])
            stability.to_csv(args.output_dir / "feature_stability.csv", index=False)
            (args.output_dir / "run_manifest.json").write_text(json.dumps({
                "data_rows": len(rows), "outer_folds": len(chromosomes),
                "positive_definition": "ZSTAT > 1.96",
                "residual_positive_definition": "heldout residual >= outer-training residual q90",
                "inner_folds": 3, "alpha_ratios": ALPHA_RATIOS,
                "baseline_ridge_alpha": BASELINE_RIDGE_ALPHA,
                "screen_top_k": TOP_K, "screen_min_nonzero": MIN_NONZERO,
                "legacy_architecture_score": "abs(sum of three per-block raw-X dot standardized-model coefficients)",
                "functional_score_excludes_genomic_baseline_at_inference": True,
            }, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
