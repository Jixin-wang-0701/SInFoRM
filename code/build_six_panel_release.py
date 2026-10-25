"""Collect the full-fit and three-quarter-coverage LOCO model behind Figure X.

Run from the project with: python outputs/pfam_rebuild_npnt/build_six_panel_release.py
Existing fitted models are not retrained. All copied artifacts are hash-listed.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


SOURCE = Path(__file__).resolve()
ROOT = SOURCE.parents[4] if SOURCE.parent.name == "code" else SOURCE.parents[2]
HERE = ROOT / "outputs" / "pfam_rebuild_npnt"
RELEASE = HERE / "six_panel_release_20261025"
DATA = ROOT / "sinform_loco_full" / "data"
LOCO = HERE / "loco_threequarter_chromosome_screen"
BASE = ROOT / "sinform_loco_full" / "results_primary"


def copy(source: Path, destination: str) -> None:
    assert source.is_file(), source
    target = RELEASE / destination
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def scores() -> None:
    rows = pd.read_csv(DATA / "rows.csv", usecols=["bed_id", "Chromosome", "ZSTAT"])
    assert len(rows) == 95_750 and rows.bed_id.is_unique
    matrices = {
        "expression": sparse.load_npz(DATA / "X_expression.npz").tocsr(),
        "go": sparse.load_npz(DATA / "X_go.npz").tocsr(),
        "pfam": sparse.csr_matrix(mmread(HERE / "X_pfam_rebuilt.mtx")),
    }
    features = {
        "expression": pd.read_csv(DATA / "features_expression.csv").feature,
        "go": pd.read_csv(DATA / "features_go.csv").feature,
        "pfam": pd.read_csv(HERE / "pfam_feature_map.csv").feature,
    }
    components = {}
    coef_frames = []
    for block in matrices:
        X = matrices[block]
        coef = np.load(HERE / f"coef_{'pfam_new' if block == 'pfam' else block}.npy")
        assert X.shape == (len(rows), len(features[block]))
        assert len(coef) == X.shape[1]
        components[block] = np.asarray(X @ coef).ravel()
        coef_frames.append(pd.DataFrame({
            "block": block,
            "feature_index_zero_based": np.arange(len(coef)),
            "feature": features[block].to_numpy(),
            "coefficient_scaled_fit_applied_to_raw_matrix": coef,
            "nonzero": coef != 0,
        }))
    pd.concat(coef_frames, ignore_index=True).to_csv(
        RELEASE / "coefficients" / "full_fit_coefficients.csv", index=False)

    full = rows.copy()
    for block, value in components.items():
        full[f"{block}_component"] = value
        full[f"{block}_abs_score"] = np.abs(value)
    full["signed_sum"] = sum(components.values())
    full["feature_only_abs_score"] = np.abs(full.signed_sum)
    full["rank_desc"] = full.feature_only_abs_score.rank(
        ascending=False, method="first").astype(int)
    prior = pd.read_csv(HERE / "scores_comparison.csv",
                        usecols=["bed_id", "all_junction_refit_new_pfam"])
    check = full[["bed_id", "feature_only_abs_score"]].merge(prior, on="bed_id",
        validate="one_to_one", sort=False)
    full_diff = float(np.max(np.abs(check.feature_only_abs_score - check.all_junction_refit_new_pfam)))
    assert full_diff < 1e-10, full_diff

    held = pd.read_csv(LOCO / "oof_predictions.csv")
    assert len(held) == len(rows) and held.bed_id.is_unique
    held = rows[["bed_id"]].merge(held, on="bed_id", validate="one_to_one", sort=False)
    assert held.Chromosome.equals(rows.Chromosome)
    assert np.allclose(held.ZSTAT, rows.ZSTAT)
    out = held[["bed_id", "Chromosome", "ZSTAT"]].copy()
    for name, source in (("expression", "expression_component"),
                         ("go", "go_component"),
                         ("pfam", "rebuilt_pfam_component")):
        out[f"{name}_component"] = held[source]
        out[f"{name}_abs_score"] = held[source].abs()
    out["signed_sum"] = (out.expression_component + out.go_component
                         + out.pfam_component)
    out["feature_only_abs_score"] = out.signed_sum.abs()
    out["rank_desc"] = out.feature_only_abs_score.rank(
        ascending=False, method="first").astype(int)
    loco_diff = float(np.max(np.abs(out.feature_only_abs_score - held.rebuilt_feature_only_score)))
    assert loco_diff < 1e-10, loco_diff

    score_dir = RELEASE / "scores"
    score_dir.mkdir(parents=True, exist_ok=True)
    columns = ["bed_id", "Chromosome", "ZSTAT",
               "expression_component", "go_component", "pfam_component",
               "expression_abs_score", "go_abs_score", "pfam_abs_score",
               "signed_sum", "feature_only_abs_score", "rank_desc"]
    full[columns].to_csv(score_dir / "full_fit_scores.csv", index=False)
    out[columns].to_csv(score_dir / "loco_scores.csv", index=False)
    sparse.save_npz(RELEASE / "matrices" / "X_pfam_rebuilt.npz", matrices["pfam"])
    (RELEASE / "score_checks.json").write_text(json.dumps({
        "rows": len(rows),
        "full_fit_max_abs_difference_vs_existing_score": full_diff,
        "full_fit_score_reference": "outputs/pfam_rebuild_npnt/scores_comparison.csv",
        "loco_max_abs_difference_vs_existing_score": loco_diff,
        "matrix_shapes": {k: list(x.shape) for k, x in matrices.items()},
        "matrix_nnz": {k: int(x.nnz) for k, x in matrices.items()},
        "score_definition": "abs(expression_component + go_component + pfam_component)",
        "component_definition": "raw sparse feature matrix dot fitted standardized-model coefficients, no intercept",
    }, indent=2) + "\n", encoding="utf-8")


def fold_coefficients() -> None:
    frames = []
    for number in range(1, 23):
        chromosome = f"chr{number}"
        selected = pd.read_csv(LOCO / f"features_{chromosome}.csv")
        assert set(selected.block) <= {"go", "pfam"}
        base = pd.read_csv(BASE / f"features_{chromosome}.csv")
        expression = base.loc[base.block.eq("expression") & base.model.eq("expression"),
                              ["block", "feature", "screen_correlation",
                               "coefficient_scaled", "nonzero", "holdout_chromosome"]]
        assert len(expression) > 0
        frames.extend([expression, selected])
    table = pd.concat(frames, ignore_index=True)
    assert set(table.holdout_chromosome) == {f"chr{i}" for i in range(1, 23)}
    table.to_csv(RELEASE / "coefficients" / "loco_fold_coefficients.csv", index=False)


def files() -> None:
    for name in ["rows.csv", "features_expression.csv", "features_go.csv",
                 "X_expression.npz", "X_go.npz"]:
        copy(DATA / name, f"matrices/{name}")
    copy(HERE / "pfam_feature_map.csv", "matrices/features_pfam_rebuilt.csv")
    copy(HERE / "X_pfam_rebuilt.mtx", "matrices/X_pfam_rebuilt.mtx")
    copy(ROOT / "df_pfam_rebuilt.rds", "matrices/df_pfam_rebuilt.rds")
    for block in ("expression", "go", "pfam_new"):
        copy(HERE / f"coef_{block}.npy", f"coefficients/coef_{block}.npy")
    copy(HERE / "result.json", "full_fit/result.json")
    copy(HERE / "scores_comparison.csv", "full_fit/scores_comparison.csv")
    copy(LOCO / "manifest.json", "loco/manifest.json")
    copy(LOCO / "oof_predictions.csv", "loco/oof_predictions.csv")
    copy(BASE / "oof_predictions.csv", "loco/expression_source_oof_predictions.csv")
    for number in range(1, 23):
        chrom = f"chr{number}"
        for prefix in ("features", "tuning", "fold"):
            extension = "json" if prefix == "tuning" else "csv"
            copy(LOCO / f"{prefix}_{chrom}.{extension}",
                 f"loco/threequarter_folds/{prefix}_{chrom}.{extension}")
        for prefix in ("features", "tuning"):
            extension = "json" if prefix == "tuning" else "csv"
            copy(BASE / f"{prefix}_{chrom}.{extension}",
                 f"loco/expression_source_folds/{prefix}_{chrom}.{extension}")
    for name in ("rank_enrichment_model_vs_sqtl_stats.csv",
                 "rank_enrichment_model_vs_sqtl_curves.csv",
                 "rank_enrichment_model_vs_sqtl_common_events.csv",
                 "rank_enrichment_model_vs_sqtl_manifest.json",
                 "global_top2_figure_stats.csv", "global_top2_figure_events.csv",
                 "global_top2_figure_manifest.json", "loco_five_power_p_sweep.csv",
                 "figure_A_sqtl_ranking_sensitivity.csv",
                 "figure_A_sqtl_ranking_sensitivity_curves.csv",
                 "figure_A_sqtl_ranking_sensitivity.md"):
        copy(HERE / name, f"figure_data/{name}")
    copy(ROOT / "figures" / "rebuilt_pfam_six_panel.png", "figure/six_panel.png")
    copy(ROOT / "figures" / "pptx" / "rebuilt_pfam_six_panel_sqtl_magma_loco_editable.pptx",
         "figure/six_panel_editable.pptx")
    copy(ROOT / "figures" / "rebuilt_pfam_six_panel_caption.md", "figure/detailed_caption.md")
    for name in ("evaluate_npnt.py", "run_strict_feature_loco.py",
                 "export_rebuilt_pfam.R", "rank_enrichment_model_vs_sqtl.py",
                 "compare_global_top2_figure.py", "compare_rebuilt_sqtl_magma.py",
                 "loco_five_power_p_sweep.py", "sweep_loco_zstat_thresholds.py",
                 "figure_A_sqtl_ranking_sensitivity.py", "build_six_panel_release.py"):
        copy(HERE / name, f"code/{name}")
    copy(ROOT / "sinform_loco_full" / "run_loco.py", "code/run_loco.py")
    copy(ROOT / "sinform_loco_full" / "prepare_full_universe.R",
         "code/prepare_full_universe.R")
    copy(HERE / "figure_build" / "make_rebuilt_pfam_six_panel.mjs",
         "code/make_rebuilt_pfam_six_panel.mjs")


def manifest() -> None:
    records = []
    for file in sorted(RELEASE.rglob("*")):
        if not file.is_file() or file.name == "file_manifest.csv":
            continue
        records.append({
            "path_in_release": file.relative_to(RELEASE).as_posix(),
            "bytes": file.stat().st_size,
            "sha256": hashlib.sha256(file.read_bytes()).hexdigest(),
        })
    pd.DataFrame(records).to_csv(RELEASE / "file_manifest.csv", index=False)
    print("Release:", RELEASE, "files:", len(records))


def main() -> None:
    for name in ("coefficients", "matrices", "scores"):
        (RELEASE / name).mkdir(parents=True, exist_ok=True)
    files()
    scores()
    fold_coefficients()
    manifest()


if __name__ == "__main__":
    main()
