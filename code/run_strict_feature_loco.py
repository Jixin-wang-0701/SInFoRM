"""Exploratory fold-local screening for GO and rebuilt Pfam LOCO.

Defaults require >=50 nonzero training rows and >=5 training chromosomes.
Optional arguments allow a fraction of the available training chromosomes;
all screening and fitting remain inside each outer/inner training split.
Expression is unchanged and its saved independent LOCO component is reused.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = HERE / "loco_strict_feature_screen"
DATA = ROOT / "sinform_loco_full/data"
sys.path.insert(0, str(ROOT / "sinform_loco_full"))
import run_loco as base  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chromosomes", nargs="*")
    parser.add_argument("--aggregate-only", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=OUT)
    parser.add_argument("--min-nonzero", type=int, default=50)
    parser.add_argument("--min-active-chromosomes", type=int, default=5)
    parser.add_argument("--min-active-fraction", type=float)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    if args.min_nonzero < 1 or args.min_active_chromosomes < 1:
        parser.error("Minimum nonzero and chromosome counts must be positive")
    if args.min_active_fraction is not None and not (0 < args.min_active_fraction <= 1):
        parser.error("--min-active-fraction must be in (0, 1]")
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "manifest.json"
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        requested = {
            "go_pfam_min_nonzero_training_rows": args.min_nonzero,
            "go_pfam_min_active_training_chromosomes": (
                args.min_active_chromosomes if args.min_active_fraction is None else None),
            "go_pfam_min_active_training_chromosome_fraction": args.min_active_fraction,
        }
        if any(previous.get(key) != value for key, value in requested.items()):
            parser.error("Output directory already contains folds from different screening settings")
    rows = pd.read_csv(DATA / "rows.csv", usecols=["bed_id", "Chromosome", "ZSTAT"])
    old_oof = pd.read_csv(ROOT / "sinform_loco_full/results_primary/oof_predictions.csv",
                          usecols=["bed_id", "expression_component"])
    assert len(rows) == len(old_oof) == 95750
    assert rows.bed_id.is_unique and old_oof.bed_id.is_unique
    chromosome = rows.Chromosome.to_numpy()
    chromosomes = sorted(set(chromosome), key=lambda name: int(name[3:]))
    chrom_index = {name: i for i, name in enumerate(chromosomes)}
    y = rows.ZSTAT.to_numpy(float)
    matrices = {
        "go": sparse.load_npz(DATA / "X_go.npz").tocsr(),
        "pfam": sparse.csr_matrix(mmread(HERE / "X_pfam_rebuilt.mtx")),
    }
    names = {
        "go": pd.read_csv(DATA / "features_go.csv").feature.to_numpy(),
        "pfam": pd.read_csv(HERE / "pfam_feature_map.csv").feature.to_numpy(),
    }
    by_chromosome = {
        block: np.stack([
            np.asarray(X[chromosome == name].getnnz(axis=0)).ravel()
            for name in chromosomes
        ])
        for block, X in matrices.items()
    }

    def strict_screen_block(X: sparse.csr_matrix, train: np.ndarray,
                            target: np.ndarray, block: str):
        Xtr = X[train]
        n = len(train)
        sx = np.asarray(Xtr.sum(axis=0)).ravel().astype(float)
        sx2 = np.asarray(Xtr.power(2).sum(axis=0)).ravel().astype(float)
        sy = float(target.sum())
        vary = float(target @ target - sy * sy / n)
        sxy = np.asarray(Xtr.T @ target).ravel().astype(float)
        varx = sx2 - sx * sx / n
        correlation = (sxy - sx * sy / n) / np.sqrt(np.maximum(varx * vary, 1e-20))
        nonzero = np.asarray(Xtr.getnnz(axis=0)).ravel()
        training_chromosomes = [chrom_index[name] for name in set(chromosome[train])]
        breadth = (by_chromosome[block][training_chromosomes] > 0).sum(axis=0)
        required_chromosomes = (math.ceil(len(training_chromosomes) * args.min_active_fraction)
                                if args.min_active_fraction is not None
                                else args.min_active_chromosomes)
        eligible = (np.isfinite(correlation) & (varx > 1e-12)
                    & (nonzero >= args.min_nonzero)
                    & (nonzero <= n - args.min_nonzero)
                    & (breadth >= required_chromosomes))
        indices = np.flatnonzero(eligible)
        indices = indices[np.argsort(-np.abs(correlation[indices]), kind="stable")]
        indices = indices[:100]
        return indices, correlation[indices]

    # run_loco.fit_candidate performs fold-local screening in both inner CV
    # and outer fitting. Replace only its screen function for this process.
    base.screen_block = strict_screen_block
    chosen = args.chromosomes if args.chromosomes is not None else chromosomes
    assert set(chosen) <= set(chromosomes)
    for held_out in ([] if args.aggregate_only else chosen):
        fold_file = output_dir / f"fold_{held_out}.csv"
        if fold_file.exists():
            print(f"Reusing {held_out}", flush=True)
            continue
        train = np.flatnonzero(chromosome != held_out)
        test = np.flatnonzero(chromosome == held_out)
        fold = rows.loc[test, ["bed_id", "Chromosome", "ZSTAT"]].copy()
        diagnostics = {}
        selected_rows = []
        for block in ("go", "pfam"):
            _, component, _, selected, tuning = base.fit_candidate(
                {block: matrices[block]}, np.empty((len(rows), 0)), y, chromosome,
                train, test, (block,), "raw", 1.0, {block: names[block]}
            )
            column = "rebuilt_pfam_component" if block == "pfam" else "go_component"
            fold[column] = component
            diagnostics[block] = tuning
            selected_rows.extend({**item, "holdout_chromosome": held_out}
                                 for item in selected)
        fold.to_csv(fold_file, index=False)
        pd.DataFrame(selected_rows).to_csv(output_dir / f"features_{held_out}.csv", index=False)
        (output_dir / f"tuning_{held_out}.json").write_text(
            json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8"
        )
        print(held_out, len(test),
              {block: {"selected": value["selected_features"],
                       "nonzero_coef": value["nonzero_coefficients"],
                       "alpha_ratio": value["alpha_ratio"]}
               for block, value in diagnostics.items()}, flush=True)

    if not all((output_dir / f"fold_{name}.csv").exists() for name in chromosomes):
        print("Awaiting remaining chromosome folds", flush=True)
        return
    strict = pd.concat([pd.read_csv(output_dir / f"fold_{name}.csv")
                        for name in chromosomes], ignore_index=True)
    assert len(strict) == len(rows) and strict.bed_id.is_unique
    output = strict.merge(old_oof, on="bed_id", validate="one_to_one")
    output["rebuilt_feature_only_score"] = np.abs(
        output.expression_component + output.go_component
        + output.rebuilt_pfam_component
    )
    output.to_csv(output_dir / "oof_predictions.csv", index=False)
    (output_dir / "manifest.json").write_text(json.dumps({
        "outer_folds": 22, "inner_folds": 3,
        "go_pfam_min_nonzero_training_rows": args.min_nonzero,
        "go_pfam_min_active_training_chromosomes": (
            args.min_active_chromosomes if args.min_active_fraction is None else None),
        "go_pfam_min_active_training_chromosome_fraction": args.min_active_fraction,
        "go_pfam_top_k_by_training_absolute_correlation": 100,
        "expression_component": "reused from independent original nested LOCO; unchanged filter",
        "alpha_ratios": base.ALPHA_RATIOS,
        "status": "exploratory post-hoc threshold sensitivity; not independent model selection",
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(output)} out-of-fold rows", flush=True)


if __name__ == "__main__":
    main()
