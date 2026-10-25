"""Evaluate existing LOCO scores at P < 5×10^-4 ... 5×10^-8.

This changes only the MAGMA-derived positive labels. No model is refit.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm

from sweep_loco_zstat_thresholds import HERE, SOURCES, metrics


CUTOFFS = (5e-4, 5e-5, 5e-6, 5e-7, 5e-8)
MODEL = "threequarter_training_chromosomes"


def main() -> None:
    source, pfam_col, combined_col = SOURCES[MODEL]
    frame = pd.read_csv(source)
    assert len(frame) == 95_750 and frame.bed_id.is_unique
    scores = {
        "expression": np.abs(frame.expression_component.to_numpy(float)),
        "GO": np.abs(frame.go_component.to_numpy(float)),
        "Pfam": np.abs(frame[pfam_col].to_numpy(float)),
        "combined": frame[combined_col].to_numpy(float),
    }
    pvalues = norm.sf(frame.ZSTAT.to_numpy(float))
    rows = []
    for cutoff in CUTOFFS:
        positive = pvalues < cutoff
        for block, score in scores.items():
            rows.append({
                "model": MODEL,
                "magma_p_cutoff": cutoff,
                "equivalent_z_cutoff": float(norm.isf(cutoff)),
                "block": block,
                **metrics(positive, score),
            })
    out = pd.DataFrame(rows)
    target = HERE / "loco_five_power_p_sweep.csv"
    out.to_csv(target, index=False)
    print(out.loc[:, ["magma_p_cutoff", "block", "positive_n", "positive_prevalence",
                      "auprc", "auprc_over_prevalence",
                      "sensitivity_at_10_expected_random_ties"]].to_string(index=False))
    print("Wrote", target)


if __name__ == "__main__":
    main()
