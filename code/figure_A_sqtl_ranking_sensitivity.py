"""Post-hoc comparison of four sQTL ranking rules on fixed coloc common sets.

This does not alter the prespecified lead-|beta| six-panel figure. The same
posterior > 0.9 hits and positive genomic blocks are used throughout.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from compare_rebuilt_sqtl_magma import HERE, ROOT, TISSUE
from rank_enrichment_model_vs_sqtl import (
    block_signflip_and_ci,
    curve_knots,
    score_metrics,
)


METHODS = {
    "lead_minP_abs_beta": ("lead_abs_beta", False),
    "max_abs_beta": ("sqtl_max_abs_slope", False),
    "min_nominal_P": ("sqtl_min_pval_nominal", True),
    "max_abs_beta_over_SE": ("sqtl_max_abs_z", False),
}


def holm(p_values: pd.Series) -> np.ndarray:
    values = p_values.to_numpy(float)
    order = np.argsort(values)
    adjusted_sorted = np.maximum.accumulate(
        values[order] * np.arange(len(values), 0, -1)
    ).clip(max=1)
    adjusted = np.empty_like(values)
    adjusted[order] = adjusted_sorted
    return adjusted


def main() -> None:
    common = pd.read_csv(HERE / "rank_enrichment_model_vs_sqtl_common_events.csv")
    extra = pd.read_csv(
        ROOT / "outputs/new_coloc_review/sqtl_baseline_sensitivity_20261017/junction_baseline_scores.csv",
        usecols=["event_id", "tissue_set", "sqtl_max_abs_z"],
    )
    extra = extra.loc[extra.tissue_set.eq(TISSUE)].drop(columns="tissue_set")
    assert extra.event_id.is_unique
    common = common.merge(extra, on="event_id", how="left", validate="many_to_one")
    assert common.sqtl_max_abs_z.notna().all()
    rows, curves = [], []
    for endpoint_index, (endpoint, frame) in enumerate(common.groupby("endpoint", sort=True)):
        frame = frame.copy()
        frame["positive"] = frame.posterior.gt(.9)
        y = frame.positive.to_numpy(bool)
        model = score_metrics(frame, y, "all_junction_refit_new_pfam", False)
        for method_index, (method, (column, ascending)) in enumerate(METHODS.items()):
            base = score_metrics(frame, y, column, ascending)
            p, ci_low, ci_high, nblocks = block_signflip_and_ci(
                frame, "all_junction_refit_new_pfam", column, y,
                171024 + 100 * endpoint_index + 10 * method_index,
            )
            rows.append({
                "endpoint": endpoint,
                "sqtl_ranking": method,
                "common_junctions": len(frame),
                "coloc_hits_PP_gt_0p9": int(y.sum()),
                "positive_genomic_blocks": nblocks,
                "sinform_AUPRC": model["auprc"],
                "sqtl_AUPRC": base["auprc"],
                "sinform_curve_area": model["cumulative_area"],
                "sqtl_curve_area": base["cumulative_area"],
                "curve_area_difference": model["cumulative_area"] - base["cumulative_area"],
                "region_signflip_P_one_sided": p,
                "block_bootstrap_CI_low": ci_low,
                "block_bootstrap_CI_high": ci_high,
            })
            curves.extend(curve_knots(frame, column, ascending, method, endpoint, .9))
    out = pd.DataFrame(rows)
    for endpoint in out.endpoint.unique():
        idx = out.endpoint.eq(endpoint)
        out.loc[idx, "holm_P_four_rankings_within_endpoint"] = holm(
            out.loc[idx, "region_signflip_P_one_sided"]
        )
    out["holm_P_eight_comparisons"] = holm(out.region_signflip_P_one_sided)
    out.to_csv(HERE / "figure_A_sqtl_ranking_sensitivity.csv", index=False)
    pd.DataFrame(curves).to_csv(HERE / "figure_A_sqtl_ranking_sensitivity_curves.csv", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
