"""Paired rank-enrichment analysis on the same coloc-evaluable junctions.

The primary contrast is the full-data rebuilt-Pfam score versus the absolute
beta of the minimum-P GTEx lead sQTL pair, with posterior > 0.9. The 0.8
posterior rule, alternative sQTL rankings, and LOCO scores are sensitivities.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, roc_auc_score

from compare_rebuilt_sqtl_magma import HERE, ROOT, TISSUE, reference_frames, score_frame


SEED = 171024
BOOTSTRAPS = 20_000
METRICS = {
    "rebuilt_full": ("all_junction_refit_new_pfam", False),
    "rebuilt_loco_fixed_alpha": ("loco_full_feature_fixed_alpha", False),
    "lead_sqtl_abs_beta": ("lead_abs_beta", False),
    "max_sqtl_abs_beta": ("sqtl_max_abs_slope", False),
    "sqtl_min_nominal_p": ("sqtl_min_pval_nominal", True),
}
CONTRASTS = (
    ("rebuilt_full", "lead_sqtl_abs_beta", "primary_full_vs_lead_beta"),
    ("rebuilt_loco_fixed_alpha", "lead_sqtl_abs_beta", "sensitivity_loco_vs_lead_beta"),
    ("rebuilt_full", "max_sqtl_abs_beta", "sensitivity_full_vs_max_beta"),
    ("rebuilt_full", "sqtl_min_nominal_p", "sensitivity_full_vs_min_p"),
)


def ordered(frame: pd.DataFrame, key: str, ascending: bool) -> pd.DataFrame:
    return frame.sort_values([key, "event_id"], ascending=[ascending, True], kind="stable")


def score_metrics(frame: pd.DataFrame, y: np.ndarray, key: str, ascending: bool) -> dict:
    values = frame[key].to_numpy(float)
    strength = -values if ascending else values
    order = ordered(frame, key, ascending)
    ranks = rankdata(-strength, method="average")
    hit_ranks = ranks[y]
    n, h = len(frame), int(y.sum())
    assert h > 0 and h < n
    # Area under the cumulative recall curve: E[1 - hit percentile rank].
    area = float(np.mean(1.0 - hit_ranks / n))
    out = {"auprc": float(average_precision_score(y, strength)),
           "auroc": float(roc_auc_score(y, strength)),
           "cumulative_area": area,
           "mean_hit_percentile": float(np.mean(hit_ranks / n)),
           "top_2pct_n": int(np.ceil(n * .02)),
           "top_10pct_n": int(np.ceil(n * .10)),
           "top_20pct_n": int(np.ceil(n * .20))}
    for frac, name in ((.02, "top_2pct"), (.10, "top_10pct"), (.20, "top_20pct")):
        k = int(np.ceil(n * frac))
        out[f"{name}_hits"] = int(order.head(k).positive.sum())
    return out


def block_signflip_and_ci(frame: pd.DataFrame, a: str, b: str,
                          y: np.ndarray, seed: int) -> tuple[float, float, float, int]:
    """Paired rank-position difference, with genomic-block sign-flip P.

    H0 is exchangeability of the two rankings' positive-event rank positions
    at the block level. This does not test model calibration or independence.
    """
    n = len(frame)
    ar = rankdata(-frame[a].to_numpy(float), method="average")
    bv = frame[b].to_numpy(float)
    br = rankdata(bv if b == "sqtl_min_pval_nominal" else -bv, method="average")
    pos = frame.loc[y, ["block"]].copy()
    pos["advantage"] = (br[y] - ar[y]) / n
    by_block = pos.groupby("block", dropna=False).agg(
        contribution=("advantage", "sum"), hits=("advantage", "size")
    )
    c = by_block.contribution.to_numpy(float)
    h = by_block.hits.to_numpy(int)
    observed = c.sum() / h.sum()
    if len(c) <= 20:
        codes = np.arange(1 << len(c), dtype=np.uint32)[:, None]
        signs = 2 * ((codes >> np.arange(len(c), dtype=np.uint32)) & 1).astype(np.int8) - 1
        null = signs @ c / h.sum()
        p = float(np.mean(null >= observed - 1e-12))
    else:
        rng = np.random.default_rng(seed)
        signs = rng.choice((-1, 1), size=(200_000, len(c)))
        null = signs @ c / h.sum()
        p = float((1 + np.count_nonzero(null >= observed - 1e-12)) / (len(null) + 1))
    rng = np.random.default_rng(seed + 1)
    draws = rng.integers(0, len(c), size=(BOOTSTRAPS, len(c)))
    boot = c[draws].sum(axis=1) / h[draws].sum(axis=1)
    low, high = np.quantile(boot, [.025, .975])
    assert np.isclose(observed,
                      score_metrics(frame, y, a, False)["cumulative_area"] -
                      score_metrics(frame, y, b, b == "sqtl_min_pval_nominal")["cumulative_area"])
    return p, float(low), float(high), len(c)


def curve_knots(frame: pd.DataFrame, key: str, ascending: bool,
                method: str, endpoint: str, threshold: float) -> list[dict]:
    o = ordered(frame, key, ascending)
    hit_rank = np.flatnonzero(o.positive.to_numpy(bool)) + 1
    n, h = len(o), len(hit_rank)
    knots = [{"endpoint": endpoint, "posterior_threshold": threshold,
              "method": method, "rank_fraction": 0.0, "cumulative_recall": 0.0}]
    for i, rank in enumerate(hit_rank, start=1):
        x = rank / n
        knots.append({"endpoint": endpoint, "posterior_threshold": threshold,
                      "method": method, "rank_fraction": x,
                      "cumulative_recall": (i - 1) / h})
        knots.append({"endpoint": endpoint, "posterior_threshold": threshold,
                      "method": method, "rank_fraction": x,
                      "cumulative_recall": i / h})
    knots.append({"endpoint": endpoint, "posterior_threshold": threshold,
                  "method": method, "rank_fraction": 1.0, "cumulative_recall": 1.0})
    return knots


def holm_two(p1: float, p2: float) -> tuple[float, float]:
    p = np.array([p1, p2])
    order = np.argsort(p)
    out = np.empty(2)
    out[order[0]] = min(1.0, 2 * p[order[0]])
    out[order[1]] = max(out[order[0]], p[order[1]])
    return float(out[0]), float(out[1])


def main() -> None:
    score = score_frame()
    lead = pd.read_csv(
        ROOT / "outputs/new_coloc_review/sqtl_baseline_sensitivity_20261017/junction_baseline_scores.csv",
        usecols=["event_id", "tissue_set", "lead_abs_beta", "sqtl_max_abs_slope",
                 "sqtl_min_pval_nominal"],
    )
    lead = lead.loc[lead.tissue_set.eq(TISSUE)].drop(columns="tissue_set")
    assert lead.event_id.is_unique
    rows, curves, details = [], [], []
    for endpoint_index, (endpoint, ref) in enumerate(reference_frames().items()):
        frame = ref.merge(score, on="event_id", how="inner", validate="one_to_one")
        frame = frame.merge(lead, on="event_id", how="inner", validate="one_to_one")
        assert frame.event_id.is_unique and frame.block.notna().all()
        assert frame[[v[0] for v in METRICS.values()]].notna().all().all()
        for threshold in (.9, .8):
            frame["positive"] = frame.posterior.gt(threshold)
            y = frame.positive.to_numpy(bool)
            base = {"endpoint": endpoint, "posterior_threshold": threshold,
                    "common_n": len(frame), "positive_n": int(y.sum()),
                    "prevalence": float(y.mean()),
                    "positive_blocks_n": int(frame.loc[y, "block"].nunique())}
            metric_rows = {}
            for method, (key, ascending) in METRICS.items():
                metric_rows[method] = score_metrics(frame, y, key, ascending)
                rows.append({**base, "kind": "ranking", "method": method,
                             **metric_rows[method]})
                if threshold == .9 and method in (
                    "rebuilt_full", "rebuilt_loco_fixed_alpha", "lead_sqtl_abs_beta"
                ):
                    curves.extend(curve_knots(frame, key, ascending, method, endpoint, threshold))
            for contrast_index, (model, comparator, label) in enumerate(CONTRASTS):
                p, low, high, nblocks = block_signflip_and_ci(
                    frame, METRICS[model][0], METRICS[comparator][0], y,
                    SEED + 100 * endpoint_index + 10 * contrast_index + int(threshold * 10),
                )
                rows.append({**base, "kind": "contrast", "method": label,
                             "model_method": model, "comparator_method": comparator,
                             "delta_cumulative_area": metric_rows[model]["cumulative_area"] -
                                                      metric_rows[comparator]["cumulative_area"],
                             "delta_auprc": metric_rows[model]["auprc"] -
                                            metric_rows[comparator]["auprc"],
                             "block_signflip_p_one_sided": p,
                             "block_bootstrap_ci_low": low,
                             "block_bootstrap_ci_high": high,
                             "tested_positive_blocks_n": nblocks})
        frame["endpoint"] = endpoint
        details.append(frame)
    out = pd.DataFrame(rows)
    primary = out.loc[(out.kind == "contrast") &
                      (out.method == "primary_full_vs_lead_beta") &
                      (out.posterior_threshold == .9)]
    assert len(primary) == 2
    q = holm_two(*primary.block_signflip_p_one_sided.to_list())
    for endpoint, adjusted in zip(primary.endpoint, q):
        mask = ((out.kind == "contrast") & (out.method == "primary_full_vs_lead_beta") &
                (out.posterior_threshold == .9) & (out.endpoint == endpoint))
        out.loc[mask, "holm_p_two_endpoints"] = adjusted
    out.to_csv(HERE / "rank_enrichment_model_vs_sqtl_stats.csv", index=False)
    pd.DataFrame(curves).to_csv(HERE / "rank_enrichment_model_vs_sqtl_curves.csv", index=False)
    pd.concat(details, ignore_index=True).to_csv(
        HERE / "rank_enrichment_model_vs_sqtl_common_events.csv", index=False
    )
    manifest = {
        "ranking_population": "same posterior-, model-, and GTEx-sQTL-evaluable junctions per endpoint",
        "model_score": "full-data GO + expression + rebuilt Pfam absolute score",
        "primary_sqtl_score": "absolute beta of minimum-nominal-P GTEx lead pair",
        "posterior_primary": ">0.9",
        "test": "paired one-sided sign flip of positive-junction rank-percentile advantage at regional-block level",
        "interval": f"{BOOTSTRAPS}-resample positive-block bootstrap 95% percentile interval",
        "primary_multiplicity": "Holm adjustment across two endpoints",
        "sensitivity": ["posterior >0.8", "max absolute sQTL slope", "minimum nominal sQTL P",
                        "fixed-alpha LOCO model score"],
        "seed": SEED,
    }
    (HERE / "rank_enrichment_model_vs_sqtl_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
