"""Rebuilt-Pfam score versus lead GTEx sQTL and conditional on MAGMA top 10%.

The primary score is the full-junction refit. LOCO scores are reported separately.
No score is refitted or tuned against the external posterior endpoints here.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from plot_new_coloc_four_panel_benchmark import contingency, matched_stats  # noqa: E402
from compare_sqtl_baseline_sensitivity import regional_swap  # noqa: E402


TISSUE = "Lung_plus_Whole_Blood_plus_Brain"
POSTERIOR_CUTOFF = 0.9


def score_frame() -> pd.DataFrame:
    base = pd.read_csv(HERE / "scores_comparison.csv",
                       usecols=["bed_id", "all_junction_refit_new_pfam"])
    z = pd.read_csv(ROOT / "sinform_loco_full/data/rows.csv", usecols=["bed_id", "ZSTAT"])
    assert len(base) == len(z) == 95750 and base.bed_id.equals(z.bed_id)
    base["ZSTAT"] = z.ZSTAT
    for label, path, column in (
        ("loco_full_feature_fixed_alpha", HERE / "loco_full_feature_fixed_alpha/oof_predictions.csv", "rebuilt_score"),
        ("loco_fold_selected_min30", HERE / "loco_rebuilt_pfam/oof_predictions.csv", "rebuilt_feature_only_score"),
        ("loco_threequarter_coverage", HERE / "loco_threequarter_chromosome_screen/oof_predictions.csv", "rebuilt_feature_only_score"),
    ):
        f = pd.read_csv(path, usecols=["bed_id", column])
        assert len(f) == len(base) and f.bed_id.equals(base.bed_id)
        base[label] = f[column].to_numpy(float)
    k = math.ceil(len(base) * .1)
    base["magma_top10"] = False
    base.loc[base.ZSTAT.nlargest(k, keep="first").index, "magma_top10"] = True
    for name in ["all_junction_refit_new_pfam", "loco_full_feature_fixed_alpha",
                 "loco_fold_selected_min30", "loco_threequarter_coverage"]:
        base[f"top10_{name}"] = False
        base.loc[base[name].nlargest(k, keep="first").index, f"top10_{name}"] = True
    return base.rename(columns={"bed_id": "event_id"})


def reference_frames() -> dict[str, pd.DataFrame]:
    ltrc = pd.read_csv(ROOT / "outputs/new_coloc_review/moloc_external/ranking_comparison/LTRC_GWAS_significant_matched_events.csv",
                       usecols=["event_id", "coloc_ppas", "block", "gene", "gwas_p"])
    ltrc = ltrc.loc[ltrc.gwas_p.lt(5e-8) & ltrc.coloc_ppas.notna()].copy()
    ltrc = ltrc.rename(columns={"coloc_ppas": "posterior", "gene": "gene_id"})
    fev = pd.read_csv(ROOT / "outputs/new_coloc_review/FEV1FVC_four_panel_matched_events.csv",
                      usecols=["event_id", "PP.H4.abf", "locus_id", "gene_id"])
    blocks = pd.read_csv(ROOT / "outputs/new_coloc_review/sqtl_baseline_sensitivity_20261017/FEV1FVC_regional_blocks.csv",
                         usecols=["locus_id", "block"])
    fev = fev.merge(blocks, on="locus_id", how="left", validate="many_to_one")
    fev = fev.rename(columns={"PP.H4.abf": "posterior"})
    assert fev.block.notna().all()
    assert ltrc.event_id.is_unique and fev.event_id.is_unique
    return {"ICGC_COPD_LTRC": ltrc, "FEV1FVC_GTEx": fev}


def main() -> None:
    score = score_frame()
    lead = pd.read_csv(ROOT / "outputs/new_coloc_review/sqtl_baseline_sensitivity_20261017/junction_baseline_scores.csv",
                       usecols=["event_id", "tissue_set", "lead_abs_beta", "lead_p", "lead_tissue"])
    lead = lead.loc[lead.tissue_set.eq(TISSUE)].drop(columns="tissue_set")
    assert lead.event_id.is_unique and lead.lead_abs_beta.notna().all()
    variants = ["all_junction_refit_new_pfam", "loco_full_feature_fixed_alpha",
                "loco_fold_selected_min30", "loco_threequarter_coverage"]
    ab, cd, matched_rows, parent_rows = [], [], [], []
    manifest = {"n_scored": len(score), "top10_k": math.ceil(len(score) * .1),
                "posterior_rule": ">0.9", "baseline": "GTEx min nominal P pair |beta| across lung, blood, 13 brain tissues",
                "variants": variants, "endpoints": {}}
    for endpoint, ref in reference_frames().items():
        base = ref.merge(score, on="event_id", how="inner", validate="one_to_one")
        common = base.merge(lead, on="event_id", how="inner", validate="one_to_one")
        parent = base.loc[base.magma_top10].copy()
        assert base.posterior.notna().all() and base.block.notna().all()
        assert common.event_id.is_unique and parent.event_id.is_unique
        manifest["endpoints"][endpoint] = {"scored_tested_n": len(base),
                                            "with_lead_beta_n": len(common),
                                            "magma_top10_evaluable_n": len(parent),
                                            "posterior_at_0p9_n": int(base.posterior.eq(.9).sum())}
        for variant in variants:
            selected_col = f"top10_{variant}"
            k = int(common[selected_col].sum())
            comparator_ids = set(common.sort_values(["lead_abs_beta", "event_id"],
                                                    ascending=[False, True], kind="stable").head(k).event_id)
            shared = common.copy()
            shared["model_selected"] = shared[selected_col].astype(bool)
            shared["sqtl_selected"] = shared.event_id.isin(comparator_ids)
            shared["positive"] = shared.posterior.gt(POSTERIOR_CUTOFF)
            assert shared.model_selected.sum() == shared.sqtl_selected.sum() == k
            model_only = int((shared.positive & shared.model_selected & ~shared.sqtl_selected).sum())
            sqtl_only = int((shared.positive & ~shared.model_selected & shared.sqtl_selected).sum())
            odds, low, high, p = matched_stats(model_only, sqtl_only)
            region_p, discordant_regions = regional_swap(shared, shared.model_selected,
                                                           shared.sqtl_selected, shared.positive)
            ab.append({"endpoint": endpoint, "score_variant": variant,
                       "common_n": len(shared), "selected_n_each": k,
                       "positive_n": int(shared.positive.sum()),
                       "model_hits": int((shared.model_selected & shared.positive).sum()),
                       "sqtl_hits": int((shared.sqtl_selected & shared.positive).sum()),
                       "model_only_hits": model_only, "sqtl_only_hits": sqtl_only,
                       "both_hits": int((shared.model_selected & shared.sqtl_selected & shared.positive).sum()),
                       "matched_or": odds, "ci_low": low, "ci_high": high,
                       "paired_event_p_one_sided": p, "region_swap_p_one_sided": region_p,
                       "discordant_regions": discordant_regions})
            shared["endpoint"] = endpoint
            shared["score_variant"] = variant
            matched_rows.append(shared)

            within = parent.copy()
            within["model_selected"] = within[selected_col].astype(bool)
            within["positive"] = within.posterior.gt(POSTERIOR_CUTOFF)
            a = int((within.model_selected & within.positive).sum())
            b = int((within.model_selected & ~within.positive).sum())
            c = int((~within.model_selected & within.positive).sum())
            d = int((~within.model_selected & ~within.positive).sum())
            s = contingency(a, b, c, d)
            # A region selected by a model has a positive outcome only if one
            # of its selected junctions has high shared posterior.
            within["selected_positive"] = within.model_selected & within.positive
            group = within.groupby("block").agg(selected=("model_selected", "any"),
                                                selected_positive=("selected_positive", "any"),
                                                any_positive=("positive", "any"))
            group["positive"] = np.where(group.selected, group.selected_positive, group.any_positive)
            ra = int((group.selected & group.positive).sum())
            rb = int((group.selected & ~group.positive).sum())
            rc = int((~group.selected & group.positive).sum())
            rd = int((~group.selected & ~group.positive).sum())
            region = contingency(ra, rb, rc, rd) if (ra + rb and rc + rd) else {}
            genes = within.groupby("gene_id", dropna=True).agg(
                selected=("model_selected", "any"),
                selected_positive=("selected_positive", "any"),
                any_positive=("positive", "any"),
            )
            genes["positive"] = np.where(genes.selected, genes.selected_positive,
                                          genes.any_positive)
            ga = int((genes.selected & genes.positive).sum())
            gb = int((genes.selected & ~genes.positive).sum())
            gc = int((~genes.selected & genes.positive).sum())
            gd = int((~genes.selected & ~genes.positive).sum())
            gene = contingency(ga, gb, gc, gd) if (ga + gb and gc + gd) else {}
            cd.append({"endpoint": endpoint, "score_variant": variant,
                       "magma_parent_n": len(within), **s,
                       "selected_positive_regions": int(group.loc[group.selected & group.positive].shape[0]),
                       "other_positive_regions": int(group.loc[~group.selected & group.positive].shape[0]),
                       "region_selected_n": ra + rb, "region_selected_hits": ra,
                       "region_other_n": rc + rd, "region_other_hits": rc,
                       "region_or": region.get("fisher_or", math.nan),
                       "region_p_one_sided": region.get("fisher_p_greater", math.nan),
                       "gene_selected_n": ga + gb, "gene_selected_hits": ga,
                       "gene_other_n": gc + gd, "gene_other_hits": gc,
                       "gene_or": gene.get("fisher_or", math.nan),
                       "gene_p_one_sided": gene.get("fisher_p_greater", math.nan)})
            within["endpoint"] = endpoint
            within["score_variant"] = variant
            parent_rows.append(within)
    pd.DataFrame(ab).to_csv(HERE / "rebuilt_vs_lead_sqtl_pp_gt_0p9.csv", index=False)
    pd.DataFrame(cd).to_csv(HERE / "rebuilt_within_magma_pp_gt_0p9.csv", index=False)
    pd.concat(matched_rows, ignore_index=True).to_csv(HERE / "rebuilt_vs_lead_sqtl_matched_events.csv", index=False)
    pd.concat(parent_rows, ignore_index=True).to_csv(HERE / "rebuilt_within_magma_events.csv", index=False)
    (HERE / "rebuilt_sqtl_magma_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("A/B\n", pd.DataFrame(ab).to_string(index=False))
    print("C/D\n", pd.DataFrame(cd).to_string(index=False))
    print("MANIFEST\n", json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
