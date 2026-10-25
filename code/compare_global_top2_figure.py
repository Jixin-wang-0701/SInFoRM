"""Global top-2% comparison for the four-panel sensitivity figure.

Model and MAGMA are ranked over all 95,750 scoreable junctions. The lead-sQTL
reference is ranked over the 9,685 scoreable junctions with an available GTEx
lead-pair beta. These independently ranked global sets are intersected with
each external colocalization panel only after selection.
"""

from __future__ import annotations

import json
import math

import pandas as pd

from compare_rebuilt_sqtl_magma import (
    HERE,
    ROOT,
    TISSUE,
    contingency,
    reference_frames,
    score_frame,
)


def grouped_stats(parent: pd.DataFrame, level: str) -> dict:
    group = parent.groupby(level, dropna=True).agg(
        selected=("model_selected", "any"),
        selected_positive=("selected_positive", "any"),
        any_positive=("positive", "any"),
    )
    positive = group.selected_positive.where(group.selected, group.any_positive)
    a = int((group.selected & positive).sum())
    b = int((group.selected & ~positive).sum())
    c = int((~group.selected & positive).sum())
    d = int((~group.selected & ~positive).sum())
    return contingency(a, b, c, d) if a + b and c + d else {}


def main() -> None:
    score = score_frame()
    model_k = math.ceil(len(score) * .02)
    model_ids = set(score.nlargest(model_k, "all_junction_refit_new_pfam", keep="first").event_id)
    magma_ids = set(score.nlargest(model_k, "ZSTAT", keep="first").event_id)
    assert score.all_junction_refit_new_pfam.ge(0).all()

    lead = pd.read_csv(
        ROOT / "outputs/new_coloc_review/sqtl_baseline_sensitivity_20261017/junction_baseline_scores.csv",
        usecols=["event_id", "tissue_set", "lead_abs_beta"],
    )
    lead = lead.loc[lead.tissue_set.eq(TISSUE), ["event_id", "lead_abs_beta"]]
    lead = lead.loc[lead.event_id.isin(score.event_id)].copy()
    assert lead.event_id.is_unique and lead.lead_abs_beta.notna().all()
    lead_k = math.ceil(len(lead) * .02)
    lead_ids = set(lead.sort_values(
        ["lead_abs_beta", "event_id"], ascending=[False, True], kind="stable"
    ).head(lead_k).event_id)

    rows, detail = [], []
    manifest = {"model_and_magma_global_universe_n": len(score),
                "model_and_magma_global_top2_n": model_k,
                "lead_sqtl_global_universe_n": len(lead),
                "lead_sqtl_global_top2_n": lead_k,
                "model_selection": "global top 2% of full-junction refit GO + expression + rebuilt Pfam absolute score",
                "magma_selection": "global top 2% MAGMA ZSTAT on same 95,750 junctions",
                "lead_selection": "global top 2% |beta| of minimum nominal-P GTEx lead pair over scoreable junctions with lead beta",
                "posterior_rules": [">0.9", ">0.8"],
                "endpoints": {}}
    for endpoint, ref in reference_frames().items():
        tested = ref.merge(score, on="event_id", how="inner", validate="one_to_one")
        common = tested.merge(lead, on="event_id", how="inner", validate="one_to_one")
        parent = tested.loc[tested.event_id.isin(magma_ids)].copy()
        manifest["endpoints"][endpoint] = {
            "all_posterior_evaluable_n": len(tested),
            "common_with_lead_beta_n": len(common),
            "magma_top2_posterior_evaluable_n": len(parent),
            "at_exact_0p8": int(tested.posterior.eq(.8).sum()),
            "at_exact_0p9": int(tested.posterior.eq(.9).sum()),
        }
        for threshold in (.8, .9):
            common["positive"] = common.posterior.gt(threshold)
            m = common.event_id.isin(model_ids)
            s = common.event_id.isin(lead_ids)
            a = {"endpoint": endpoint, "posterior_threshold": threshold,
                 "panel": "A_B_global_model_vs_lead", "universe_n": len(common),
                 "model_selected_n": int(m.sum()),
                 "model_positive_n": int((m & common.positive).sum()),
                 "sqtl_selected_n": int(s.sum()),
                 "sqtl_positive_n": int((s & common.positive).sum()),
                 "selection_overlap_n": int((m & s).sum())}
            rows.append(a)

            parent["positive"] = parent.posterior.gt(threshold)
            parent["model_selected"] = parent.event_id.isin(model_ids)
            parent["selected_positive"] = parent.model_selected & parent.positive
            x = parent.model_selected
            y = parent.positive
            c = contingency(int((x & y).sum()), int((x & ~y).sum()),
                            int((~x & y).sum()), int((~x & ~y).sum()))
            region = grouped_stats(parent, "block")
            gene = grouped_stats(parent, "gene_id")
            rows.append({"endpoint": endpoint, "posterior_threshold": threshold,
                         "panel": "C_D_within_MAGMA_global_top2", "universe_n": len(parent),
                         "model_selected_n": c["selected_n"],
                         "model_positive_n": c["selected_hits"],
                         "other_magma_n": c["other_n"],
                         "other_magma_positive_n": c["other_hits"],
                         "event_or": c["fisher_or"],
                         "event_or_ci_low": c["or_ci_low"],
                         "event_or_ci_high": c["or_ci_high"],
                         "event_p_one_sided": c["fisher_p_greater"],
                         "region_selected_n": region.get("selected_n"),
                         "region_selected_hits": region.get("selected_hits"),
                         "region_other_n": region.get("other_n"),
                         "region_other_hits": region.get("other_hits"),
                         "region_or": region.get("fisher_or"),
                         "region_p_one_sided": region.get("fisher_p_greater"),
                         "gene_selected_n": gene.get("selected_n"),
                         "gene_selected_hits": gene.get("selected_hits"),
                         "gene_other_n": gene.get("other_n"),
                         "gene_other_hits": gene.get("other_hits"),
                         "gene_or": gene.get("fisher_or"),
                         "gene_p_one_sided": gene.get("fisher_p_greater")})
            if threshold == .9:
                common["model_selected"] = m
                common["sqtl_selected"] = s
                common["magma_selected"] = common.event_id.isin(magma_ids)
                common["endpoint"] = endpoint
                common["stratum"] = "common_with_lead"
                parent["sqtl_selected"] = parent.event_id.isin(lead_ids)
                parent["magma_selected"] = True
                parent["endpoint"] = endpoint
                parent["stratum"] = "magma_top2"
                detail.extend((common, parent))

    result = pd.DataFrame(rows)
    result.to_csv(HERE / "global_top2_figure_stats.csv", index=False)
    pd.concat(detail, ignore_index=True).to_csv(HERE / "global_top2_figure_events.csv", index=False)
    (HERE / "global_top2_figure_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(result.to_string(index=False))
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
