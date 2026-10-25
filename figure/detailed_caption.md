# Rebuilt Pfam six-panel figure: caption and analysis notes

Panels A–D retain the previously reviewed comparisons: A–B plot cumulative colocalization hits for the full-data expression + GO + rebuilt Pfam score against lead sQTL absolute effect size in each external evaluable set; C–D evaluate colocalization enrichment among junctions already in the global MAGMA top 2%. Their analysis sets and denominators differ. A–B use region-level tests for cumulative-curve area. D displays both event-level and region-level tests, which lead to different conclusions.

Panels E–F use 22-fold leave-one-chromosome-out (LOCO) predictions for all 95,750 junctions. In each outer fold, GO/Pfam candidates must appear on at least 16 of the 21 training chromosomes and at least 30 training junctions; feature ranking and LASSO tuning occur within the training chromosomes. This three-quarter coverage rule was chosen after looking at earlier results, so these panels are **exploratory sensitivity analyses**. They do not establish an independently validated improvement. The displayed scores contain no separately added genetic or MAGMA term at prediction time; the LASSO coefficients were learned from MAGMA ZSTAT in training folds.

For E–F, positives are defined by one-sided `P = norm.sf(ZSTAT)` below 5×10⁻⁴, 5×10⁻⁵, 5×10⁻⁶, 5×10⁻⁷, or 5×10⁻⁸. The counts are 1,777, 970, 595, 384, and 253 positive junctions, respectively. The existing out-of-fold scores are held fixed while labels change. E shows pooled out-of-fold average precision divided by the positive prevalence at each threshold. The values for the combined score are 1.24×, 1.26×, 1.44×, 1.57×, and 1.51× random ranking. These are descriptive ratios; confidence intervals and significance tests for the threshold trend have not been calculated. Neighboring junctions and genetic loci are correlated.

F shows the fraction of all positives recovered among the top 9,575 junctions (10% of 95,750) by each score. The random-ranking expectation is 10%. Expression, GO, and rebuilt Pfam are scored separately using the absolute value of their existing LOCO component contributions; the combination takes the absolute value of the signed sum. Ties at the top-10% boundary receive expected random allocation. F is a **component-only score comparison**, not a drop-one-feature-block ablation with model refitting. At 5×10⁻⁸, sensitivity is 11.5% for expression, 16.2% for GO, 9.9% for Pfam, and 17.0% for the combination. GO alone also has higher AUPRC than the combination at that threshold (1.68× versus 1.51× random prevalence), so the figure should not be interpreted as showing uniform gain from every block.

The omitted “fixed alpha” comparator refers to the all-feature LOCO run in which each held-out chromosome receives newly fitted scalers and LASSO coefficients, while the penalty parameter `alpha` stays at values selected before the outer folds, including an expression alpha selected with full-data cross-validation. Thus that comparator is not fully nested hyperparameter validation. It has been removed from this six-panel figure as requested; the earlier full-feature results remain in the analysis files.

## Why these plot forms

- [Saito and Rehmsmeier, *PLOS ONE* (2015)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0118432) motivate precision–recall based evaluation for rare positive classes.
- [Brabec et al., *ICCS* (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7303692/) show why the class imbalance ratio must accompany precision–recall comparisons. Accordingly, E reports AUPRC divided by each threshold's positive prevalence and labels the positive count.
- [Weeks et al., *Nature Genetics* (2023), PoPS](https://www.nature.com/articles/s41588-023-01443-6) use leave-one-chromosome-out scoring in a gene-prioritization setting.
- [Ding et al., *Human Genetics* (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12185218/) show feature-group AUPRC comparisons in a genomic prediction study, but their formal ablation removes groups and retrains. Panel F is explicitly labeled as a component-only comparison because its models are not retrained.

## Local data and reproduction

- `outputs/pfam_rebuild_npnt/loco_five_power_p_sweep.csv`
- `outputs/pfam_rebuild_npnt/loco_five_power_p_sweep.py`
- `outputs/pfam_rebuild_npnt/figure_build/make_rebuilt_pfam_six_panel.mjs`
- `outputs/pfam_rebuild_npnt/rank_enrichment_model_vs_sqtl_stats.csv`
- `outputs/pfam_rebuild_npnt/rank_enrichment_model_vs_sqtl_curves.csv`
- `outputs/pfam_rebuild_npnt/global_top2_figure_stats.csv`
- `outputs/pfam_rebuild_npnt/run_full_feature_loco.py`
