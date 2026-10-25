# Text for the six-panel figure

## Figure legend

**Figure X. Colocalization enrichment and chromosome-held-out performance of the expression + GO + rebuilt-Pfam score.** (A, B) Cumulative recall of colocalization hits (posterior >0.9), ranked by the full-data functional score or maximum GTEx sQTL |β| in their common eligible sets. The maximum-|β| comparison was selected post hoc; regional P values test curve area without adjustment for comparator selection. (C, D) Hit rates and odds ratios within the global MAGMA top 2%; event- and region-level P values are distinguished. (E) Pooled leave-one-chromosome-out (LOCO) AUPRC relative to MAGMA-positive prevalence across one-sided P thresholds. (F) Top-10% capture sensitivity for the three feature blocks and their combined score. E–F use held-out feature-only predictions. The ≥3/4 training-chromosome GO/Pfam coverage rule was chosen post hoc; panel denominators differ.

## Results

In an exploratory comparison with maximum sQTL |β|, the full-data score yielded AUPRC 0.0638 versus 0.0331 for ICGC COPD × LTRC (regional P=0.031) and 0.0150 versus 0.0068 for FEV1/FVC × GTEx lung (P=0.021). These P values do not account for choosing this comparator after inspecting the lead-sQTL analysis. Within the MAGMA top 2%, the ICGC hit rate was 9/81 versus 4/157 (event P=0.0085; region P=0.341). Across five MAGMA thresholds, exploratory combined-score LOCO AUPRC was 1.24–1.57 times positive prevalence. At P<5×10⁻⁸, it was 0.00399 versus a 0.00264 baseline, with 17.0% of positives in the top 10%. GO alone reached AUPRC 0.00444, so combining blocks did not uniformly improve ranking.

## Response to reviewer

Thank you. We now report 22-chromosome LOCO AUPRC and top-decile sensitivity for an annotation-only score and its expression, GO, and rebuilt-Pfam components. Feature screening and tuning excluded each held-out chromosome. No MAGMA statistic or genomic baseline term was added to the held-out score. At MAGMA P<5×10⁻⁸ (253 positives), the combined score yielded AUPRC 0.00399 versus a prevalence baseline of 0.00264 and recovered 43/253 positives (17.0%) in its top decile. At the original P<0.025 definition, AUPRC was 0.09955 versus prevalence 0.09542. We clarify that MAGMA ZSTAT was the training response: these analyses measure held-out functional-annotation ranking of a GWAS-derived outcome, not information independent of GWAS. A separate, full-data conditional colocalization check within MAGMA top 2% was not significant at the genomic-region level (P=0.341 and 0.265 for the two endpoints). We therefore describe incremental evidence as exploratory rather than claiming a validated gain beyond GWAS. The GO/Pfam chromosome-coverage screen was also selected post hoc.

## Editorial checks before submission

- Replace “Figure X” with the final figure number.
- Retain the post hoc screening and differing-denominator statements if the figure is shortened.
