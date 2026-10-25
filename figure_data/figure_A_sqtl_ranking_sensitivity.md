# Figure A: sQTL ranking sensitivity (exploratory)

The original panel ranks the 525 ICGC COPD × LTRC lung common junctions against lead-sQTL `|β|` from the minimum-P GTEx pair. There are 13 PP>0.9 colocalization hits in five genomic blocks. All alternatives below use exactly the same junctions, hits, SInFoRM full-data score and one-sided regional sign-flip test of cumulative-curve area.

| sQTL ordering | sQTL AUPRC | SInFoRM − sQTL curve area | Regional P | Holm P across four tried orderings |
| --- | ---: | ---: | ---: | ---: |
| Lead-pair `|β|` (current figure) | 0.02455 | 0.23824 | 0.06250 | 0.18750 |
| Maximum `|β|` among sQTL pairs | 0.03314 | 0.22095 | **0.03125** | 0.12500 |
| Minimum nominal sQTL P | 0.02904 | 0.26095 | 0.09375 | 0.18750 |
| Maximum `|β/SE|` | 0.02753 | 0.26110 | 0.09375 | 0.18750 |

The SInFoRM AUPRC is 0.06377 and its cumulative-curve area is 0.72601 in every row. `−log10(P)` gives the same order as minimum P. Because only five positive genomic blocks contribute, the exact sign-flip P values have 1/32 resolution; 0.03125 is the smallest attainable one-sided P here.

The minimum-P ranking therefore does **not** make panel A significant. Maximum `|β|` reaches nominal P=0.03125, but this ordering was examined after seeing the original result. Its four-method Holm P=0.125, so it should remain a post-hoc sensitivity analysis rather than replace the primary comparator. In panel B, minimum-P sQTL ranking has curve area 0.70002 versus SInFoRM's 0.66558 (one-sided P=0.69866); switching both panels to P-based sQTL ordering would reverse the displayed advantage in B.

The test P values concern cumulative-curve area, not the difference in AUPRC. Full numeric results and curve points are in `figure_A_sqtl_ranking_sensitivity.csv` and `figure_A_sqtl_ranking_sensitivity_curves.csv`; `figure_A_sqtl_ranking_sensitivity.py` reproduces them.
