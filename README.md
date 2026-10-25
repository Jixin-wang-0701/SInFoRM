# Rebuilt-Pfam six-panel figure: model and analysis package

This folder contains the data, model coefficients, per-junction scores, figure, and analysis scripts for the six-panel figure. `manuscript_text.md` contains a short English legend, Results paragraph, and reviewer response.

## Two score sets in the figure

| Figure panels | Score source | Interpretation |
| --- | --- | --- |
| A–D | `scores/full_fit_scores.csv` / Excel **Full fit** | Expression, GO and rebuilt-Pfam block fits on all 95,750 junctions; these are **in-sample** descriptive ranks for the colocalization comparisons. |
| E–F | `scores/loco_scores.csv` / Excel **LOCO** | Each junction was predicted with its chromosome held out. GO/Pfam selection and tuning used the other 21 chromosomes; expression predictions came from the original nested LOCO run. The ≥3/4 training-chromosome GO/Pfam coverage screen was selected **post hoc**. |

The Excel workbook is `scores/junction_scores_full_fit_and_loco.xlsx`. Each sheet has one row per junction and includes three signed block contributions, each block's absolute contribution, their signed sum, the absolute sum used for ranking, and descending rank. **The combined score is `abs(expression + GO + Pfam)`, not the sum of the three absolute values.** No genetic score is added at prediction time. MAGMA ZSTAT is retained as an outcome/reference column and was used for fitting; the score is therefore not genetically independent.

## Folder contents

| Path | Contents |
| --- | --- |
| `matrices/` | Aligned `rows.csv`, original expression and GO sparse feature matrices (`.npz`), rebuilt Pfam sparse matrix (`.npz` and Matrix Market `.mtx`), feature names, and rebuilt Pfam `.rds`. Matrix rows follow `rows.csv`; columns follow the respective feature-name CSV. Shapes are 95,750 × 17, 95,750 × 15,716, and 95,750 × 5,450. |
| `coefficients/` | Full-fit coefficients in readable CSV and original NumPy arrays; consolidated expression, GO and Pfam selected coefficients for each LOCO fold. The full-fit coefficients were learned after variance scaling and then applied to raw matrix columns to reproduce the original scoring convention. |
| `scores/` | Full-fit and LOCO scores in CSV and Excel. |
| `input_output/` | Manuscript-ready plain-language input/output table, a 95,750-row junction-level index, and the full-data fitted SInFoRM **top 2%** junction list in CSV and Excel. |
| `loco/` | Out-of-fold predictions, fold-level selected features and tuning records, including provenance of the reused expression component. |
| `figure_data/` | Figure statistics, hit/event tables, rank curves and threshold sweep, plus an exploratory comparison of four alternative sQTL ordering rules for panel A. |
| `figure/` | Six-panel PNG, editable PowerPoint and detailed figure caption. |
| `code/` | Original fitting, LOCO, colocalization-comparison and figure scripts; package and Excel builders; verification script. |
| `score_checks.json`, `file_manifest.csv` | Score reconstruction checks and SHA-256 file inventory. |

## Reproduce and check

From this folder, run `python code/verify_release.py` to check the score CSVs against Excel and the stored source predictions. Run `python code/build_score_workbook.py` to regenerate the workbook, `python code/build_input_output_table.py` to regenerate the junction-level input/output index, or the bundled Node runtime with `code/build_top2_junction_table.mjs` to regenerate the top-2% table. The original fitting and figure scripts in `code/` document the full workflow but retain their original project-relative or external data paths; external colocalization and sQTL source files are represented here by the derived event/figure tables.

The top-2% table selects exactly 1,915 junctions by the **full-fit** absolute SInFoRM score (last included score 0.10428332064126088). LOCO scores/ranks in that table are reference columns and do not determine inclusion. This full-fit ranking is in-sample.

The colocalization panels and the LOCO MAGMA panels have different eligible-junction denominators. Event-level colocalization tests should not be read as independent region-level replication. At MAGMA P<5×10⁻⁸, the combined LOCO AUPRC is 0.003985 versus the 0.002642 prevalence baseline; GO alone has AUPRC 0.004441. These exploratory results do not show that combining the blocks uniformly improves ranking.
