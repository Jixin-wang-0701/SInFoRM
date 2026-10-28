# model and analysis

This folder contains the data, model coefficients, per-junction scores, figure, and analysis scripts for the six-panel figure. `manuscript_text.md` contains a short English legend, Results paragraph, and reviewer response.

## Folder contents

| Path | Contents |
| --- | --- |
| `matrices/` | Aligned `rows.csv`, original expression and GO sparse feature matrices (`.npz`), rebuilt Pfam sparse matrix (`.npz` and Matrix Market `.mtx`), feature names, and rebuilt Pfam `.rds`. Matrix rows follow `rows.csv`; columns follow the respective feature-name CSV. Shapes are 95,750 × 17, 95,750 × 15,716, and 95,750 × 5,450. |
| `coefficients/` | Full-fit coefficients in readable CSV and original NumPy arrays; consolidated expression, GO and Pfam selected coefficients for each LOCO fold. The full-fit coefficients were learned after variance scaling and then applied to raw matrix columns to reproduce the original scoring convention. |
| `scores/` | Full-fit and LOCO scores in CSV and Excel. |
| `input_output/` | Manuscript-ready plain-language input/output table, a 95,750-row junction-level index, and the full-data fitted SInFoRM **top 2%** junction list in CSV and Excel. |
| `figure_data/` | Figure statistics, hit/event tables, rank curves and threshold sweep, plus an exploratory comparison of four alternative sQTL ordering rules for panel A. |
| `code/` | Original fitting, LOCO, colocalization-comparison and figure scripts; package and Excel builders; verification script. |



