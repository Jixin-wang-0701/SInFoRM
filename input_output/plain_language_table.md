# Plain-language input/output table for the rebuilt-Pfam SInFoRM model

**Suggested manuscript table.** Each data row represents one splice junction; the row below defines its fields.

| Analysis unit | Training response | Sequence / junction information | Single-cell predictor | Pathway predictor | Domain predictor | Model output |
| --- | --- | --- | --- | --- | --- | --- |
| One junction (`chr:start:end`) | COPD event-level MAGMA `ZSTAT` | Coordinates identify the junction; **no separate DNA sequence or splice-motif predictor** was fitted in this three-block model | Lung cell-type expression values (17 candidate features) | GO term membership (15,716 candidate features) | Rebuilt Pfam protein-domain membership (5,450 candidate features) | SInFoRM score = `abs(expression contribution + GO contribution + Pfam contribution)`; junctions ranked from highest to lowest |

**Table note.** Expression features are quantitative; GO and Pfam features are binary. MAGMA `ZSTAT` is the fitting response, not a prediction-time input. Junction coordinates, raw sQTL effect sizes, colocalization probabilities, and SMR results are not predictors in this model. The full-data score is in-sample; the LOCO score is predicted with each junction's chromosome held out. The GO/Pfam LOCO feature screen uses the training chromosomes only.

**Reviewer response.** We added a one-row-per-junction input/output table that identifies event-level MAGMA `ZSTAT` as the training response, cell-type expression, GO and rebuilt Pfam annotations as the functional predictors, and the ranked SInFoRM score as the output. We clarify that the current model has no separate DNA sequence or splice-motif feature block; genomic coordinates identify junctions but are not used to calculate the score. Raw sQTL, colocalization and SMR results are also excluded from the predictors.

The accompanying `junction_input_output_table.csv` contains all 95,750 junctions, their response, per-block nonzero annotation counts, and both full-fit and LOCO scores/ranks. The counts summarize the predictors; complete feature values remain in `../matrices/`.
