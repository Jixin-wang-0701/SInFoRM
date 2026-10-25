"""Build a one-row-per-junction, plain-language model input/output index.

Counts summarize each annotation block; complete feature values remain in
matrices/*.npz with matching feature-name CSV files.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
MATRICES = ROOT / "matrices"
SCORES = ROOT / "scores"
OUTPUT = ROOT / "input_output" / "junction_input_output_table.csv"


def main() -> None:
    rows = pd.read_csv(MATRICES / "rows.csv")
    full = pd.read_csv(SCORES / "full_fit_scores.csv")
    loco = pd.read_csv(SCORES / "loco_scores.csv")
    assert len(rows) == len(full) == len(loco) == 95_750
    assert rows.bed_id.is_unique
    assert rows.bed_id.equals(full.bed_id) and rows.bed_id.equals(loco.bed_id)
    assert np.allclose(rows.ZSTAT, full.ZSTAT) and np.allclose(rows.ZSTAT, loco.ZSTAT)

    x_expression = sparse.load_npz(MATRICES / "X_expression.npz")
    x_go = sparse.load_npz(MATRICES / "X_go.npz")
    x_pfam = sparse.load_npz(MATRICES / "X_pfam_rebuilt.npz")
    assert x_expression.shape == (len(rows), 17)
    assert x_go.shape == (len(rows), 15_716)
    assert x_pfam.shape == (len(rows), 5_450)

    table = pd.DataFrame({
        "junction_id": rows.bed_id,
        "chromosome_identifier_only": rows.Chromosome,
        "junction_start_identifier_only": rows.j_start,
        "junction_end_identifier_only": rows.j_end,
        "MAGMA_ZSTAT_training_response": rows.ZSTAT,
        "DNA_sequence_or_motif_predictor_in_this_model": False,
        "single_cell_expression_annotation_available": rows.has_expression,
        "single_cell_expression_nonzero_features": x_expression.getnnz(axis=1),
        "GO_pathway_annotation_available": rows.has_go,
        "GO_pathway_nonzero_features": x_go.getnnz(axis=1),
        "Pfam_domain_nonzero_features": x_pfam.getnnz(axis=1),
        "SInFoRM_full_fit_feature_only_score": full.feature_only_abs_score,
        "SInFoRM_full_fit_rank": full.rank_desc,
        "SInFoRM_LOCO_feature_only_score": loco.feature_only_abs_score,
        "SInFoRM_LOCO_rank": loco.rank_desc,
    })
    assert table.junction_id.is_unique and not table.isna().any().any()
    assert (table.single_cell_expression_nonzero_features <= 17).all()
    assert (table.GO_pathway_nonzero_features <= 15_716).all()
    assert (table.Pfam_domain_nonzero_features <= 5_450).all()
    OUTPUT.parent.mkdir(exist_ok=True)
    table.to_csv(OUTPUT, index=False)
    print(OUTPUT, "rows:", len(table), "columns:", len(table.columns))


if __name__ == "__main__":
    main()
