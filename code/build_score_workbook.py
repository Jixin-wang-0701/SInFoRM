"""Stream the two verified score CSV files into one Excel workbook.

The artifact-tool attempt exhausted an 8 GiB Node heap at XLSX export with
2.3 million populated cells; this constant-memory fallback preserves all rows.
"""

from __future__ import annotations

import csv
from pathlib import Path

import xlsxwriter


ROOT = Path(__file__).resolve().parents[1]
SCORES = ROOT / "scores"
OUTPUT = SCORES / "junction_scores_full_fit_and_loco.xlsx"
EXPECTED = [
    "bed_id", "Chromosome", "ZSTAT", "expression_component", "go_component",
    "pfam_component", "expression_abs_score", "go_abs_score", "pfam_abs_score",
    "signed_sum", "feature_only_abs_score", "rank_desc",
]


def main() -> None:
    wb = xlsxwriter.Workbook(OUTPUT, {"constant_memory": True})
    header = wb.add_format({"bold": True, "font_color": "white", "bg_color": "#1F4E78"})
    numeric = wb.add_format({"num_format": "0.0000000000"})
    integer = wb.add_format({"num_format": "0"})
    for tab, filename in (("Full fit", "full_fit_scores.csv"), ("LOCO", "loco_scores.csv")):
        sheet = wb.add_worksheet(tab)
        sheet.freeze_panes(1, 0)
        sheet.set_column(0, 0, 31)
        sheet.set_column(1, 1, 15)
        sheet.set_column(2, 10, 21, numeric)
        sheet.set_column(11, 11, 12, integer)
        sheet.set_row(0, 29)
        sheet.write_row(0, 0, EXPECTED, header)
        with (SCORES / filename).open(newline="", encoding="utf-8") as handle:
            reader = csv.reader(handle)
            assert next(reader) == EXPECTED, filename
            count = 0
            for count, values in enumerate(reader, start=1):
                assert len(values) == len(EXPECTED), (filename, count)
                sheet.write_string(count, 0, values[0])
                sheet.write_string(count, 1, values[1])
                for col in range(2, 11):
                    sheet.write_number(count, col, float(values[col]), numeric)
                sheet.write_number(count, 11, int(values[11]), integer)
            assert count == 95_750, (filename, count)
        sheet.autofilter(0, 0, count, len(EXPECTED) - 1)
        print(tab, count)

    notes = wb.add_worksheet("ReadMe")
    notes.set_column(0, 0, 26)
    notes.set_column(1, 1, 115)
    notes.write_row(0, 0, ["Field", "Meaning"], header)
    descriptions = [
        ("Full fit", "All 95,750 junctions used to fit separate expression, GO and rebuilt Pfam models; in-sample scores for panels A-D."),
        ("LOCO", "Each row predicted with its chromosome held out; GO/Pfam screened in 21 training chromosomes; panels E-F."),
        ("component", "Signed feature-block contribution."),
        ("*_abs_score", "Absolute value of an individual block contribution."),
        ("signed_sum", "Expression + GO + Pfam signed contributions."),
        ("feature_only_abs_score", "Absolute value of signed_sum; no genetic component at prediction time."),
        ("rank_desc", "Descending feature-only score rank; ties retain source row order."),
        ("Sources", "See matrices/, coefficients/, loco/, score_checks.json and README.md."),
    ]
    for row, record in enumerate(descriptions, start=1):
        notes.write_row(row, 0, record)
    wb.close()
    print(OUTPUT)


if __name__ == "__main__":
    main()
