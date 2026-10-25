// Export the top 2% of the full-data fitted feature-only SInFoRM ranking.
// LOCO columns are reference values; they do not determine membership.
import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { createRequire } from 'node:module';

const runtimeModules = 'C:/Users/18553/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules';
const runtimeRequire = createRequire(path.join(runtimeModules, '_resolver.js'));
const { Workbook, SpreadsheetFile } = runtimeRequire('@oai/artifact-tool');
const release = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const sourceDir = path.join(release, 'input_output');
const scoreFile = path.join(release, 'scores', 'full_fit_scores.csv');
const indexFile = path.join(sourceDir, 'junction_input_output_table.csv');
const outputCsv = path.join(sourceDir, 'top2_sinform_junctions.csv');
const outputXlsx = path.join(sourceDir, 'top2_sinform_junctions.xlsx');

function readSimpleCsv(text) {
  const lines = text.trimEnd().split(/\r?\n/);
  const headers = lines[0].split(',');
  return lines.slice(1).map((line, n) => {
    const values = line.split(',');
    if (values.length !== headers.length) throw new Error(`Unexpected CSV width on row ${n + 2}`);
    return Object.fromEntries(headers.map((name, i) => [name, values[i]]));
  });
}

const scoreRows = readSimpleCsv(await fs.readFile(scoreFile, 'utf8'));
const indexRows = readSimpleCsv(await fs.readFile(indexFile, 'utf8'));
if (scoreRows.length !== 95750 || indexRows.length !== scoreRows.length) throw new Error('Unexpected source row count');
const indexById = new Map(indexRows.map((row) => [row.junction_id, row]));
if (indexById.size !== scoreRows.length) throw new Error('Duplicate junction ID');
const topN = Math.ceil(scoreRows.length * 0.02);
const selected = scoreRows.filter((row) => Number(row.rank_desc) <= topN);
selected.sort((a, b) => Number(a.rank_desc) - Number(b.rank_desc));
if (selected.length !== topN || selected.some((row, i) => Number(row.rank_desc) !== i + 1)) {
  throw new Error('Top-2% selection does not have exactly ranks 1..1915');
}

const headers = [
  'full_fit_rank', 'junction_id', 'chromosome', 'junction_start', 'junction_end',
  'MAGMA_ZSTAT_training_response', 'expression_component', 'GO_component',
  'Pfam_component', 'signed_sum', 'SInFoRM_abs_score', 'expression_abs_score',
  'GO_abs_score', 'Pfam_abs_score', 'LOCO_abs_score_reference',
  'LOCO_rank_reference', 'expression_nonzero_features', 'GO_nonzero_features',
  'Pfam_nonzero_features',
];
const values = selected.map((score) => {
  const input = indexById.get(score.bed_id);
  if (!input) throw new Error(`Missing input row: ${score.bed_id}`);
  const row = [
    Number(score.rank_desc), score.bed_id, score.Chromosome,
    Number(input.junction_start_identifier_only), Number(input.junction_end_identifier_only),
    Number(score.ZSTAT), Number(score.expression_component), Number(score.go_component),
    Number(score.pfam_component), Number(score.signed_sum),
    Number(score.feature_only_abs_score), Number(score.expression_abs_score),
    Number(score.go_abs_score), Number(score.pfam_abs_score),
    Number(input.SInFoRM_LOCO_feature_only_score), Number(input.SInFoRM_LOCO_rank),
    Number(input.single_cell_expression_nonzero_features),
    Number(input.GO_pathway_nonzero_features), Number(input.Pfam_domain_nonzero_features),
  ];
  if (Math.abs(Math.abs(row[9]) - row[10]) > 1e-12) throw new Error(`Score mismatch: ${score.bed_id}`);
  return row;
});
const cutoff = values.at(-1)[10];
if (values.some((row, i) => i > 0 && row[10] > values[i - 1][10] + 1e-15)) {
  throw new Error('Scores are not descending');
}

const csvEscape = (value) => {
  const item = String(value);
  return /[",\r\n]/.test(item) ? `"${item.replaceAll('"', '""')}"` : item;
};
await fs.writeFile(outputCsv,
  [headers, ...values].map((row) => row.map(csvEscape).join(',')).join('\n') + '\n', 'utf8');

const workbook = Workbook.create();
const sheet = workbook.worksheets.add('Top 2%');
sheet.getRange('A1:S1').values = [headers];
sheet.getRangeByIndexes(1, 0, values.length, headers.length).values = values;
sheet.freezePanes.freezeRows(1);
sheet.showGridLines = false;
sheet.getRange('A1:S1').format.fill = '#1F4E78';
sheet.getRange('A1:S1').format.font.color = '#FFFFFF';
sheet.getRange('A1:S1').format.font.bold = true;
sheet.getRange('A1:S1').format.wrapText = true;
sheet.getRange('A1:S1').format.rowHeight = 48;
sheet.getRange('A:A').format.columnWidth = 14;
sheet.getRange('B:B').format.columnWidth = 31;
sheet.getRange('C:C').format.columnWidth = 15;
sheet.getRange('D:E').format.columnWidth = 17;
sheet.getRange('F:S').format.columnWidth = 22;

const notes = workbook.worksheets.add('ReadMe');
notes.showGridLines = false;
notes.getRange('A1:B1').values = [['Item', 'Definition']];
notes.getRange('A1:B1').format.fill = '#1F4E78';
notes.getRange('A1:B1').format.font.color = '#FFFFFF';
notes.getRange('A1:B1').format.font.bold = true;
notes.getRange('A2:B8').values = [
  ['Selection', `Top ${topN} of ${scoreRows.length} junctions (2%) by full-data fitted SInFoRM_abs_score; ties use source row order.`],
  ['Score cutoff', cutoff],
  ['SInFoRM score', 'Absolute value of the signed sum of expression, GO and rebuilt Pfam contributions.'],
  ['Full fit', 'Model fitted on all 95,750 junctions; this top-2% ranking is descriptive and in-sample.'],
  ['LOCO reference', 'Held-out score/rank included for comparison only; it did not determine top-2% membership.'],
  ['MAGMA ZSTAT', 'Training response retained for interpretation; it is not added to the score at prediction time.'],
  ['Feature counts', 'Numbers of nonzero features per block. Complete values and names are in matrices/.'],
];
notes.getRange('A:A').format.columnWidth = 24;
notes.getRange('B:B').format.columnWidth = 115;

const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(outputXlsx);
process.stdout.write(`Top ${topN}/${scoreRows.length}; cutoff ${cutoff}; saved ${outputCsv} and ${outputXlsx}\n`);
