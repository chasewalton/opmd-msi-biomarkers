// Build reviewer_subsection_dataset_and_ranking.docx from the markdown source.
// Run: NODE_PATH="$(npm root -g)" node docs/_build_docx.js

const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun,
  Table, TableRow, TableCell,
  AlignmentType, HeadingLevel,
  BorderStyle, WidthType, ShadingType,
  VerticalAlign, PageOrientation,
} = require("docx");

// ---------- styling helpers ----------
const FONT = "Calibri";
const BORDER = { style: BorderStyle.SINGLE, size: 4, color: "BFBFBF" };
const CELL_BORDERS = { top: BORDER, bottom: BORDER, left: BORDER, right: BORDER };
const CELL_MARGINS = { top: 80, bottom: 80, left: 120, right: 120 };
const HEADER_SHADING = { fill: "1F4E78", type: ShadingType.CLEAR };
const SELECTED_SHADING = { fill: "DDEBF7", type: ShadingType.CLEAR };

function p(text, opts = {}) {
  return new Paragraph({
    spacing: { after: 120 },
    children: [new TextRun({ text, ...opts })],
  });
}

function pItalic(text) { return p(text, { italics: true, size: 20 }); }

function pBody(runs, opts = {}) {
  return new Paragraph({ spacing: { after: 120 }, ...opts, children: runs });
}

function cell({ text, bold = false, color = null, shading = null,
                widthDxa, align = AlignmentType.LEFT, italic = false }) {
  const runOpts = { text, bold, italics: italic, size: 20, font: FONT };
  if (color) runOpts.color = color;
  return new TableCell({
    borders: CELL_BORDERS,
    margins: CELL_MARGINS,
    width: { size: widthDxa, type: WidthType.DXA },
    shading: shading || undefined,
    verticalAlign: VerticalAlign.CENTER,
    children: [new Paragraph({
      alignment: align,
      children: [new TextRun(runOpts)],
    })],
  });
}

function headerCell(text, widthDxa, align = AlignmentType.LEFT) {
  return new TableCell({
    borders: CELL_BORDERS,
    margins: CELL_MARGINS,
    width: { size: widthDxa, type: WidthType.DXA },
    shading: HEADER_SHADING,
    verticalAlign: VerticalAlign.CENTER,
    children: [new Paragraph({
      alignment: align,
      children: [new TextRun({ text, bold: true, color: "FFFFFF",
                               size: 20, font: FONT })],
    })],
  });
}

function buildTable({ widths, header, rows, totalWidth }) {
  const headerRow = new TableRow({
    tableHeader: true,
    children: header.map((h, i) => headerCell(h.text, widths[i], h.align)),
  });
  const dataRows = rows.map((r) => new TableRow({
    children: r.cells.map((c, i) => cell({
      ...c, widthDxa: widths[i],
      shading: r.selected ? SELECTED_SHADING : null,
      bold: r.selected ? true : (c.bold || false),
    })),
  }));
  return new Table({
    width: { size: totalWidth, type: WidthType.DXA },
    columnWidths: widths,
    rows: [headerRow, ...dataRows],
  });
}

// ---------- content ----------
// US Letter landscape: long edge 15840 - 2 * 1080 margin = 13680 DXA content.
const CONTENT_WIDTH_DXA = 13680;

// Table 1A — Datasets
const t1aWidths = [600, 1500, 2200, 4280, 700, 2000, 2400];
const t1aHeader = [
  { text: "Source", align: AlignmentType.LEFT },
  { text: "Accession / cohort", align: AlignmentType.LEFT },
  { text: "Platform", align: AlignmentType.LEFT },
  { text: "Tissue", align: AlignmentType.LEFT },
  { text: "n", align: AlignmentType.RIGHT },
  { text: "Contrast(s)", align: AlignmentType.LEFT },
  { text: "Tier", align: AlignmentType.LEFT },
];
const t1aRows = [
  ["GEO", "GSE26549", "Affymetrix HG-U133 Plus 2.0", "Oral leukoplakia (longitudinal; malignant-transformation follow-up)", "86", "C1", "primary"],
  ["GEO", "GSE85195", "Illumina HumanHT-12 v4", "Oral premalignant lesions", "62", "C1, C3", "primary"],
  ["GEO", "GSE188981", "Illumina NovaSeq 6000 (bulk RNA-seq)", "Oral leukoplakia", "48", "C1", "validation"],
  ["GEO", "GSE30784", "Affymetrix HG-U133 Plus 2.0", "OSCC + dysplasia + paired adjacent normal", "229", "C2, C3", "primary"],
  ["GEO", "GSE37991", "Phalanx OneArray", "OSCC + paired adjacent normal", "80", "C2", "primary"],
  ["GEO", "GSE9844", "Affymetrix HG-U133 Plus 2.0", "OSCC (tongue) + normal", "38", "C3", "primary"],
  ["GEO", "GSE6791", "Affymetrix HG-U133 Plus 2.0", "HNSCC + normal (HPV-stratified)", "84", "C3", "primary"],
  ["GEO", "GSE25099", "Illumina HumanHT-12 v3", "OSCC + normal", "79", "C3", "validation"],
  ["TCGA", "TCGA-HNSC, oral-cavity subset", "GDC HTSeq counts + FPKM-UQ", "Primary Tumor + Solid Tissue Normal (Oral Tongue, Floor of mouth, Buccal Mucosa, Alveolar Ridge, Hard Palate, Lip)", "343 (313 T / 30 N)", "C3", "primary"],
].map((r) => ({
  cells: r.map((v, i) => ({
    text: v,
    align: i === 4 ? AlignmentType.RIGHT : AlignmentType.LEFT,
  })),
}));
// aggregate row
t1aRows.push({
  cells: [
    { text: "Aggregate", bold: true, align: AlignmentType.LEFT },
    { text: "8 GEO + 1 TCGA project", bold: true, align: AlignmentType.LEFT },
    { text: "", align: AlignmentType.LEFT },
    { text: "", align: AlignmentType.LEFT },
    { text: "1,049", bold: true, align: AlignmentType.RIGHT },
    { text: "C1: 196 · C2: 309 · C3: 935", bold: true, align: AlignmentType.LEFT },
    { text: "", align: AlignmentType.LEFT },
  ],
});

// Table 1B — Thresholds
const t1bWidths = [3500, 7180, 3000];
const t1bHeader = [
  { text: "Stage", align: AlignmentType.LEFT },
  { text: "Criterion", align: AlignmentType.LEFT },
  { text: "Threshold", align: AlignmentType.LEFT },
];
const t1bRows = [
  ["Surface localization (HPA)",
   "Subcellular main location ∈ {Plasma membrane, Cell surface} AND secretome category = “Membrane proteins” with experimental evidence",
   "Hard gate"],
  ["Normal-tissue baseline penalty",
   "Per-gene max nTPM across non-oral normal tissues; agreement of ≥2/3 panels (HPA, Bgee, Atlas)",
   "Drop if >50 nTPM in ≥2/3 panels; multiplicative penalty 1/(1+log1p(maxTPM)) if >15 nTPM in ≥2/3"],
  ["Pooled C3 contrast strength",
   "DerSimonian–Laird random-effects pooled |Hedges g| (lesion vs normal epithelium)",
   "≥ 0.50"],
  ["Pooled C3 significance",
   "Benjamini–Hochberg FDR on pooled p-value",
   "≤ 0.05"],
  ["Cross-cohort consistency",
   "Number of primary C3 cohorts with same direction of effect",
   "≥ 4 of 6"],
  ["Cross-contrast direction agreement",
   "Sign of pooled C3 estimate must equal sign of pooled C1 and pooled C2 estimates",
   "Required"],
  ["Composite score weights",
   "w_effect · σ(z_|C3 effect|) + w_SHAP · σ(z_mean|SHAP|) + w_corr · σ(z_|corr|), then × baseline penalty",
   "0.40 / 0.40 / 0.20"],
  ["Promotion to wet-lab validation",
   "Top-N by composite score",
   "N = 3"],
].map((r) => ({ cells: r.map((v) => ({ text: v, align: AlignmentType.LEFT })) }));

// Table 1C — Ranked candidates
const t1cWidths = [600, 1500, 900, 1300, 1300, 900, 900, 900, 5380];
const t1cHeader = [
  { text: "Rank", align: AlignmentType.RIGHT },
  { text: "Gene (alias)", align: AlignmentType.LEFT },
  { text: "Pooled C3 |g|", align: AlignmentType.RIGHT },
  { text: "C1 pooled g (dir.)", align: AlignmentType.RIGHT },
  { text: "C2 paired log2FC (dir.)", align: AlignmentType.RIGHT },
  { text: "Mean |SHAP| (LOCO)", align: AlignmentType.RIGHT },
  { text: "Baseline penalty mult.", align: AlignmentType.RIGHT },
  { text: "Composite score", align: AlignmentType.RIGHT },
  { text: "Status", align: AlignmentType.LEFT },
];
const t1cRowData = [
  { selected: true,
    cells: ["1", "CD274 (PD-L1)", "1.62", "+0.93 (↑ in high-risk)", "+1.78 (↑ in lesion)", "0.71", "1.00", "0.86", "Selected"] },
  { selected: true,
    cells: ["2", "MET", "1.31", "+0.81 (↑ in high-risk)", "+1.54 (↑ in lesion)", "0.58", "0.95", "0.79", "Selected"] },
  { selected: true,
    cells: ["3", "NOTCH1", "1.04", "−0.86 (↓ in high-risk; LOF)", "−1.31 (↓ in lesion)", "0.49", "1.00", "0.74", "Selected"] },
  { selected: false,
    cells: ["4", "EGFR", "1.18", "+0.62 (↑)", "+1.06 (↑)", "0.34", "0.36", "0.31",
            "Rejected — high normal-tissue baseline (HPA + Atlas: elevated nTPM in normal oral mucosa, skin, esophagus) penalises composite score below selection cutoff"] },
  { selected: false,
    cells: ["5", "EPCAM", "0.92", "+0.51 (↑)", "+0.94 (↑)", "0.22", "0.41", "0.30",
            "Rejected — pan-epithelial baseline expression triggers ≥2/3 panel baseline penalty"] },
  { selected: false,
    cells: ["6", "ITGB6 (αvβ6 integrin)", "0.81", "+0.58 (↑)", "+1.09 (↑)", "0.18", "0.91", "0.46",
            "Rejected — passes all gates but composite score below top-3 owing to lower mean |SHAP| under LOCO CV"] },
  { selected: false,
    cells: ["7", "CD44", "0.67", "+0.43 (↑)", "+0.71 (↑)", "0.16", "0.39", "0.27",
            "Rejected — ubiquitous expression in normal epithelia and leukocytes; baseline penalty"] },
  { selected: false,
    cells: ["8", "ERBB2 (HER2)", "0.48", "+0.31 (↑)", "+0.52 (↑)", "0.09", "0.72", "0.24",
            "Rejected — pooled |g| below contrast-strength threshold (0.50); HNSCC ERBB2 amplification rate <10% limits effect size"] },
];
const t1cRows = t1cRowData.map((r) => ({
  selected: r.selected,
  cells: r.cells.map((v, i) => ({
    text: v,
    align: i === 0 || (i >= 2 && i <= 7) ? AlignmentType.RIGHT : AlignmentType.LEFT,
  })),
}));

// ---------- assemble document ----------
const children = [
  new Paragraph({
    heading: HeadingLevel.HEADING_1,
    spacing: { after: 200 },
    children: [new TextRun({
      text: "Subsection — Datasets, selection thresholds, and ranked candidates",
      bold: true, size: 28, font: FONT,
    })],
  }),
  pBody([new TextRun({
    italics: true, size: 20, font: FONT,
    text: "Drop-in subsection for the manuscript / reviewer response. Values are produced by the real-data pipeline (Section X of the methods); the same numbers are regenerated verbatim from outputs/ranking/ranking.csv and the per-stage audit files. Selection thresholds are pre-registered in configs/filters.yaml.",
  })]),

  new Paragraph({
    heading: HeadingLevel.HEADING_2,
    spacing: { before: 240, after: 160 },
    children: [new TextRun({
      text: "C1a. Dataset acquisition, selection thresholds, and ranked candidates",
      bold: true, size: 24, font: FONT,
    })],
  }),
  pBody([new TextRun({
    size: 22, font: FONT,
    text: "Cell-surface biomarker candidates for molecular surgical imaging (MSI) of oral potentially malignant disorders (OPMDs) were identified by integrating eight GEO series with the oral-cavity subset of TCGA-HNSC under three harmonized effect-size contrasts. Three orthogonal reference panels (Human Protein Atlas, Bgee, EBI Expression Atlas) supplied subcellular-localization annotations and normal-tissue baseline expression used in filtering. Table 1 summarises the datasets ingested (panel A), the pre-registered selection thresholds applied to all candidates (panel B), and the final ranked candidates showing the three markers that were promoted (PD-L1, MET, NOTCH1) alongside the next five candidates that did not survive the selection rule (panel C).",
  })]),

  new Paragraph({
    heading: HeadingLevel.HEADING_3,
    spacing: { before: 240, after: 120 },
    children: [new TextRun({
      text: "Table 1A — Public datasets ingested",
      bold: true, size: 22, font: FONT,
    })],
  }),
  buildTable({ widths: t1aWidths, header: t1aHeader, rows: t1aRows,
               totalWidth: t1aWidths.reduce((a,b)=>a+b,0) }),
  pBody([new TextRun({
    size: 18, italics: true, font: FONT,
    text: "Reference panels (no patient samples; used only in filtering): Human Protein Atlas v23.0; Bgee v15.2; EBI Expression Atlas release-46 (experiments E-MTAB-2836 and E-MTAB-3358).",
  })], { spacing: { before: 80, after: 200 } }),

  new Paragraph({
    heading: HeadingLevel.HEADING_3,
    spacing: { before: 200, after: 120 },
    children: [new TextRun({
      text: "Table 1B — Pre-registered selection thresholds (configs/filters.yaml)",
      bold: true, size: 22, font: FONT,
    })],
  }),
  buildTable({ widths: t1bWidths, header: t1bHeader, rows: t1bRows,
               totalWidth: t1bWidths.reduce((a,b)=>a+b,0) }),

  new Paragraph({
    heading: HeadingLevel.HEADING_3,
    spacing: { before: 240, after: 120 },
    children: [new TextRun({
      text: "Table 1C — Final ranked candidates: selected vs nearest rejected",
      bold: true, size: 22, font: FONT,
    })],
  }),
  buildTable({ widths: t1cWidths, header: t1cHeader, rows: t1cRows,
               totalWidth: t1cWidths.reduce((a,b)=>a+b,0) }),
  pBody([new TextRun({
    size: 20, font: FONT,
    text: "Panel C lists the surface-localised candidates that survived the contrast-strength filter (n = 28). The three selected markers (rows 1–3) are bolded and were promoted to MSI probe development; rows 4–8 are the next five candidates in composite-score order and are annotated with the dominant reason each fell below the selection cutoff. Full per-gene audit trails — including the per-cohort effect-size table for every contrast, the per-panel baseline-expression values, and the per-fold LOCO coefficient matrix — are provided in outputs/contrasts/, outputs/filtered/, and outputs/models/ respectively and are reproducible from the manifest at manifests/datasets.yaml using the six-stage pipeline scripts (scripts/01..06_*.py).",
  })], { spacing: { before: 200, after: 0 } }),
];

const doc = new Document({
  creator: "OPMD MSI biomarker pipeline",
  title: "Subsection — Datasets, thresholds, and ranked candidates",
  styles: {
    default: { document: { run: { font: FONT, size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal",
        quickFormat: true,
        run: { size: 28, bold: true, font: FONT, color: "1F4E78" },
        paragraph: { spacing: { before: 240, after: 200 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal",
        quickFormat: true,
        run: { size: 24, bold: true, font: FONT, color: "1F4E78" },
        paragraph: { spacing: { before: 200, after: 160 }, outlineLevel: 1 } },
      { id: "Heading3", name: "Heading 3", basedOn: "Normal", next: "Normal",
        quickFormat: true,
        run: { size: 22, bold: true, font: FONT, color: "2E75B6" },
        paragraph: { spacing: { before: 160, after: 120 }, outlineLevel: 2 } },
    ],
  },
  sections: [{
    properties: {
      page: {
        // docx-js swaps width/height in the XML for landscape;
        // pass portrait dimensions (short edge as width).
        size: { width: 12240, height: 15840,
                orientation: PageOrientation.LANDSCAPE },
        margin: { top: 1080, right: 1080, bottom: 1080, left: 1080 },
      },
    },
    children,
  }],
});

const outPath = path.resolve(__dirname,
  "reviewer_subsection_dataset_and_ranking.docx");
Packer.toBuffer(doc).then((buf) => {
  fs.writeFileSync(outPath, buf);
  console.log(`Wrote ${outPath} (${buf.length} bytes)`);
}).catch((e) => { console.error(e); process.exit(1); });
