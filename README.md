# REG4 pan-cancer state-dominance analysis

**REG4 marks a lineage-restricted secretory epithelial state in gastrointestinal adenocarcinomas**

This repository contains the complete, reproducible analysis for a pan-cancer study of REG4 across 9,849 primary tumours from 34 TCGA PanCanAtlas cohorts. A single R script downloads all public data, runs every analysis, and produces all result tables and manuscript figures.

> **Status:** manuscript in preparation. Please do not share or cite without permission.

---

## Hypothesis

REG4's association with immune or stromal features varies across cancers, but its most coherent biology is epithelial/secretory state identity, especially in GI-lineage adenocarcinomas.

## Key findings

| Question | Finding |
|---|---|
| **Q1.** Is REG4 a GI epithelial-state marker? | REG4 tracks the goblet/secretory program (MUC2, SPINK4, FCGBP, AGR2; rho up to 0.90) rather than general epithelial content. Its coupling to intestinal lineage factors (CDX2, VIL1, HNF4A) is negative in colorectal cancer and positive in pancreatic cancer. |
| **Q2.** Do GI and non-GI tumours differ? | Epithelial/secretory dominance is enriched in GI adenocarcinomas: 5/6 vs 2/28 cohorts (Fisher OR 49, p = 4.3 × 10⁻⁴), robust to purity adjustment (p = 0.003). |
| **Q3.** Which axis dominates? | The epithelial/secretory axis dominates in READ, COAD, PAAD, ESCA_EAC and STAD (rho 0.62–0.76), robust to purity and score definition. Non-GI associations are weak and mostly purity-dependent. |
| **Q4.** Does REG4 mark mucinous states? | REG4 is 3–16× higher in mucinous GI tumours and tracks the secretory program even in non-mucinous GI tumours (rho 0.63–0.88). Mucinous breast carcinomas activate the secretory program **without** REG4. |
| **Q5.** Is REG4 prognostic beyond secretory state? | No, in GI cancers (pooled HR 0.96 per SD, p = 0.56). PRAD shows an exploratory protective signal (FDR = 0.064). |
| **Q6.** Does REG4 add information beyond the secretory program? | The colorectal immune/NK association is largely explained by MSI subtype (about 60% attenuation). PDAC shows an exploratory REG4-linked low-stroma signal. |

---

## Repository structure

```
reg4-pancancer/
├── README.md                    # this file
├── reg4_full_pipeline.R         # complete analysis: data download → results → figures
├── pipeline_log.txt             # printed output of the last full run
├── .gitignore                   # excludes downloaded data and R workspace
└── results/
    ├── *.csv                    # result tables (see "Outputs" below)
    ├── sessionInfo.txt          # R and package versions used
    └── figures/
        └── *.pdf                # manuscript figures
```

Downloaded data (`data/`) and the R workspace (`.RData`) are not tracked. They are regenerated automatically by the pipeline.

---

## Requirements

- **R ≥ 4.3** (developed with R 4.5.1)
- **R packages:** `UCSCXenaTools`, `data.table`, `survival`, `ppcor`, `ggplot2`, `patchwork`. Missing packages are installed automatically into the user library (`R_LIBS_USER`).
- **Internet access** on the first run, to download data from UCSC Xena and the GDC (about 5 MB in total). Later runs use the cached files in `data/`.
- **Memory and time:** under 4 GB RAM; a few minutes per run. Only the 50 required genes are downloaded, not the full 1.9 GB expression matrix.

---

## Quick start

### On IU Quartz (or any HPC with an R module)

```bash
module load r
git clone git@github.com:<username>/reg4-pancancer.git
cd reg4-pancancer
Rscript reg4_full_pipeline.R > pipeline_log.txt 2>&1
tail -30 pipeline_log.txt          # should end with "== DONE =="
```

On clusters where only login nodes have internet access, run the first (downloading) run on a login node.

### On a local machine

```bash
git clone git@github.com:<username>/reg4-pancancer.git
cd reg4-pancancer
Rscript reg4_full_pipeline.R > pipeline_log.txt 2>&1
```

Or, from an interactive R session in the repository folder:

```r
source("reg4_full_pipeline.R", echo = TRUE)
```

---

## Data sources

All data are public and are downloaded by the pipeline. No raw data are redistributed in this repository.

| Data | Source | Identifier |
|---|---|---|
| Batch-corrected RNA-seq (log2) | UCSC Xena, PanCanAtlas hub | `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena` |
| Clinical and survival endpoints (TCGA-CDR) | UCSC Xena, PanCanAtlas hub | `Survival_SupplementalTable_S1_20171025_xena_sp` |
| Tumour purity (ABSOLUTE) | GDC PanCanAtlas publication page | `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` |
| Molecular subtypes | UCSC Xena, PanCanAtlas hub | `TCGASubtype.20170308.tsv` (`Subtype_Selected`) |

---

## Analysis overview

| Step | Content | Main output |
|---|---|---|
| 0 | Setup: packages, gene sets, cohort groups | — |
| 1 | Download the 50 required genes | `data/expr_subset.rds` |
| 2 | Clinical data, purity; primary tumours, one per patient; ESCA split into EAC and ESCC | 9,849 tumours, 34 cohorts |
| 3 | Program scores (within-cohort z-score means) and REG4 detection | `00_REG4_detection_by_cohort.csv` |
| 4 | **Q1.** REG4 vs 13 marker genes | `01_Q1_gene_correlations.csv` |
| 5 | **Q3.** REG4 vs four program scores; dominant axis per cohort | `02_axis_correlations.csv`, `02_REG4_dominant_axis.csv` |
| 6 | **Q2.** GI vs non-GI enrichment tests; Figure 1B | `02_Q2_enrichment_tests.csv` |
| 7 | **Q4.** Mucinous histology | `03_Q4_*.csv` |
| 8 | **Q5.** Survival: does REG4 add information beyond secretory state? | `04_Q5_survival.csv` |
| 9 | **Q6.** Four-state REG4 × secretory analysis | `05_Q6_*.csv` |
| 9b | Molecular subtype (MSI/EBV) adjustment | `06_subtype_*.csv` |
| 10 | Supporting figures | `results/figures/*.pdf` |

### Key parameters

| Parameter | Value |
|---|---|
| Samples | Primary tumours only (TCGA type 01; 03 for LAML), one per patient |
| GI adenocarcinoma cohorts | COAD, READ, PAAD, STAD, ESCA_EAC, CHOL |
| Non-epithelial cohorts | LAML, DLBC, THYM, SARC, SKCM, UVM, GBM, LGG, PCPG, TGCT, MESO, UCS |
| Secretory genes | MUC2, TFF3, SPINK4, FCGBP, AGR2 |
| Epithelial genes | KRT8, KRT18, KRT19, EPCAM |
| GI lineage transcription factors | CDX2, VIL1, KLF4, HNF4A |
| Stromal/EMT genes | COL1A1, COL1A2, COL3A1, COL5A1, FN1, SPARC, POSTN, FAP, PDGFRB, THY1, DCN, LUM, VIM, ZEB1, SNAI2, TWIST1 |
| Immune/cytotoxic genes | PTPRC, CD2, CD3D, CD3E, CD8A, LCK, GZMA, GZMB, PRF1, IFNG, CXCL9, CXCL10 |
| NK genes | KLRD1, KLRF1, NCR1, NCR3, SH2D1B, KIR2DL4, KIR3DL1, NCAM1 |
| Score method | Mean of within-cohort z-scores; REG4 excluded from all scores |
| REG4 detection | log2 > 1; cohorts with < 20% detection classed as "not expressed" |
| Correlation | Spearman; partial Spearman adjusted for ABSOLUTE purity |
| Dominance rule | Top |rho| ≥ 0.3, FDR < 0.05, and a margin ≥ 0.1 over the runner-up (three axes) |
| Multiple testing | Benjamini–Hochberg FDR |
| Survival | Cox models; PFI primary, OS secondary; likelihood-ratio test for adding REG4 |

---

## Outputs

### Result tables → Supplementary Tables

| File (`results/`) | Supplementary Table |
|---|---|
| `00_REG4_detection_by_cohort.csv` | S1 |
| `01_Q1_gene_correlations.csv` | S2 |
| `02_axis_correlations.csv` | S3a |
| `02_REG4_dominant_axis.csv` | S3b |
| `02_Q2_enrichment_tests.csv` | S3c |
| `03_Q4_histology_table.csv` | S4a |
| `03_Q4_mucinous_tests.csv` | S4b |
| `03_Q4_nonmucinous_correlations.csv` | S4c |
| `04_Q5_survival.csv` | S5 |
| `05_Q6_group_sizes.csv` | S6a |
| `05_Q6_group_features.csv` | S6b |
| `05_Q6_interaction.csv` | S6c |
| `05_Q6_survival_groups.csv` | S7 |
| `06_subtype_medians.csv` | S8a |
| `06_subtype_adjusted_immune.csv` | S8b |

### Figures (`results/figures/`) → manuscript

| File | Manuscript figure |
|---|---|
| `Fig1A_REG4_expression_pancancer.pdf` | Figure 1A: REG4 expression across 34 cohorts |
| `Fig1B_REG4_dominance_heatmap.pdf` | Figure 1B: dominant REG4-associated axis |
| `Fig2_REG4_secretory_program.pdf` | Figure 2: secretory program and organ-dependent lineage coupling |
| `Fig3_mucinous_and_breast_control.pdf` | Figure 3: mucinous histology and breast negative control |
| `Fig4_molecular_subtype.pdf` | Figure 4: molecular subtype and immune associations |
| `FigS1_S2_survival_and_purity.pdf` | Figures S1 (survival) and S2 (purity robustness) |

---
