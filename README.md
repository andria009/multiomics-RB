# Transcriptomic and DNA methylation markers of retinoblastoma de-differentiation, tested against histological high-risk features in a held-out cohort

Analysis code for:

> Amiruddin PO, Mazaya M, Arisal A, Irfani I, Kartiwa A, Kartasasmita AS. Transcriptomic and DNA methylation markers of retinoblastoma de-differentiation and their association with histological high-risk features in a held-out cohort: a multi-omic meta-analysis. (submitted)

One script run reproduces every result, table and figure of the paper from public data.

## What the analysis does

1. **Expression discovery (steps 01–11).** Five public contrasts in 195 primary retinoblastomas are combined by random-effects meta-analysis (REML, Hedges' g). Genes are then filtered (F1–F7), ranked by the lower 95% CI bound of the pooled effect and checked for stability by a tumour bootstrap. The candidates are described in single cells (13 tumours), in immune context and against published subtype signatures. The five contrasts are: histology (GSE97508), expression clusters (GSE29683, GSE59983), anaplasia grade (GSE110811) and methylation subtype (GSE268136).
2. **Methylation discovery (steps 13–17).** GSE267015 (450K and EPIC arrays) is the methylation data of the GSE268136 tumours. In the tumours that also have RNA-seq, the analysis finds CpGs and regions associated with an expression-based phenotype score. It then tests whether candidate genes' promoter methylation tracks their expression, builds a methylation-supported shortlist and selects a 10-CpG panel for FFPE assays.
3. **Held-out validation (step 18).** All markers are tested against high-risk histology and against subtype in the Institut Curie cohort of Liu et al. 2021 (GSE58780 expression, GSE58783 methylation, per-tumour histology from the paper's Supplementary Data 1). This cohort is never used for selection. Two high-risk definitions were fixed before the cohort was analysed:
   - strict: postlaminar optic nerve, deep choroidal or scleral invasion;
   - broad: the strict features plus "extended" choroidal invasion.

Ryl's methylation clusters (A/B/C) were derived from the GSE267015 data. Comparing Cluster C with A+B within those data is therefore circular and is reported as descriptive only. The phenotype used in the methylation analysis is instead the expression-based phenotype score of each tumour's paired RNA-seq. Its gene sets come from the four other contrasts, so it uses neither the methylation data nor the GSE268136 labels.

## Main results

| | |
|---|---|
| Expression candidates / passing F1–F7 | 3,298 / 497 |
| Expression shortlist (core genes in bold) | **SSTR2**, RCOR2, NKAIN1, BCL11A, EBF3 (up); **CAMK1D**, ARL6IP5, ANKRD33, SORBS2, CA2 (down) |
| Exploratory, subtype-independent tier | MYLK, DNM3 |
| Methylation (47 paired tumours) | 19,836 DMPs (4,491 higher, 15,345 lower in aggressive tumours); 4,543 DMRs |
| Promoter methylation coupled to expression | 5.5% of candidate genes vs 2.3% of other genes (OR 2.5) |
| Multi-omic shortlist | GRIK3, SMAD9, TFF1, C20orf27, TSPAN11 (up); SV2B, PARD6B, PLEKHH1, KLF3, PROM1 (down) |
| CpG panel (10 sites) | cross-validated ρ 0.71 with the phenotype score; ρ 0.92 with Liu's 9-CpG predictor |
| Liu/Curie cohort, subtype 2 | AUC 0.87–1.00 (phenotype score, shortlisted genes, CpG panel) |
| Liu/Curie cohort, strict high-risk histology | AUC 0.51–0.71 (expression), 0.57 (CpG panel); none significant after BH correction; the SSTR2 association (AUC 0.71) is explained by subtype; DNM3 subtype-adjusted OR 2.3 (0.86–6.2) |

The malignant-cell calls are stable between runs on the same machine, but the CopyKAT cross-check varies between runs (10–12 of 13 tumours separated in repeated runs). The shortlists do not depend on it.

## Repository layout

```
code/                  analysis scripts, numbered in run order; run_all.sh runs them all
code/R/                config.R (paths, parameters), utils.R (helpers), gene_lists.R
data/download_raw_data.sh   fetches all GEO files into data/raw/
data/supplementary/    published supplementary tables (download by hand, see below)
data/cibersortx/       optional CIBERSORTx web-tool results (see below)
results/               all result tables (created by the run); results/supporting_info/ = S1–S12 Tables
figures/               diagnostic figures from the R scripts (created by the run)
logs/run_all.log       full log of the run
```

## Pipeline

| Step | Script | What it does | Output |
|---|---|---|---|
| 01 | `01_download_geo.R` | Loads (or downloads) the four microarray cohorts | `results/01_bulk/` |
| 02 | `02_normalize_microarray.R` | Primary tumours only; probe-to-gene collapse; quantile normalisation | `results/01_bulk/` |
| 03 | `03_rnaseq_gse268136.R` | GSE268136 gene matrix; Cluster C vs A+B (limma) | `results/01_bulk/` |
| 04 | `04_invasion_meta_analysis.R` | Per-contrast limma and Hedges' g; REML meta-analysis, HKSJ, prediction intervals, Stouffer. Also run twice with `INVASION_EXCLUDE` for the sensitivity analyses | `results/02_invasion*/` |
| 05 | `05_invasion_signature.R` | Exploratory LASSO classifier (not used for selection) | `results/02_invasion/` |
| 06 | `06_scrna_integration.R` | scRNA-seq QC, Harmony integration, cell types, arm-level CNV malignant calls | `results/03_scrna/` |
| 06b | `06b_copykat_validation.R` | CopyKAT check of the malignant calls | `results/03_scrna/` |
| 07 | `07_scrna_validation.R` | Per-tumour expression in malignant and tumour-lineage cells | `results/02_invasion/`, `results/03_scrna/` |
| 08 | `08_immune_deconvolution.R` | MCP-counter, ESTIMATE; optional CIBERSORTx summary | `results/04_immune/` |
| 09 | `09_invasion_immune.R` | Leave-one-dataset-out phenotype score; immune co-variation | `results/05_invasion_immune/` |
| 10 | `10_subtype_benchmark.R` | Liu subtype-2 and Kooi photoreceptorness scores; subtype-adjusted effects | `results/06_subtype/` |
| 11 | `11_priority_ranking.R` | Filters F1–F7, lower-bound ranking, bootstrap, alternative scenarios | `results/07_priority/` |
| 13 | `13_methylation_preprocess.R` | GSE267015 beta matrix, tumour selection, probe filters, annotation, PCA | `results/08_methylation/` |
| 14 | `14_methylation_phenotype.R` | DMPs and DMRs associated with the phenotype score; Cluster C vs A+B (circular, descriptive) | `results/08_methylation/` |
| 15 | `15_methylation_expression.R` | Promoter and gene-body methylation vs expression and phenotype score; methylation support for candidates | `results/09_multiomic/` |
| 16 | `16_cpg_marker_panel.R` | 10-CpG panel: ranking, bootstrap, cross-validation with selection inside each fold | `results/08_methylation/` |
| 17 | `17_multiomic_priority.R` | Expression + methylation evidence for genes passing F1–F7; multi-omic shortlist | `results/09_multiomic/` |
| 18 | `18_liu_validation.R` | Held-out Liu/Curie cohort: all markers vs high-risk histology (with and without subtype adjustment) and vs subtype; subtype vs histology | `results/10_liu_validation/` |
| 19 | `19_supporting_info.R` | S1–S12 Tables, package versions, session info | `results/supporting_info/` |

There is no step 12. `00_scout_methylation_cohorts.R` is a one-off inspection of GSE57362 and GSE58783 and is not part of the pipeline.

## Data

All expression, methylation and single-cell data are public on GEO. `data/download_raw_data.sh` fetches them into `data/raw/` (about 7 GB). It skips complete files and re-downloads truncated `.gz` files.

Three published supplementary tables are read from `data/supplementary/` and must be downloaded by hand from the journals:

| File | Source | Used for |
|---|---|---|
| `41467_2021_25792_MOESM4_ESM.xlsx` | Liu et al. 2021, *Nat Commun* 12:5578, Supplementary Data 1 | Liu/Curie histology, subtype assignment, 9-CpG predictor (step 18) |
| `41467_2021_25792_MOESM6_ESM.xlsx` | Liu et al. 2021, Supplementary Data 3 | Liu subtype-2 signature (step 10) |
| `mmc10.xlsx` | Kooi et al. 2015, *EBioMedicine* 2:660–670, supplementary table of gene clusters | Kooi photoreceptorness signature (step 10) |

CIBERSORTx (optional, reported in Methods only): step 08 writes mixture files to `results/04_immune/cibersortx_input/`. Run them on https://cibersortx.stanford.edu (LM22, 100 permutations) and save the results as `data/cibersortx/CIBERSORTx_<GSE>.txt`. Without these files the summary is skipped.

## Setup (once)

R 4.6.1 from CRAN (https://cran.r-project.org) was used. In R:

```r
install.packages(c("BiocManager", "remotes", "here", "dplyr", "tidyr", "tibble", "purrr", "readr", "ggplot2",
                   "ggpubr", "ggrepel", "circlize", "glmnet", "pROC", "readxl", "openxlsx", "data.table",
                   "Matrix", "Seurat", "harmony"))
BiocManager::install(c("GEOquery", "Biobase", "limma", "preprocessCore", "AnnotationDbi", "org.Hs.eg.db",
                       "ComplexHeatmap", "minfi", "DMRcate", "IlluminaHumanMethylation450kanno.ilmn12.hg19"))
remotes::install_github(c("navinlabcode/copykat", "omnideconv/immunedeconv"))
```

The exact versions used are in `results/supporting_info/package_versions.tsv`.

## Run

From the repository root:

```bash
bash data/download_raw_data.sh   # once
bash code/run_all.sh             # about 1 h; stops at the first failing step; log in logs/run_all.log
```

The random seed is 42 throughout.

## Citation and licence

Please cite the paper above and the data sets' original publications. The code is released under the MIT licence (`LICENSE`).
