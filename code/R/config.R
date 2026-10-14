# =============================================================================
# config.R — paths and parameters (transcriptome 01-11, methylation 13-17, validation 18, copy number 20, MOFA 21)
# Sourced at the top of every script:  source(here::here("code", "R", "config.R"))
# The project root is marked by the empty file `.here`.
# =============================================================================

library(here)

# -----------------------------------------------------------------------------
# Directories
# -----------------------------------------------------------------------------
PROJ_ROOT      <- here::here()
DIR_RAW_MICRO  <- here::here("data", "raw", "microarray")
DIR_RAW_RNASEQ <- here::here("data", "raw", "rnaseq")
DIR_RAW_SCRNA  <- here::here("data", "raw", "scrna")   # RAW.tar files; extracted in place by step 06
DIR_RAW_METH   <- here::here("data", "raw", "methylation")
DIR_RAW_LIU    <- here::here("data", "raw", "liu")          # GSE58780 / GSE58783 validation cohort

DIR_RES  <- here::here("results")
DIR_FIG  <- here::here("figures")
DIR_BULK <- file.path(DIR_RES, "01_bulk")            # steps 01-03: normalised matrices, RNA-seq contrast
DIR_INV  <- file.path(DIR_RES, "02_invasion")        # steps 04, 05, 07: meta-analysis, signature, scRNA check
DIR_SC   <- file.path(DIR_RES, "03_scrna")           # step 06: integrated single-cell object
DIR_IMM  <- file.path(DIR_RES, "04_immune")          # step 08: MCP-counter deconvolution
DIR_INT  <- file.path(DIR_RES, "05_invasion_immune") # step 09: phenotype score, immune co-variation
DIR_SUB  <- file.path(DIR_RES, "06_subtype")         # step 10: published subtype / photoreceptorness scores
DIR_PRI  <- file.path(DIR_RES, "07_priority")        # step 11: priority ranking
DIR_METH <- file.path(DIR_RES, "08_methylation")     # steps 13-14, 16: methylation preprocessing, DMPs, CpG panel
DIR_MO   <- file.path(DIR_RES, "09_multiomic")       # steps 15, 17: methylation-expression integration, priority
DIR_VAL  <- file.path(DIR_RES, "10_liu_validation")  # step 18: held-out Liu/Curie cohort, histology
DIR_SI   <- file.path(DIR_RES, "supporting_info")    # step 19: S tables + session info
DIR_CNV  <- file.path(DIR_RES, "11_cnv")             # step 20: copy number from the methylation arrays
DIR_MOFA <- file.path(DIR_RES, "12_mofa")            # step 21: multi-omics factor analysis

# Sensitivity analysis: INVASION_EXCLUDE=GSE29683 drops that dataset from the
# invasion meta-analysis and training (steps 04-05) and writes to separate
# result/figure folders so the main results are not overwritten.
INVASION_EXCLUDE <- Filter(nzchar, strsplit(Sys.getenv("INVASION_EXCLUDE", ""), ",")[[1]])
if (length(INVASION_EXCLUDE) > 0) {
  DIR_INV <- paste0(DIR_INV, "_excl_", paste(INVASION_EXCLUDE, collapse = "_"))
  DIR_FIG <- file.path(DIR_FIG, basename(DIR_INV))
}

for (d in c(DIR_BULK, DIR_INV, DIR_SC, DIR_IMM, DIR_INT, DIR_SUB, DIR_PRI, DIR_METH, DIR_MO, DIR_VAL, DIR_SI, DIR_CNV, DIR_MOFA, DIR_FIG)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Datasets
# -----------------------------------------------------------------------------
GEO_INVASION_MICRO <- c("GSE97508", "GSE29683", "GSE59983", "GSE110811")
GEO_RNASEQ         <- "GSE268136"   # Ryl 2024; log2(TPM+1) only, no counts in GEO
GEO_SCRNA          <- c("GSE168434", "GSE249995", "GSE166173")
GEO_METH           <- "GSE267015"   # Ryl 2024 450K + EPIC methylation; same tumours as GSE268136

# Primary tumours only (step 02): pattern matched against the GEO title/characteristics.
# Drops normal/adult retina, retinocytoma, cell lines and xenografts.
PRIMARY_PATTERN <- c(GSE97508  = "tissue: (non-)?invasive retinoblastoma",
                     GSE29683  = "cell type: primary tumor",
                     GSE59983  = "tissue: primary Rb tissue",
                     GSE110811 = "tissue: Retinoblastoma")

# What each contrast's label actually measures (reported in Table 1 and S1 Table)
LABEL_TYPE <- c(GSE97508  = "histology (submitter sample titles)",
                GSE29683  = "expression cluster (this study: k-means on marker score)",
                GSE59983  = "expression cluster (Kooi 2015 unsupervised clustering)",
                GSE110811 = "anaplasia grade (Hudson 2018)",
                GSE268136_ClusterC = "methylation cluster (Ryl 2024)")
EXPRESSION_LABELLED <- c("GSE29683", "GSE59983")

# -----------------------------------------------------------------------------
# Thresholds and parameters
# -----------------------------------------------------------------------------
DEG_LOG2FC_CUTOFF <- 1.0     # RNA-seq summary count only (step 03)
DEG_ADJP_CUTOFF   <- 0.05
INV_LOG2FC_CUTOFF <- 0.585   # per-contrast volcano lines only (step 04)
INV_ADJP_CUTOFF   <- 0.05
RNASEQ_MIN_TPM    <- 1       # GSE268136 filter: TPM >= 1 in >= smallest-group samples

# Random-effects meta-analysis (Hedges' g, REML tau^2, Wald z for discovery)
META_ADJP   <- 0.05   # BH-adjusted p of the pooled effect
META_G_MIN  <- 0.5    # |pooled g| (medium standardised effect)

# scRNA-seq
SCRNA_MIN_CELLS    <- 3
SCRNA_MIN_FEATURES <- 200
SCRNA_MAX_FEATURES <- 6000
SCRNA_MAX_MT_PCT   <- 20
SCRNA_N_DIMS       <- 30
SCRNA_RESOLUTION   <- 0.5
SCRNA_N_HVG        <- 3000

# Exploratory classifier (step 05)
ML_CV_FOLDS <- 10
LASSO_ALPHA <- 1.0   # 1 = LASSO
N_BOOT      <- 200   # bootstrap resamples for stability selection and AUC CIs

# Priority filters (step 11); see the step 11 header for which choices were revised and why
PRI_MIN_CONTRASTS <- 4   # measured in >= 4 of 5 contrasts
PRI_MIN_AGREE     <- 4   # same sign as the pooled effect in >= 4 contrasts
PRI_N_PER_DIR     <- 5   # shortlist: top 5 up + top 5 down by the lower 95% CI bound of |g|
PRI_SC_MIN_PCT    <- 10  # F7: detected in >= 10% of malignant cells ...
PRI_SC_MIN_TUMOURS <- 3  #     ... in >= 3 single-cell tumours (added in revision 2)
PRI_CORE_PROB     <- 0.5 # core genes: in the top PRI_N_PER_DIR in >= 50% of bootstrap resamples
SENS_EXCLUDES     <- c("GSE29683", "GSE29683,GSE59983")  # sensitivity runs of step 04

# Methylation (steps 13-17)
METH_MAX_NA       <- 0      # drop CpGs with any missing beta value
METH_BETA_CLIP    <- 0.001  # beta clipped to [clip, 1 - clip] before the M-value transform
METH_DMP_ADJP     <- 0.05   # BH threshold for phenotype-associated CpGs
METH_DMP_DB_MIN   <- 0.05   # |delta beta per SD of phenotype score| for a DMP
METH_RHO_MIN      <- 0.3    # |Spearman rho| for methylation-expression coupling (step 15)
METH_PROMOTER     <- c("TSS1500", "TSS200", "5'UTR", "1stExon")   # promoter CpG groups
CPG_PANEL_N       <- 5      # CpG panel: top regions per direction (hyper, hypo) (step 16)
CPG_MIN_RANGE     <- 0.3    # CpG panel: beta range (10th-90th percentile) across tumours
CPG_MERGE_BP      <- 1000   # CpGs within this distance and correlated are one region
CPG_MERGE_R       <- 0.8

# Validation cohort (step 18). Both definitions were fixed before the cohort was analysed; in the
# revision the broad definition was made primary because, read with the IRSWG consensus, "extended"
# choroidal invasion corresponds to massive choroidal invasion (a standard high-risk feature).
LIU_SUPP_FILE  <- "41467_2021_25792_MOESM4_ESM.xlsx"   # Liu et al. 2021 Supplementary Data 1 (Table S1)
HR_STRICT_ON   <- "postlaminar"                          # optic nerve: high risk
HR_STRICT_CH   <- c("choroid_deep", "choroid_extended_and_sclera")   # choroid/sclera: strict (secondary)
HR_BROAD_CH    <- c(HR_STRICT_CH, "choroid_extended")    # broad (primary since the revision)
HR_POWER_AUC   <- c(0.65, 0.70, 0.75, 0.80)              # true AUCs for the power table (step 18)

# Copy number from methylation arrays (step 20; conumee with CopyNeutralIMA normal references)
CNV_ARM_THR   <- 0.10   # arm called gained / lost when its length-weighted segment mean is beyond +/- this (log2 ratio)
CNV_SEG_THR   <- 0.15   # segments beyond +/- this count towards the fraction of genome altered (FGA)
CNV_AMP_THR   <- 0.50   # focal amplification call for MYCN (mean of overlapping bins)
CNV_MAX_NOISE <- 1.00   # arrays with conumee fit noise above this are excluded (QC)
CNV_KEY_ARMS  <- c("1q", "2p", "6p", "16q")   # recurrent RB copy-number changes
# Multi-omics factor analysis (step 21; MOFA2)
MOFA_N_GENES   <- 2000  # most variable expressed genes (RNA view)
MOFA_N_CPGS    <- 5000  # most variable CpGs (methylation view)
MOFA_N_FACTORS <- 10
MOFA_SEEDS     <- 1:5   # extra seeds for the stability check of the main factor

# Round-2 revision
COPYKAT_SEEDS  <- c(42, 1:4)   # step 06b: CopyKAT repeated with these seeds (first = reported run)
METH_EPIC_ARTEFACT_CHR <- c(paste0("chr", 3:8), paste0("chr", 17:22))   # step 23: chromosomes with the EPIC copy-number artefact

# -----------------------------------------------------------------------------
# Reproducibility
# -----------------------------------------------------------------------------
GLOBAL_SEED <- 42
set.seed(GLOBAL_SEED)

# -----------------------------------------------------------------------------
# Plot defaults
# -----------------------------------------------------------------------------
PLOT_WIDTH_WIDE <- 12
PLOT_HEIGHT     <- 8
PLOT_DPI        <- 300
PLOT_THEME      <- ggplot2::theme_bw(base_size = 12)

COL_UPREGULATED   <- "#E41A1C"
COL_DOWNREGULATED <- "#377EB8"
COL_NONSIG        <- "grey70"
COL_INVASIVE      <- "#C0392B"
COL_NONINVASIVE   <- "#2980B9"
COL_NORMAL        <- "#27AE60"

message("[config.R] Configuration loaded. Project root: ", PROJ_ROOT)
