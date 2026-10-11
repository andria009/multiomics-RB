# =============================================================================
# Step 03 — GSE268136 RNA-seq (Ryl et al. 2024): gene matrix + Cluster C vs A+B
# Data  : log2(TPM+1) of sleuth-normalised TPM (GEO supplementary; the file is
#         named "log1ptpm" but values are base 2). Raw counts are not in GEO, so
#         limma (trend, robust) is used on log-scale data.
# Steps : 1. Phenotype (Ryl cluster, RB1 status) from the GPL24676 series matrix
#         2. Transcripts -> genes (sum of transcript TPMs, then log2(TPM+1))
#         3. Align matrix samples to phenotype by title
#         4. limma: Cluster C (MYCN-driven) vs Cluster A+B, RB1-wild-type excluded
#            (including RB_13, flagged by the submitters as possible retina
#            contamination), genes kept if TPM >= 1 in >= the smaller group size
# Output: results/01_bulk/rnaseq_log2tpm_GSE268136.rds, rnaseq_pheno_GSE268136.rds,
#         limma_clusterC_vs_AB_GSE268136.rds/.tsv
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))

load_pkg("limma",    bioc = TRUE)
load_pkg("GEOquery", bioc = TRUE)
load_pkg("dplyr")
load_pkg("readr")
load_pkg("tibble")

# =============================================================================
# 1. Load pheno (58 primary RBs from GPL24676 sub-series)
# =============================================================================
message("Loading GSE268136 pheno ...")

pheno_cache <- file.path(DIR_RAW_RNASEQ, "GSE268136_rnaseq.rds")
if (!file.exists(pheno_cache)) {
  # Rebuild the cache from the local GPL24676 series matrix (data/download_raw_data.sh)
  sm <- file.path(DIR_RAW_RNASEQ, "GSE268136-GPL24676_series_matrix.txt.gz")
  if (!file.exists(sm)) stop("Missing: ", sm, "\nRun data/download_raw_data.sh first.")
  eset <- GEOquery::getGEO(filename = sm, getGPL = FALSE)
  saveRDS(list(gse_id = "GSE268136", pheno = Biobase::pData(eset)), pheno_cache)
}
pheno <- readRDS(pheno_cache)$pheno
message("  Pheno samples: ", nrow(pheno))

# =============================================================================
# 2. Annotate Ryl clusters from characteristics
#    Raw value: "cluster ryl_et_al: A" / "cluster ryl_et_al: B" / "cluster ryl_et_al: C"
#    WT samples: "genotype: WT" in another characteristics column → exclude
# =============================================================================
char_cols <- grep("characteristics_ch", names(pheno), value = TRUE)

# Extract cluster letter directly from "cluster ryl_et_al: X" pattern
pheno$ryl_cluster <- NA_character_
for (col in char_cols) {
  vals <- as.character(pheno[[col]])
  extracted <- sub(".*cluster ryl_et_al:\\s*([ABC]).*", "\\1",
                   vals, ignore.case = TRUE, perl = TRUE)
  # Only update rows where extraction succeeded (extracted != original string)
  updated <- (extracted != vals) & !is.na(extracted)
  pheno$ryl_cluster[updated] <- extracted[updated]
}

# Flag RB1-wild-type tumours (separate biology, mostly MYCN-amplified) and the
# "genotype: WT" samples; "RB1wt/potential_retina_contamination" (RB_13) is caught too.
pheno$is_wt <- apply(
  sapply(char_cols, function(col) grepl("genotype: *WT|rb1 status: *RB1wt", pheno[[col]],
                                         ignore.case = TRUE)),
  1, any
)

message("Ryl cluster distribution (before WT exclusion):")
print(table(pheno$ryl_cluster, useNA = "always"))
message("WT samples: ", sum(pheno$is_wt))

# =============================================================================
# 3. Load log2(TPM+1) matrix from GEO supplementary
#    File: GSE268136_ryl_et_al_rb_tumors_log1ptpm.txt.gz
#    Columns: target_id (ENST), sample columns, ext_gene (gene symbol)
# =============================================================================
message("\nLoading log2(TPM+1) matrix ...")

tpm_file <- file.path(DIR_RAW_RNASEQ, "GSE268136_ryl_et_al_rb_tumors_log1ptpm.txt.gz")
if (!file.exists(tpm_file))
  stop("Missing: ", tpm_file)

tpm_raw <- readr::read_tsv(tpm_file, show_col_types = FALSE)
message("  Raw matrix: ", nrow(tpm_raw), " transcripts x ",
        ncol(tpm_raw) - 2, " samples")  # -2 for target_id and ext_gene

# Separate annotation from expression
enst_ids  <- tpm_raw$target_id
ext_gene  <- tpm_raw$ext_gene
sample_cols <- setdiff(colnames(tpm_raw), c("target_id", "ext_gene"))

expr_mat  <- as.matrix(tpm_raw[, sample_cols])
rownames(expr_mat) <- enst_ids

message("  Samples in matrix: ", length(sample_cols))
message("  Genes with ext_gene annotation: ",
        sum(!is.na(ext_gene) & ext_gene != ""))

# =============================================================================
# 4. Collapse ENST transcripts → gene symbols
#    Gene TPM is the SUM of its transcripts' TPMs; averaging log values instead
#    would pull down genes with many lowly expressed isoforms.
# =============================================================================
message("\nCollapsing transcripts to gene level ...")

# Use ext_gene; drop rows without gene annotation
has_gene  <- !is.na(ext_gene) & ext_gene != ""
expr_gene <- expr_mat[has_gene, , drop = FALSE]
gene_ids  <- ext_gene[has_gene]

# Aggregate by gene: sum transcript TPMs on the linear scale, then back to log2(TPM+1).
# Negative TPMs (sleuth normalisation artefacts) are treated as 0.
gene_tpm <- rowsum(pmax(2^expr_gene - 1, 0), gene_ids, na.rm = TRUE)
gene_mat <- log2(gene_tpm + 1)
stopifnot(!anyNA(gene_mat), min(gene_mat) >= 0)

# Remove genes with zero variance (uninformative for limma)
gene_sd  <- apply(gene_mat, 1, sd, na.rm = TRUE)
gene_mat <- gene_mat[gene_sd > 0.01, , drop = FALSE]

message("  Unique genes after collapse: ", nrow(gene_mat))

# =============================================================================
# 5. Align matrix samples to pheno via title matching
#    Matrix columns: "RB_3", "INF_R_1024", etc.
#    Pheno titles:   "RB_3", "RB_4", etc.
# =============================================================================
message("\nAligning samples ...")

mat_names   <- colnames(gene_mat)
title_match <- match(mat_names, pheno$title)
matched_idx <- which(!is.na(title_match))

if (length(matched_idx) < 10)
  stop("Fewer than 10 samples matched between matrix and pheno. Check title format.")

gene_mat_sub <- gene_mat[, matched_idx, drop = FALSE]
pheno_sub    <- pheno[title_match[matched_idx], , drop = FALSE]
colnames(gene_mat_sub) <- rownames(pheno_sub)

message("  Matched samples: ", ncol(gene_mat_sub),
        " of ", ncol(gene_mat), " matrix samples")
message("  Unmatched matrix samples: ",
        paste(mat_names[is.na(title_match)], collapse = ", "))

# Cluster distribution in matched set
message("Matched sample cluster distribution:")
print(table(pheno_sub$ryl_cluster, useNA = "always"))
message("Matched WT samples: ", sum(pheno_sub$is_wt))

# =============================================================================
# 6. limma DEG: Cluster C vs Cluster A+B
#    Exclude WT samples (separate biology — RB1wt/MYCN-amplified)
# =============================================================================
message("\n=== limma: Cluster C vs Cluster A+B (non-WT samples) ===")

pheno_sub$cluster_comp <- dplyr::case_when(
  pheno_sub$ryl_cluster == "C" & !pheno_sub$is_wt  ~ "ClusterC",
  pheno_sub$ryl_cluster %in% c("A", "B") & !pheno_sub$is_wt ~ "ClusterAB",
  TRUE ~ NA_character_
)

keep_clust  <- !is.na(pheno_sub$cluster_comp)
mat_clust   <- gene_mat_sub[, keep_clust, drop = FALSE]
pheno_clust <- pheno_sub[keep_clust, , drop = FALSE]

n_c  <- sum(pheno_clust$cluster_comp == "ClusterC")
n_ab <- sum(pheno_clust$cluster_comp == "ClusterAB")
message("  Cluster C (MYCN-driven): n = ", n_c)
message("  Cluster A+B:             n = ", n_ab)

if (n_c < 2 || n_ab < 2)
  stop("Insufficient samples in one group for limma. Check cluster annotation.")

# Remove genes with any NA, then low-expression genes (limma-trend on unfiltered
# log-TPM ranks noisy low-count pseudogenes at the top)
mat_clust <- mat_clust[complete.cases(mat_clust), , drop = FALSE]
min_n     <- min(n_c, n_ab)
expressed <- rowSums((2^mat_clust - 1) >= RNASEQ_MIN_TPM) >= min_n
message("  Expressed genes (TPM >= ", RNASEQ_MIN_TPM, " in >= ", min_n, " samples): ",
        sum(expressed), " of ", nrow(mat_clust))
mat_clust <- mat_clust[expressed, , drop = FALSE]

pheno_clust$cluster_comp <- factor(pheno_clust$cluster_comp,
                                    levels = c("ClusterAB", "ClusterC"))
design <- model.matrix(~ cluster_comp, data = pheno_clust)

fit <- limma::lmFit(mat_clust, design)
fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

res_clusterC <- limma::topTable(fit, coef = "cluster_compClusterC",
                                 number = Inf, sort.by = "P") |>
  tibble::rownames_to_column("gene") |>
  dplyr::rename(log2FC = logFC) |>
  dplyr::arrange(adj.P.Val)
# Ordinary (unmoderated) t for Hedges' g in step 04; moderated t is kept for testing
t_ord <- fit$coefficients[, "cluster_compClusterC"] / (fit$stdev.unscaled[, "cluster_compClusterC"] * fit$sigma)
res_clusterC$t_ord <- unname(t_ord[res_clusterC$gene])
attr(res_clusterC, "n1") <- n_c    # group sizes for effect sizes (step 04)
attr(res_clusterC, "n2") <- n_ab

n_sig <- sum(abs(res_clusterC$log2FC) >= DEG_LOG2FC_CUTOFF &
             res_clusterC$adj.P.Val   <= DEG_ADJP_CUTOFF, na.rm = TRUE)
message("  Significant DEGs (|log2FC|>=", DEG_LOG2FC_CUTOFF,
        ", padj<=", DEG_ADJP_CUTOFF, "): ", n_sig)

message("\n  Top 10 upregulated in Cluster C:")
print(head(res_clusterC |> dplyr::filter(log2FC > 0) |>
             dplyr::select(gene, log2FC, adj.P.Val), 10))
message("\n  Top 10 downregulated in Cluster C:")
print(head(res_clusterC |> dplyr::filter(log2FC < 0) |>
             dplyr::select(gene, log2FC, adj.P.Val), 10))

# =============================================================================
# 7. Save
# =============================================================================
# Expression matrix (gene-level log2(TPM+1), all matched samples including WT).
save_rds(gene_mat_sub,  file.path(DIR_BULK, "rnaseq_log2tpm_GSE268136.rds"),
         "GSE268136 gene-level log2(TPM+1) matrix")
save_rds(pheno_sub,     file.path(DIR_BULK, "rnaseq_pheno_GSE268136.rds"),
         "GSE268136 pheno (matched to matrix)")
save_rds(res_clusterC,  file.path(DIR_BULK, "limma_clusterC_vs_AB_GSE268136.rds"),
         "limma Cluster C vs A+B (GSE268136)")
save_tsv(res_clusterC,  file.path(DIR_BULK, "limma_clusterC_vs_AB_GSE268136.tsv"))

message("\n[03_rnaseq_gse268136.R] Complete.")
message("  Cluster C (n=", n_c, ") vs A+B (n=", n_ab, ") | Sig DEGs: ", n_sig)
