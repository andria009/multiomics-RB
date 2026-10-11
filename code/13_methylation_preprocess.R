# =============================================================================
# Step 13 — DNA methylation: GSE267015 (Ryl et al. 2024) preprocessing
# Input : data/raw/methylation/GSE267015_ryl_et_al_beta_matrix_tumors.txt.gz
#           (submitter-processed beta values, probes common to 450K and EPIC)
#         data/raw/methylation/GSE267015-GPL13534/-GPL21145_series_matrix.txt.gz
#           (sample annotation of the 450K and EPIC sub-series)
#         results/01_bulk/rnaseq_pheno_GSE268136.rds (step 03; RNA-seq of the same tumours)
# Output: results/08_methylation/
#   meth_samples.tsv          every sample: platform, Ryl cluster, RB1/MYCN status, used?, paired RNA
#   meth_beta.rds             filtered beta matrix (CpGs x used tumours)
#   meth_probe_annotation.rds probe annotation (hg19): chr, position, island context, genes and regions
#   meth_probe_gene.rds       long table probe - gene - region (one row per pair)
#   figures/13_meth_pca
# Samples used: primary tumours with an RB1 mutation or RB1 mutation not detected,
# i.e. the same rule as the RNA-seq contrast (RB1-wild-type tumours, relapses and cell
# lines excluded; RB_13 excluded as possible retina contamination, as in step 03).
# Probes: autosomal; SNP-affected (MAF >= 0.05 within 2 bp of the CpG) and cross-hybridising
# probes removed (DMRcate::rmSNPandCH); CpGs with missing values removed.
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))

load_pkg("data.table", bioc = FALSE)
load_pkg("dplyr")
load_pkg("ggplot2")
load_pkg("DMRcate", bioc = TRUE)
load_pkg("IlluminaHumanMethylation450kanno.ilmn12.hg19", bioc = TRUE)

beta_file <- file.path(DIR_RAW_METH, "GSE267015_ryl_et_al_beta_matrix_tumors.txt.gz")
mat_files <- file.path(DIR_RAW_METH, paste0("GSE267015-", c("GPL13534", "GPL21145"), "_series_matrix.txt.gz"))
for (f in c(beta_file, mat_files)) if (!file.exists(f)) stop("Missing ", f, " (run data/download_raw_data.sh)")

# -----------------------------------------------------------------------------
# 1. Sample annotation from the two sub-series matrices (no data tables needed)
# -----------------------------------------------------------------------------
pheno <- dplyr::bind_rows(lapply(mat_files, read_series_pheno))[, c("title", "gsm", "platform", "source",
                                                               "cluster_ryl_et_al", "mycn_amplification", "rb1_status")]
pheno$platform <- ifelse(pheno$platform == "GPL13534", "450K", ifelse(pheno$platform == "GPL21145", "EPIC", pheno$platform))
message("Sub-series samples: ", nrow(pheno)); print(table(pheno$source, pheno$platform))

# -----------------------------------------------------------------------------
# 2. Beta matrix; match columns to sample titles
# -----------------------------------------------------------------------------
dt <- data.table::fread(beta_file, data.table = FALSE)
beta <- as.matrix(dt[, -1]); rownames(beta) <- dt[[1]]; rm(dt); invisible(gc())
title_of <- function(col) ifelse(grepl("^INF_R_", col), sub("^INF_R_", "RB_INF_", col), col)
colnames(beta) <- title_of(colnames(beta))
unmatched <- setdiff(colnames(beta), pheno$title)
if (length(unmatched)) stop("Beta-matrix columns without annotation: ", paste(unmatched, collapse = ", "))
message("Beta matrix: ", nrow(beta), " CpGs x ", ncol(beta), " samples")

# -----------------------------------------------------------------------------
# 3. Sample selection (same rule as the RNA-seq contrast) and pairing with RNA-seq
# -----------------------------------------------------------------------------
p_rna <- readRDS(file.path(DIR_BULK, "rnaseq_pheno_GSE268136.rds"))
pheno <- pheno |>
  dplyr::mutate(in_matrix = title %in% colnames(beta),
                primary   = grepl("tumor, primary", source),
                rb1_wt    = grepl("^RB1wt", rb1_status),
                used      = in_matrix & primary & !rb1_wt & title != "RB_13",
                rnaseq_sample = rownames(p_rna)[match(title, p_rna$title)],
                paired_rna = !is.na(rnaseq_sample))
print(table(used = pheno$used, platform = pheno$platform, cluster = pheno$cluster_ryl_et_al))
message("Used tumours: ", sum(pheno$used), " (", sum(pheno$used & pheno$paired_rna), " with RNA-seq)")
save_tsv(pheno, file.path(DIR_METH, "meth_samples.tsv"))
beta <- beta[, pheno$title[pheno$used], drop = FALSE]

# -----------------------------------------------------------------------------
# 4. Probe filtering
# -----------------------------------------------------------------------------
n0 <- nrow(beta)
beta <- beta[rowSums(is.na(beta)) <= METH_MAX_NA, , drop = FALSE]
n_na <- nrow(beta)
clip <- function(b) pmin(pmax(b, METH_BETA_CLIP), 1 - METH_BETA_CLIP)
M <- log2(clip(beta) / (1 - clip(beta)))
M <- DMRcate::rmSNPandCH(M, dist = 2, mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY = TRUE)
beta <- beta[rownames(M), , drop = FALSE]
message("CpGs: ", n0, " -> ", n_na, " without missing values -> ", nrow(beta),
        " after removing SNP, cross-hybridising and sex-chromosome probes")

# -----------------------------------------------------------------------------
# 5. Annotation (hg19, Illumina 450K manifest; all retained probes are on both arrays)
# -----------------------------------------------------------------------------
ann <- as.data.frame(minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19))
ann <- ann[rownames(beta), c("chr", "pos", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Group")]
stopifnot(!anyNA(ann$chr))
ann$cpg <- rownames(ann)
gn <- strsplit(ann$UCSC_RefGene_Name, ";"); rg <- strsplit(ann$UCSC_RefGene_Group, ";")
stopifnot(all(lengths(gn) == lengths(rg)))
pg <- data.frame(cpg = rep(ann$cpg, lengths(gn)), gene = unlist(gn), region = unlist(rg)) |> dplyr::distinct()
pg$promoter <- pg$region %in% METH_PROMOTER
message("Probe-gene pairs: ", nrow(pg), " (", length(unique(pg$gene)), " genes; ",
        length(unique(pg$gene[pg$promoter])), " with promoter CpGs)")

# -----------------------------------------------------------------------------
# 6. QC: PCA of the most variable CpGs, coloured by platform and Ryl cluster
# -----------------------------------------------------------------------------
smp <- pheno[match(colnames(beta), pheno$title), ]
top <- head(order(apply(M[rownames(beta), ], 1, stats::var), decreasing = TRUE), 10000)
pc <- prcomp(t(M[rownames(beta)[top], ]), scale. = FALSE)
pcs <- data.frame(title = colnames(beta), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                  platform = smp$platform, cluster = smp$cluster_ryl_et_al)
ve <- round(100 * summary(pc)$importance[2, 1:2], 1)
r2_platform <- summary(stats::lm(pc$x[, 1] ~ smp$platform))$r.squared
message("PC1 (", ve[1], "%) variance explained by platform: R2 = ", round(r2_platform, 2))
p <- ggplot2::ggplot(pcs, ggplot2::aes(PC1, PC2, colour = cluster, shape = platform)) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"),
                title = "GSE267015: methylation PCA (10,000 most variable CpGs)") + PLOT_THEME
save_plot(p, file.path(DIR_FIG, "13_meth_pca"), width = 7, height = 5.5)

save_rds(beta, file.path(DIR_METH, "meth_beta.rds"), "filtered beta matrix")
save_rds(ann, file.path(DIR_METH, "meth_probe_annotation.rds"), "probe annotation")
save_rds(pg, file.path(DIR_METH, "meth_probe_gene.rds"), "probe-gene-region table")
save_tsv(pcs, file.path(DIR_METH, "meth_pca.tsv"))
message("\n[13_methylation_preprocess.R] Complete.")
