# =============================================================================
# Step 14 — CpG methylation associated with the aggressive (de-differentiated) phenotype
# Input : results/08_methylation/ (step 13), results/05_invasion_immune/invasion_scores.tsv
#         (step 09), results/01_bulk/rnaseq_pheno_GSE268136.rds (step 03)
# Output: results/08_methylation/
#   meth_phenotype_samples.tsv  tumours with methylation, RNA-seq and phenotype score
#   dmp_phenotype.tsv.gz        every CpG: delta beta per SD of the phenotype score, t, p, BH
#   dmp_clusterC.tsv.gz         every CpG: Cluster C vs A+B (CIRCULAR, descriptive only)
#   dmr_phenotype.tsv           differentially methylated regions (DMRcate)
#   dmp_context.tsv             DMPs by genomic context vs all tested CpGs
#   figures/14_dmp_volcano
#
# Phenotype. Ryl's clusters A/B/C were defined from these methylation data, so a
# Cluster C vs A+B comparison of the same data is circular. The primary phenotype is
# therefore the expression-based phenotype score of each tumour's paired RNA-seq
# (step 09): mean z of up- minus down-regulated candidates, with gene sets derived
# from the four other contrasts only. It uses neither the methylation data nor the
# GSE268136 labels. Model per CpG: M-value ~ score (per SD) + platform (450K/EPIC),
# limma; the effect is reported on the beta scale from the same model fitted to beta.
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
if (length(INVASION_EXCLUDE) > 0) stop("Unset INVASION_EXCLUDE.")

load_pkg("limma", bioc = TRUE)
load_pkg("dplyr")
load_pkg("ggplot2")
load_pkg("DMRcate", bioc = TRUE)

beta <- readRDS(file.path(DIR_METH, "meth_beta.rds"))
ann  <- readRDS(file.path(DIR_METH, "meth_probe_annotation.rds"))
smp  <- read.delim(file.path(DIR_METH, "meth_samples.tsv"), stringsAsFactors = FALSE)
smp  <- smp[match(colnames(beta), smp$title), ]
sco  <- read.delim(file.path(DIR_INT, "invasion_scores.tsv"), stringsAsFactors = FALSE)
sco  <- sco[sco$dataset == "GSE268136_ClusterC", ]
clip <- function(b) pmin(pmax(b, METH_BETA_CLIP), 1 - METH_BETA_CLIP)
to_M <- function(b) log2(clip(b) / (1 - clip(b)))

# -----------------------------------------------------------------------------
# 1. Tumours with a phenotype score (paired RNA-seq)
# -----------------------------------------------------------------------------
smp$score <- sco$invasion_score[match(smp$rnaseq_sample, sco$sample_id)]
ps <- smp[!is.na(smp$score), ]
ps$score_sd <- as.numeric(scale(ps$score))
message("Tumours with methylation and phenotype score: ", nrow(ps), " (",
        paste(names(table(ps$platform)), table(ps$platform), collapse = ", "), "); Ryl clusters: ",
        paste(names(table(ps$cluster_ryl_et_al)), table(ps$cluster_ryl_et_al), collapse = ", "))
stopifnot(nrow(ps) >= 20)
save_tsv(ps, file.path(DIR_METH, "meth_phenotype_samples.tsv"))
b_ps <- beta[, ps$title, drop = FALSE]

fit_cpg <- function(b, design, coef) {
  fM <- limma::eBayes(limma::lmFit(to_M(b), design))
  fB <- limma::lmFit(b, design)
  tt <- limma::topTable(fM, coef = coef, number = Inf, sort.by = "none")
  se_b <- fB$stdev.unscaled[, coef] * fB$sigma
  data.frame(cpg = rownames(b), delta_beta = fB$coefficients[, coef], se_delta_beta = se_b,
             t = tt$t, p = tt$P.Value, padj = tt$adj.P.Val, mean_beta = rowMeans(b))
}

# -----------------------------------------------------------------------------
# 2. DMPs: methylation ~ phenotype score + platform
# -----------------------------------------------------------------------------
design <- model.matrix(~ score_sd + platform, data = ps)
dmp <- fit_cpg(b_ps, design, "score_sd") |>
  dplyr::left_join(ann[, c("cpg", "chr", "pos", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Group")], by = "cpg")
dmp$dmp <- dmp$padj < METH_DMP_ADJP & abs(dmp$delta_beta) >= METH_DMP_DB_MIN
dmp$direction <- ifelse(dmp$delta_beta > 0, "hyper_in_aggressive", "hypo_in_aggressive")
message("DMPs (BH < ", METH_DMP_ADJP, ", |delta beta| >= ", METH_DMP_DB_MIN, " per SD): ", sum(dmp$dmp),
        " (", sum(dmp$dmp & dmp$delta_beta > 0), " hyper, ", sum(dmp$dmp & dmp$delta_beta < 0),
        " hypo) of ", nrow(dmp), " CpGs")
data.table::fwrite(dmp, file.path(DIR_METH, "dmp_phenotype.tsv.gz"), sep = "\t")

# Genomic context of DMPs vs all tested CpGs
ctx <- function(d) {
  grp <- ifelse(d$UCSC_RefGene_Group == "" | is.na(d$UCSC_RefGene_Group), "intergenic",
                ifelse(grepl("TSS1500|TSS200|5'UTR|1stExon", d$UCSC_RefGene_Group), "promoter", "gene body/3'UTR"))
  isl <- sub("^[NS]_", "", d$Relation_to_Island)
  data.frame(region = grp, island = isl)
}
cx <- cbind(dmp[, c("dmp", "direction")], ctx(dmp))
context <- dplyr::bind_rows(lapply(c("region", "island"), function(v) {
  tab <- cx |> dplyr::group_by(level = .data[[v]]) |>
    dplyr::summarise(all = dplyr::n(), hyper = sum(dmp & direction == "hyper_in_aggressive"),
                     hypo = sum(dmp & direction == "hypo_in_aggressive"), .groups = "drop")
  tab$variable <- v
  tab$pct_all <- 100 * tab$all / sum(tab$all)
  tab$pct_hyper <- 100 * tab$hyper / max(1, sum(tab$hyper)); tab$pct_hypo <- 100 * tab$hypo / max(1, sum(tab$hypo))
  tab
}))
save_tsv(context, file.path(DIR_METH, "dmp_context.tsv")); print(context)

p_vol <- ggplot2::ggplot(dmp, ggplot2::aes(delta_beta, -log10(p), colour = dmp)) +
  ggplot2::geom_point(size = 0.3, alpha = 0.4) +
  ggplot2::scale_colour_manual(values = c(`TRUE` = COL_INVASIVE, `FALSE` = "grey70"), guide = "none") +
  ggplot2::labs(x = "Delta beta per SD of the phenotype score", y = "-log10 p",
                title = paste0("GSE267015: CpG methylation vs expression-based phenotype score (n = ", nrow(ps), ")")) +
  PLOT_THEME
save_plot(p_vol, file.path(DIR_FIG, "14_dmp_volcano"), width = 7, height = 5.5)

# -----------------------------------------------------------------------------
# 3. Cluster C vs A+B on all used tumours (circular: the clusters come from these data)
# -----------------------------------------------------------------------------
cl <- smp[smp$cluster_ryl_et_al %in% c("A", "B", "C"), ]
cl$clusterC <- as.integer(cl$cluster_ryl_et_al == "C")
if (sum(cl$clusterC) >= 3) {
  dC <- fit_cpg(beta[, cl$title, drop = FALSE], model.matrix(~ clusterC + platform, data = cl), "clusterC")
  data.table::fwrite(dC, file.path(DIR_METH, "dmp_clusterC.tsv.gz"), sep = "\t")
  message("Cluster C vs A+B (", sum(cl$clusterC), " vs ", sum(!cl$clusterC), "; circular): ",
          sum(dC$padj < METH_DMP_ADJP & abs(dC$delta_beta) >= 0.1), " CpGs with BH < 0.05 and |delta beta| >= 0.1")
  message("Correlation of phenotype-score and Cluster C effects across CpGs: r = ",
          round(stats::cor(dmp$t, dC$t[match(dmp$cpg, dC$cpg)]), 2))
}

# -----------------------------------------------------------------------------
# 4. DMRs (DMRcate) for the phenotype score
# -----------------------------------------------------------------------------
dmr <- tryCatch({
  an <- DMRcate::cpg.annotate(datatype = "array", object = to_M(b_ps), what = "M", arraytype = "450K",
                              analysis.type = "differential", design = design, coef = "score_sd",
                              fdr = METH_DMP_ADJP)
  dm <- DMRcate::dmrcate(an, lambda = 1000, C = 2)
  rg <- DMRcate::extractRanges(dm, genome = "hg19")
  as.data.frame(rg)
}, error = function(e) { message("  DMRcate failed: ", conditionMessage(e)); NULL })
if (!is.null(dmr)) {
  dmr$overlapping.genes <- as.character(dmr$overlapping.genes)
  save_tsv(dmr, file.path(DIR_METH, "dmr_phenotype.tsv"))
  message("DMRs: ", nrow(dmr), " (min Fisher FDR < 0.05: ", sum(dmr$min_smoothed_fdr < 0.05, na.rm = TRUE), ")")
}

message("\n[14_methylation_phenotype.R] Complete.")
