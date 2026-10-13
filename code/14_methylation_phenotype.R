# =============================================================================
# Step 14 — CpG methylation associated with the aggressive (de-differentiated) phenotype
# Input : results/08_methylation/ (step 13), results/05_invasion_immune/invasion_scores.tsv
#         (step 09), results/01_bulk/rnaseq_pheno_GSE268136.rds (step 03)
#         results/04_immune/estimate_purity.tsv (step 08), results/01_bulk/rnaseq_log2tpm_GSE268136.rds
# Output: results/08_methylation/
#   meth_phenotype_samples.tsv  tumours with methylation, RNA-seq and phenotype score
#   dmp_phenotype.tsv.gz        every CpG: delta beta per SD of the phenotype score, t, p, BH
#   dmp_clusterC.tsv.gz         every CpG: Cluster C vs A+B (CIRCULAR, descriptive only)
#   dmr_phenotype.tsv           differentially methylated regions (DMRcate)
#   dmp_context.tsv             DMPs by genomic context vs all tested CpGs
#   dmp_sensitivity.tsv         DMP counts under platform, composition and cluster sensitivity models
#   dmp_sensitivity_cpg.tsv.gz  every main DMP under each sensitivity model
#   dmr_phenotype_composition.tsv  DMRs with composition covariates
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
source(here::here("code", "R", "gene_lists.R"))
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

# -----------------------------------------------------------------------------
# 5. Sensitivity analyses (revision): platform, cell composition, Ryl cluster
#    platform_450K / platform_EPIC  each array alone
#    composition    + ESTIMATE StromalScore and purity (immune + stromal) and a Muller-glia
#                   (residual-retina) score, all from the paired RNA-seq. Photoreceptor genes are
#                   left out of the retina score because tumour cells express cone genes.
#    within_cluster_A  Ryl cluster A tumours only (the only cluster large enough)
#    cluster_adjusted  + Ryl cluster
#    A main DMP is "reproduced" when the sensitivity model gives the same sign and p < 0.05.
# -----------------------------------------------------------------------------
est <- read.delim(file.path(DIR_IMM, "estimate_purity.tsv"), stringsAsFactors = FALSE)
est <- est[est$dataset == "GSE268136", ]
rna <- readRDS(file.path(DIR_BULK, "rnaseq_log2tpm_GSE268136.rds"))
mg  <- intersect(RETINA_MULLER_MARKERS, rownames(rna))
ps$stroma  <- est$stroma_score[match(ps$rnaseq_sample, est$sample_id)]
ps$purity  <- est$purity[match(ps$rnaseq_sample, est$sample_id)]
ps$muller  <- colMeans(zscore_rows(rna[mg, ps$rnaseq_sample, drop = FALSE]))
stopifnot(!anyNA(ps[, c("stroma", "purity", "muller")]), length(mg) >= 3)
message("Composition covariates vs phenotype score (Spearman): StromalScore ",
        round(stats::cor(ps$stroma, ps$score, method = "spearman"), 2), ", purity ",
        round(stats::cor(ps$purity, ps$score, method = "spearman"), 2), ", Muller-glia score ",
        round(stats::cor(ps$muller, ps$score, method = "spearman"), 2))
sens_models <- list(
  platform_450K    = list(rows = ps$platform == "450K", f = ~ score_sd),
  platform_EPIC    = list(rows = ps$platform == "EPIC", f = ~ score_sd),
  composition      = list(rows = rep(TRUE, nrow(ps)), f = ~ score_sd + platform + scale(stroma) + scale(purity) + scale(muller)),
  within_cluster_A = list(rows = ps$cluster_ryl_et_al == "A", f = ~ score_sd + platform),
  cluster_adjusted = list(rows = rep(TRUE, nrow(ps)), f = ~ score_sd + platform + cluster_ryl_et_al))
main_dmp <- dmp$cpg[dmp$dmp]
sens_cpg <- data.frame(cpg = main_dmp, main_delta_beta = dmp$delta_beta[dmp$dmp])
sens <- dplyr::bind_rows(lapply(names(sens_models), function(nm) {
  m <- sens_models[[nm]]; d <- ps[m$rows, ]
  d$score_sd <- as.numeric(scale(d$score))            # SD units within the subset
  des <- model.matrix(m$f, data = d)
  r <- fit_cpg(b_ps[, d$title, drop = FALSE], des, "score_sd")
  is_dmp <- r$padj < METH_DMP_ADJP & abs(r$delta_beta) >= METH_DMP_DB_MIN
  i <- match(main_dmp, r$cpg)
  rep_ok <- sign(r$delta_beta[i]) == sign(sens_cpg$main_delta_beta) & r$p[i] < 0.05
  sens_cpg[[paste0(nm, "_delta_beta")]] <<- r$delta_beta[i]
  sens_cpg[[paste0(nm, "_padj")]] <<- r$padj[i]
  data.frame(model = nm, n_tumours = nrow(d), n_dmp = sum(is_dmp), n_dmp_hyper = sum(is_dmp & r$delta_beta > 0),
             n_dmp_hypo = sum(is_dmp & r$delta_beta < 0), main_dmp_reproduced = sum(rep_ok),
             pct_main_dmp_reproduced = 100 * mean(rep_ok), main_dmp_still_dmp = sum(is_dmp[i]),
             r_t_with_main = stats::cor(r$t, dmp$t[match(r$cpg, dmp$cpg)]))
}))
sens <- dplyr::bind_rows(data.frame(model = "main", n_tumours = nrow(ps), n_dmp = sum(dmp$dmp),
                                    n_dmp_hyper = sum(dmp$dmp & dmp$delta_beta > 0), n_dmp_hypo = sum(dmp$dmp & dmp$delta_beta < 0),
                                    main_dmp_reproduced = sum(dmp$dmp), pct_main_dmp_reproduced = 100,
                                    main_dmp_still_dmp = sum(dmp$dmp), r_t_with_main = 1), sens)
save_tsv(sens, file.path(DIR_METH, "dmp_sensitivity.tsv")); print(sens)
data.table::fwrite(sens_cpg, file.path(DIR_METH, "dmp_sensitivity_cpg.tsv.gz"), sep = "\t")
save_tsv(ps, file.path(DIR_METH, "meth_phenotype_samples.tsv"))   # now with the composition covariates

# DMRs under the composition-adjusted model
dmr_c <- tryCatch({
  des_c <- model.matrix(sens_models$composition$f, data = ps)
  an <- DMRcate::cpg.annotate(datatype = "array", object = to_M(b_ps), what = "M", arraytype = "450K",
                              analysis.type = "differential", design = des_c, coef = "score_sd", fdr = METH_DMP_ADJP)
  as.data.frame(DMRcate::extractRanges(DMRcate::dmrcate(an, lambda = 1000, C = 2), genome = "hg19"))
}, error = function(e) { message("  DMRcate (composition-adjusted) failed: ", conditionMessage(e)); NULL })
if (!is.null(dmr_c)) {
  dmr_c$overlapping.genes <- as.character(dmr_c$overlapping.genes)
  save_tsv(dmr_c, file.path(DIR_METH, "dmr_phenotype_composition.tsv"))
  message("DMRs, composition-adjusted: ", nrow(dmr_c), " (", sum(dmr_c$meandiff > 0), " hyper)")
}

message("\n[14_methylation_phenotype.R] Complete.")
