# =============================================================================
# Step 04 — Per-contrast differential expression and meta-analysis
#
# Contrasts (primary tumours only; what each label measures is in config.R LABEL_TYPE):
#   GSE97508  — 3 "invasive" vs 3 "non-invasive" (submitter titles; only histology label)
#   GSE29683  — k-means (k = 3) on a proliferation-minus-rod-marker score (this study)
#   GSE59983  — Kooi 2015 expression cluster 1 (n = 26) vs 3 (n = 46); GEO group 2 (n = 4,
#               likely retina contamination per Kooi) excluded
#   GSE110811 — severe (n = 9) vs mild (n = 6) anaplasia (Hudson 2018)
#   GSE268136 — methylation Cluster C vs A+B (Ryl 2024), from step 03
#
# Primary meta-analysis: random-effects (REML tau^2) on Hedges' g from the ordinary
# two-group t (utils.R::meta_re), with tau^2, I^2, HKSJ intervals and 95% prediction
# intervals per gene. BH over genes measured in >= 2 contrasts. Discovery uses the
# Wald z (HKSJ with k <= 5 contrasts cannot reach genome-wide FDR; reported as sensitivity).
# Secondary: signed Stouffer weighted by effective sample size (utils.R::meta_stouffer).
# Datasets in INVASION_EXCLUDE are dropped (sensitivity runs; see config.R).
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
source(here::here("code", "R", "gene_lists.R"))

load_pkg("limma",  bioc = TRUE)
load_pkg("dplyr")
load_pkg("ggplot2")
load_pkg("tibble")

set.seed(GLOBAL_SEED)
MAIN_RUN <- length(INVASION_EXCLUDE) == 0

normalized <- readRDS(file.path(DIR_BULK, "normalized_matrices.rds"))
pheno_list <- readRDS(file.path(DIR_BULK, "pheno_list.rds"))

# -----------------------------------------------------------------------------
# limma per contrast; returns gene, log2FC, t (moderated), t_ord, P.Value, adj.P.Val, g, v
# -----------------------------------------------------------------------------
run_invasion_limma <- function(mat, labels, ds_name) {
  stopifnot(ncol(mat) == length(labels))
  mat    <- mat[complete.cases(mat), , drop = FALSE]
  f      <- factor(labels, levels = c("non-invasive", "invasive"))
  design <- model.matrix(~ f)
  colnames(design) <- c("Intercept", "InvasiveVsNonInvasive")
  fit <- limma::eBayes(limma::lmFit(mat, design), trend = TRUE, robust = TRUE)
  n1 <- sum(labels == "invasive"); n2 <- sum(labels == "non-invasive")
  res <- limma::topTable(fit, coef = "InvasiveVsNonInvasive", number = Inf, sort.by = "none") |>
    tibble::rownames_to_column("gene") |>
    dplyr::rename(log2FC = logFC) |>
    dplyr::mutate(dataset = ds_name)
  t_ord <- fit$coefficients[, 2] / (fit$stdev.unscaled[, 2] * fit$sigma)   # ordinary t for g
  res$t_ord <- unname(t_ord[res$gene]); res$t_ord[!is.finite(res$t_ord)] <- NA
  res <- cbind(res, hedges_g(res$t_ord, n1, n2)) |> dplyr::arrange(P.Value)
  message("  ", ds_name, ": ", n1, " vs ", n2, " samples, ",
          sum(res$adj.P.Val < 0.05), " genes with BH p < 0.05 (",
          round(100 * mean(res$adj.P.Val < 0.05)), "%)")
  attr(res, "n1") <- n1; attr(res, "n2") <- n2
  res
}

# -----------------------------------------------------------------------------
# GSE29683 labels: k-means (k = 3) on a marker score (proliferation minus rod
# photoreceptor markers); highest cluster = "invasive". These labels are this
# study's own and are circular (derived from expression): see the sensitivity runs.
# -----------------------------------------------------------------------------
marker_kmeans <- function(mat) {
  a_genes <- intersect(SUBTYPE_A_MARKERS, rownames(mat))
  c_genes <- intersect(SUBTYPE_C_MARKERS, rownames(mat))
  score_a <- colMeans(mat[a_genes, , drop = FALSE])
  score_c <- colMeans(mat[c_genes, , drop = FALSE])
  disc    <- score_c - score_a
  set.seed(GLOBAL_SEED)
  km  <- kmeans(disc, centers = 3, nstart = 25)
  inv <- as.integer(names(which.max(tapply(disc, km$cluster, mean))))
  data.frame(sample_id = colnames(mat), score_rod = score_a, score_prolif = score_c,
             disc = disc, kmeans_cluster = km$cluster,
             label = ifelse(km$cluster == inv, "invasive", "non-invasive"))
}

# -----------------------------------------------------------------------------
# Labels for all microarray contrasts (primary tumours only, from step 02)
# -----------------------------------------------------------------------------
km29 <- marker_kmeans(normalized$GSE29683[, pheno_list$GSE29683$sample_id])
labels <- list(
  GSE97508  = setNames(label_invasion_gse97508(pheno_list$GSE97508),  pheno_list$GSE97508$sample_id),
  GSE29683  = setNames(km29$label, km29$sample_id),
  GSE59983  = setNames(label_invasion_gse59983(pheno_list$GSE59983),  pheno_list$GSE59983$sample_id),
  GSE110811 = setNames(label_invasion_gse110811(pheno_list$GSE110811), pheno_list$GSE110811$sample_id)
)

inv_results <- list()
for (ds in names(labels)) {
  lab <- labels[[ds]]; lab <- lab[!is.na(lab)]
  message("\n=== ", ds, " (", LABEL_TYPE[[ds]], ") ===")
  inv_results[[ds]] <- run_invasion_limma(normalized[[ds]][, names(lab), drop = FALSE], lab, ds)
}

# GSE268136 Cluster C vs A+B (limma from step 03; group sizes stored as attributes)
rna <- readRDS(file.path(DIR_BULK, "limma_clusterC_vs_AB_GSE268136.rds"))
n1 <- attr(rna, "n1"); n2 <- attr(rna, "n2")
if (is.null(n1) || is.null(rna$t_ord)) stop("Rerun code/03_rnaseq_gse268136.R (group sizes / t_ord missing).")
rna <- cbind(as.data.frame(rna), hedges_g(rna$t_ord, n1, n2)) |>
  dplyr::mutate(dataset = "GSE268136_ClusterC")
attr(rna, "n1") <- n1; attr(rna, "n2") <- n2
inv_results$GSE268136_ClusterC <- rna
message("\n=== GSE268136_ClusterC: ", n1, " vs ", n2, " (", LABEL_TYPE[["GSE268136_ClusterC"]], ")")

# -----------------------------------------------------------------------------
# Per-sample label table (Supporting Information, main run only)
# -----------------------------------------------------------------------------
if (MAIN_RUN) {
  p_rna <- readRDS(file.path(DIR_BULK, "rnaseq_pheno_GSE268136.rds"))
  lab_tbl <- dplyr::bind_rows(
    lapply(names(labels), function(ds) {
      p <- pheno_list[[ds]]
      data.frame(dataset = ds, sample_id = p$sample_id, title = p$title,
                 label = unname(labels[[ds]][p$sample_id]), label_type = LABEL_TYPE[[ds]])
    }),
    data.frame(dataset = "GSE268136_ClusterC", sample_id = rownames(p_rna), title = p_rna$title,
               label = dplyr::case_when(p_rna$is_wt ~ NA_character_,
                                        p_rna$ryl_cluster == "C" ~ "invasive",
                                        p_rna$ryl_cluster %in% c("A", "B") ~ "non-invasive"),
               label_type = LABEL_TYPE[["GSE268136_ClusterC"]],
               ryl_cluster = p_rna$ryl_cluster, rb1_wildtype = p_rna$is_wt)
  ) |>
    dplyr::left_join(km29 |> dplyr::transmute(sample_id, kmeans_score_rod = score_rod,
                                              kmeans_score_prolif = score_prolif,
                                              kmeans_disc = disc, kmeans_cluster), by = "sample_id") |>
    dplyr::mutate(used_in_contrast = !is.na(label))
  save_tsv(lab_tbl, file.path(DIR_INV, "sample_labels.tsv"))
}

# -----------------------------------------------------------------------------
# Meta-analysis
# -----------------------------------------------------------------------------
inv_results <- inv_results[setdiff(names(inv_results), INVASION_EXCLUDE)]
n1s <- vapply(inv_results, function(d) as.numeric(attr(d, "n1")), numeric(1))
n2s <- vapply(inv_results, function(d) as.numeric(attr(d, "n2")), numeric(1))
n_eff <- 4 / (1 / n1s + 1 / n2s)                 # effective sample size of a two-group contrast
attr(inv_results, "n")  <- n_eff                 # read by steps 05 and 09
attr(inv_results, "n1") <- n1s; attr(inv_results, "n2") <- n2s
message("\n=== Meta-analysis across ", length(inv_results), " contrasts ===")
message("  Effective n: ", paste0(names(n_eff), " (", round(n_eff, 1), ")", collapse = ", "))
if (!MAIN_RUN) message("  Excluded: ", paste(INVASION_EXCLUDE, collapse = ", "))

re  <- meta_re(inv_results)
st  <- meta_stouffer(inv_results, n_eff) |>
  dplyr::select(gene, mean_log2FC, stouffer_z = meta_z, stouffer_p = meta_p, stouffer_padj = meta_padj)
meta_inv <- re |>
  dplyr::left_join(st, by = "gene") |>
  dplyr::mutate(label_marker = gene %in% c(SUBTYPE_A_MARKERS, SUBTYPE_C_MARKERS))

invasion_candidates <- meta_inv |>
  dplyr::filter(padj_re < META_ADJP, abs(g_re) >= META_G_MIN, k >= 2)

message("  Genes evaluated: ", nrow(meta_inv), " (", sum(meta_inv$k >= 2), " in >= 2 contrasts, the BH family)")
message("  Candidates (BH p < ", META_ADJP, ", |g| >= ", META_G_MIN, ", >= 2 contrasts): ",
        nrow(invasion_candidates), " (", sum(invasion_candidates$g_re > 0), " up, ",
        sum(invasion_candidates$g_re < 0), " down)")
message("  Label-defining marker genes among candidates (flagged): ",
        paste(invasion_candidates$gene[invasion_candidates$label_marker], collapse = ", "))

# HKSJ sensitivity: how many genes pass with HKSJ inference at several FDR levels
hksj_tbl <- dplyr::bind_rows(lapply(c(2, 4, 5), function(kmin) {
  m <- meta_inv[meta_inv$k >= kmin, ]
  q <- stats::p.adjust(m$p_hksj, "BH"); big <- abs(m$g_re) >= META_G_MIN
  data.frame(min_contrasts = kmin, genes_tested = nrow(m),
             fdr_0.05 = sum(q < 0.05 & big), fdr_0.10 = sum(q < 0.10 & big),
             fdr_0.20 = sum(q < 0.20 & big), nominal_p_0.01 = sum(m$p_hksj < 0.01 & big),
             wald_fdr_0.05 = sum(m$padj_re < META_ADJP & big, na.rm = TRUE))
}))
message("  HKSJ sensitivity (|g| >= ", META_G_MIN, "):"); print(hksj_tbl)

# Per-contrast effects for every gene, wide (Supporting Information + forest plots)
per_contrast <- Reduce(function(a, b) dplyr::full_join(a, b, by = "gene"),
  lapply(names(inv_results), function(ds) {
    d <- inv_results[[ds]][, c("gene", "log2FC", "g", "v", "P.Value", "adj.P.Val")]
    names(d)[-1] <- paste0(names(d)[-1], "__", ds); d
  }))

# -----------------------------------------------------------------------------
# Figures: volcano per contrast
# -----------------------------------------------------------------------------
for (nm in names(inv_results)) {
  d <- inv_results[[nm]]
  p_vol <- plot_volcano(d, title = paste0(sub("_ClusterC", " (Cluster C vs A+B)", nm),
                                          ": ", attr(d, "n1"), " vs ", attr(d, "n2")),
                        log2fc_cut = INV_LOG2FC_CUTOFF, adjp_cut = INV_ADJP_CUTOFF, label_top = 12)
  save_plot(p_vol, file.path(DIR_FIG, paste0("04_volcano_", nm, "_invasion")), width = 8, height = 7)
}

# -----------------------------------------------------------------------------
# GSE97508 quality check (main run): PCA and a residual-retina score
# -----------------------------------------------------------------------------
if (MAIN_RUN) {
  m97 <- normalized$GSE97508[, names(labels$GSE97508)]
  top <- head(order(apply(m97, 1, var), decreasing = TRUE), 2000)
  pc  <- prcomp(t(m97[top, ]), scale. = TRUE)
  retina_genes <- intersect(c(RETINA_PHOTORECEPTOR_MARKERS, RETINA_MULLER_MARKERS), rownames(m97))
  qc97 <- data.frame(sample_id = colnames(m97), label = labels$GSE97508[colnames(m97)],
                     PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                     retina_score = colMeans(zscore_rows(m97[retina_genes, , drop = FALSE])))
  qc97$pct_var_PC1 <- round(100 * summary(pc)$importance[2, 1], 1)
  qc97$pct_genes_bh05 <- round(100 * mean(inv_results$GSE97508$adj.P.Val < 0.05), 1)
  save_tsv(qc97, file.path(DIR_INV, "gse97508_qc.tsv"))
  p_qc <- ggplot2::ggplot(qc97, ggplot2::aes(PC1, PC2, colour = label, size = retina_score)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_manual(values = c(invasive = COL_INVASIVE, "non-invasive" = COL_NONINVASIVE)) +
    ggplot2::labs(title = "GSE97508: PCA of the six tumours",
                  subtitle = "Point size = residual-retina score (photoreceptor + Müller markers)") +
    PLOT_THEME
  save_plot(p_qc, file.path(DIR_FIG, "04_gse97508_qc"), width = 7, height = 5.5)
}

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
save_rds(inv_results, file.path(DIR_INV, "invasion_deg_per_dataset.rds"),
         "per-contrast limma + Hedges' g (attributes: n, n1, n2)")
save_rds(meta_inv,            file.path(DIR_INV, "meta_invasion_all.rds"), "meta-analysis, all genes")
save_rds(invasion_candidates, file.path(DIR_INV, "invasion_candidates.rds"), "candidates")
save_tsv(meta_inv,            file.path(DIR_INV, "meta_invasion_all.tsv"))
save_tsv(invasion_candidates, file.path(DIR_INV, "invasion_candidates.tsv"))
save_tsv(per_contrast,        file.path(DIR_INV, "per_contrast_effects.tsv"))
save_tsv(hksj_tbl,            file.path(DIR_INV, "hksj_sensitivity.tsv"))
if (MAIN_RUN) {
  p29 <- pheno_list$GSE29683; p29$invasion_status <- labels$GSE29683[p29$sample_id]
  p59 <- pheno_list$GSE59983; p59$invasion_status <- labels$GSE59983[p59$sample_id]
  save_rds(p29, file.path(DIR_INV, "pheno_29683_annotated.rds"), "GSE29683 labels")
  save_rds(p59, file.path(DIR_INV, "pheno_59983_annotated.rds"), "GSE59983 labels")
}

message("\n[04_invasion_meta_analysis.R] Complete. ", nrow(invasion_candidates),
        " candidates from ", length(inv_results), " contrasts.")
