# =============================================================================
# Step 09 — Phenotype score and immune co-variation (labelled primary tumours only)
# Input : results/01_bulk/ (matrices), results/02_invasion/ (per-contrast results,
#         candidates, sample_labels.tsv), results/04_immune/ (MCP-counter, ESTIMATE)
# Output: results/05_invasion_immune/
#   invasion_scores.tsv          per-sample score (leave-one-dataset-out gene sets)
#   score_validation.tsv         score by label per contrast (Wilcoxon, n, medians)
#   invasion_immune_corr_celltype.tsv  score vs MCP-counter cell types (+ StromalScore-partial,
#                                      sensitivity; not for GSE268136, see below)
#   gene_invasion_immune_corr.tsv      per gene: rho with the score (leave-one-gene-out)
#                                      and with each MCP-counter variable, meta-analysed
#   covarying_genes.tsv                genes co-varying with the score and an immune variable
#   figures/09_invasion_score_validation
# Score: mean z(up genes) - mean z(down genes), genes z-scored within each cohort's
# labelled primary tumours; up/down sets are candidates from a meta-analysis that
# excludes the scored contrast. A gene's own correlation with the score uses a score
# computed without that gene.
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
if (length(INVASION_EXCLUDE) > 0) stop("Unset INVASION_EXCLUDE: step 09 uses the main results.")

load_pkg("dplyr")
load_pkg("tidyr")
load_pkg("ggplot2")

set.seed(GLOBAL_SEED)

normalized  <- readRDS(file.path(DIR_BULK, "normalized_matrices.rds"))
rnaseq_mat  <- readRDS(file.path(DIR_BULK, "rnaseq_log2tpm_GSE268136.rds"))
inv_results <- readRDS(file.path(DIR_INV, "invasion_deg_per_dataset.rds"))
cand        <- readRDS(file.path(DIR_INV, "invasion_candidates.rds"))
lab_tbl     <- read.delim(file.path(DIR_INV, "sample_labels.tsv"), stringsAsFactors = FALSE)
deconv_long <- readRDS(file.path(DIR_IMM, "deconv_long_all_datasets.rds"))
purity      <- read.delim(file.path(DIR_IMM, "estimate_purity.tsv"), stringsAsFactors = FALSE)

CONTRASTS <- names(inv_results)
mat_of <- function(ds) if (ds == "GSE268136_ClusterC") rnaseq_mat else normalized[[ds]]
cohort_of <- function(ds) sub("_ClusterC", "", ds)

stouffer_signed_p <- function(ps, rhos) {
  ok <- !is.na(ps) & !is.na(rhos)
  if (!any(ok)) return(NA_real_)
  z <- sign(rhos[ok]) * stats::qnorm(pmax(ps[ok], 1e-300) / 2, lower.tail = FALSE)
  as.numeric(2 * stats::pnorm(-abs(sum(z) / sqrt(length(z)))))
}
stopifnot(stouffer_signed_p(c(1e-4, 1e-4), c(0.5, -0.5)) > 0.99,
          stouffer_signed_p(c(1e-4, 1e-4), c(0.5,  0.5)) < 1e-6)

# -----------------------------------------------------------------------------
# 1. Scores on labelled primary tumours, leave-one-dataset-out gene sets
# -----------------------------------------------------------------------------
cohorts <- lapply(setNames(CONTRASTS, CONTRASTS), function(ds) {
  l <- lab_tbl[lab_tbl$dataset == ds & !is.na(lab_tbl$label), ]
  l <- l[l$sample_id %in% colnames(mat_of(ds)), ]
  m <- mat_of(ds)[, l$sample_id, drop = FALSE]
  Z <- zscore_rows(m); Z <- Z[stats::complete.cases(Z), , drop = FALSE]
  mm <- meta_re(inv_results[setdiff(CONTRASTS, ds)])
  mm <- mm[mm$padj_re < META_ADJP & abs(mm$g_re) >= META_G_MIN & mm$k >= 2, ]
  up <- intersect(mm$gene[mm$g_re > 0], rownames(Z)); down <- intersect(mm$gene[mm$g_re < 0], rownames(Z))
  message("  ", ds, ": ", ncol(Z), " labelled tumours; LODO sets up ", length(up), ", down ", length(down))
  list(Z = Z, label = setNames(l$label, l$sample_id), up = up, down = down,
       score = score_updown(Z, up, down))
})

scores <- dplyr::bind_rows(lapply(CONTRASTS, function(ds) data.frame(
  dataset = ds, sample_id = names(cohorts[[ds]]$score),
  invasion_score = unname(cohorts[[ds]]$score),
  invasion_label = unname(cohorts[[ds]]$label[names(cohorts[[ds]]$score)]))))
save_tsv(scores, file.path(DIR_INT, "invasion_scores.tsv"))

val <- scores |> dplyr::group_by(dataset) |>
  dplyr::summarise(n_inv = sum(invasion_label == "invasive"), n_noninv = sum(invasion_label == "non-invasive"),
                   median_inv = median(invasion_score[invasion_label == "invasive"]),
                   median_noninv = median(invasion_score[invasion_label == "non-invasive"]),
                   wilcox_p = stats::wilcox.test(invasion_score[invasion_label == "invasive"],
                                                 invasion_score[invasion_label == "non-invasive"])$p.value,
                   .groups = "drop")
save_tsv(val, file.path(DIR_INT, "score_validation.tsv"))
print(val)

lab_facet <- setNames(paste0(sub("_ClusterC", " (C vs A+B)", val$dataset), "\n", val$n_inv, " vs ", val$n_noninv,
                             ", p = ", signif(val$wilcox_p, 2)), val$dataset)
p_val <- ggplot2::ggplot(scores, ggplot2::aes(invasion_label, invasion_score, fill = invasion_label)) +
  ggplot2::geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.6) +
  ggplot2::geom_jitter(width = 0.12, size = 1, alpha = 0.7) +
  ggplot2::facet_wrap(~ dataset, nrow = 1, scales = "free_y", labeller = ggplot2::as_labeller(lab_facet)) +
  ggplot2::scale_fill_manual(values = c(invasive = COL_INVASIVE, "non-invasive" = COL_NONINVASIVE)) +
  ggplot2::scale_x_discrete(labels = c(invasive = "aggressive", "non-invasive" = "reference")) +
  ggplot2::labs(x = NULL, y = "Score (mean z up - mean z down)",
                title = "Phenotype score from gene sets derived without the scored cohort") +
  PLOT_THEME + ggplot2::theme(legend.position = "none")
save_plot(p_val, file.path(DIR_FIG, "09_invasion_score_validation"), width = 14, height = 4.5)

# -----------------------------------------------------------------------------
# 2. Score vs MCP-counter cell types (Spearman; StromalScore-partial Spearman as a
#    sensitivity analysis). ESTIMATE purity is not used as the covariate because its
#    ImmuneScore shares genes with MCP-counter. No partials for GSE268136 (4 vs 43;
#    ESTIMATE was built for microarrays and its purity there spans 0.61-1.00).
# -----------------------------------------------------------------------------
NO_PARTIAL <- "GSE268136_ClusterC"
CORR_SETS <- intersect(c("GSE59983", "GSE29683", "GSE110811", "GSE268136_ClusterC"), CONTRASTS)
mcp <- deconv_long |> dplyr::filter(method == "mcp_counter")
purity_of <- function(ds, smp) {   # StromalScore covariate (NA = no partial correlation)
  if (ds %in% NO_PARTIAL) return(rep(NA_real_, length(smp)))
  p <- purity[purity$dataset == cohort_of(ds), ]
  p$stroma_score[match(smp, p$sample_id)]
}
corr_ct <- dplyr::bind_rows(lapply(CORR_SETS, function(ds) {
  s  <- cohorts[[ds]]$score
  d  <- mcp |> dplyr::filter(dataset == cohort_of(ds), sample_id %in% names(s))
  pu <- purity_of(ds, names(s))
  dplyr::bind_rows(lapply(split(d, d$cell_type), function(x) {
    v  <- x$score[match(names(s), x$sample_id)]
    ct <- suppressWarnings(stats::cor.test(s, v, method = "spearman", exact = FALSE))
    pp <- partial_spearman(s, v, pu)
    data.frame(dataset = ds, cell_type = x$cell_type[1], n = sum(!is.na(v)),
               rho = unname(ct$estimate), p = ct$p.value,
               rho_partial_stromal = unname(pp["rho"]), p_partial_stromal = unname(pp["p"]))
  }))
})) |>
  dplyr::group_by(dataset) |>
  dplyr::mutate(padj = p.adjust(p, "BH"), padj_partial_stromal = p.adjust(p_partial_stromal, "BH")) |>
  dplyr::ungroup()
save_tsv(corr_ct, file.path(DIR_INT, "invasion_immune_corr_celltype.tsv"))

# -----------------------------------------------------------------------------
# 3. Gene-level correlations: candidate genes vs score (leave-one-gene-out) and
#    vs each MCP-counter variable; macrophage also StromalScore-partial
# -----------------------------------------------------------------------------
IMMUNE_VARS <- sort(unique(mcp$cell_type))
gene_corr <- dplyr::bind_rows(lapply(CORR_SETS, function(ds) {
  co <- cohorts[[ds]]; Z <- co$Z; smp <- colnames(Z)
  mu <- colMeans(Z[co$up, , drop = FALSE]); md <- colMeans(Z[co$down, , drop = FALSE])
  nu <- length(co$up); nd <- length(co$down)
  iv <- lapply(setNames(IMMUNE_VARS, IMMUNE_VARS), function(v) {
    x <- mcp[mcp$dataset == cohort_of(ds) & mcp$cell_type == v, ]; x$score[match(smp, x$sample_id)] })
  pu <- purity_of(ds, smp)
  genes <- intersect(cand$gene, rownames(Z))
  dplyr::bind_rows(lapply(genes, function(g) {
    e <- Z[g, ]
    s <- if (g %in% co$up) (mu * nu - e) / (nu - 1) - md else if (g %in% co$down) mu - (md * nd - e) / (nd - 1) else mu - md
    r <- suppressWarnings(stats::cor.test(e, s, method = "spearman", exact = FALSE))
    out <- data.frame(gene = g, dataset = ds, rho_invasion = unname(r$estimate), p_invasion = r$p.value)
    for (v in IMMUNE_VARS) {
      ri <- suppressWarnings(stats::cor.test(e, iv[[v]], method = "spearman", exact = FALSE))
      out[[paste0("rho__", v)]] <- unname(ri$estimate); out[[paste0("p__", v)]] <- ri$p.value
    }
    pp <- partial_spearman(e, iv[["Macrophage/Monocyte"]], pu)
    out$rho_partial_macrophage <- unname(pp["rho"])
    out
  }))
}))

meta_g <- gene_corr |> dplyr::group_by(gene) |>
  dplyr::summarise(n_datasets = dplyr::n(),
                   mean_rho_invasion = mean(rho_invasion, na.rm = TRUE),
                   meta_p_invasion = stouffer_signed_p(p_invasion, rho_invasion),
                   mean_rho_partial_macrophage = mean(rho_partial_macrophage, na.rm = TRUE),
                   dplyr::across(dplyr::starts_with("rho__"), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
                   .groups = "drop") |>
  dplyr::mutate(meta_padj_invasion = p.adjust(meta_p_invasion, "BH"))
for (v in IMMUNE_VARS) {
  pv <- gene_corr |> dplyr::group_by(gene) |>
    dplyr::summarise(p = stouffer_signed_p(.data[[paste0("p__", v)]], .data[[paste0("rho__", v)]]), .groups = "drop")
  meta_g[[paste0("meta_padj__", v)]] <- p.adjust(pv$p[match(meta_g$gene, pv$gene)], "BH")
}
names(meta_g) <- gsub("[^A-Za-z0-9_]", "_", names(meta_g))
save_tsv(meta_g, file.path(DIR_INT, "gene_invasion_immune_corr.tsv"))

# Co-varying genes: score criterion AND, for the same immune variable, |mean rho| >= 0.2 and padj < 0.05
vn <- gsub("[^A-Za-z0-9_]", "_", IMMUNE_VARS)
hit <- sapply(vn, function(v) abs(meta_g[[paste0("mean_rho__", v)]]) >= 0.2 & meta_g[[paste0("meta_padj__", v)]] < 0.05)
keep <- abs(meta_g$mean_rho_invasion) >= 0.2 & meta_g$meta_padj_invasion < 0.05 & rowSums(hit, na.rm = TRUE) > 0
cov_genes <- meta_g[keep, ]
cov_genes$immune_variables <- apply(hit[keep, , drop = FALSE], 1, function(h) paste(vn[which(h)], collapse = ";"))
save_tsv(cov_genes, file.path(DIR_INT, "covarying_genes.tsv"))

message("\n[09_invasion_immune.R] Complete. Co-varying genes: ", nrow(cov_genes), " of ",
        nrow(meta_g), " tested")
