# =============================================================================
# Step 11 — Priority ranking for FFPE testing
#
# Filters:
#   F1 random-effects BH p < META_ADJP (genes in >= 2 contrasts) and |pooled g| >= META_G_MIN
#   F2 measured in >= PRI_MIN_CONTRASTS contrasts
#   F3 >= PRI_MIN_AGREE contrasts agree in sign with the pooled effect
#   F4 same sign and nominal p < 0.05 for the pooled effect without GSE29683 and without
#      GSE29683 + GSE59983. History: F4 was first set at BH level; after the first full
#      run it was relaxed to nominal p. Both versions are reported (scenario "F4_BH").
#   F5 not a gene used to define the GSE29683 labels
#   F6 same direction in GSE97508, the only histology-labelled contrast (scenario without it)
#   F7 detected in >= PRI_SC_MIN_PCT % of malignant cells in >= PRI_SC_MIN_TUMOURS
#      single-cell tumours (added in revision 2: excludes genes that bulk tissue carries
#      from non-tumour cells; scenarios without it and with all tumour-lineage cells
#      instead of malignant cells, because the malignant calls are uncertain)
# Ranking: lower bound of the 95% CI of |pooled g| (REML), within each direction.
#   Robust rank aggregation (RRA) of per-contrast ranks is a sensitivity scenario.
# Shortlist: top PRI_N_PER_DIR up + top PRI_N_PER_DIR down among genes passing F1-F7.
# Stability: N_BOOT bootstrap resamples of tumours within each contrast (stratified by
#   label); per resample the per-contrast g, both meta-analyses, F1-F7 and the ranking are
#   recomputed. Core genes = in the top PRI_N_PER_DIR in >= PRI_CORE_PROB of resamples.
# Reported, not filtered: I^2, tau^2, HKSJ CI, prediction interval, subtype-adjusted
#   effects, rho with the phenotype score, StromalScore-partial macrophage rho, LASSO
#   selection frequency, GSE249995 description.
# Output: results/07_priority/; figures/11_*
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
source(here::here("code", "R", "gene_lists.R"))
if (length(INVASION_EXCLUDE) > 0) stop("Unset INVASION_EXCLUDE.")

load_pkg("dplyr")
load_pkg("tidyr")
load_pkg("ggplot2")

inv_results <- readRDS(file.path(DIR_INV, "invasion_deg_per_dataset.rds"))
meta <- readRDS(file.path(DIR_INV, "meta_invasion_all.rds"))
sens <- lapply(setNames(SENS_EXCLUDES, SENS_EXCLUDES), function(x) {
  d <- paste0(DIR_INV, "_excl_", gsub(",", "_", x))
  f <- file.path(d, "meta_invasion_all.rds")
  if (!file.exists(f)) stop("Missing sensitivity run: INVASION_EXCLUDE=", x, " Rscript code/04_invasion_meta_analysis.R")
  readRDS(f)
})
pb <- readRDS(file.path(DIR_SC, "pseudobulk_malignant.rds"))
CONTRASTS <- names(inv_results)
LABEL_MARKERS <- c(SUBTYPE_A_MARKERS, SUBTYPE_C_MARKERS)
Z975 <- stats::qnorm(0.975)

# -----------------------------------------------------------------------------
# Per-contrast effects (forest plots, F6) and RRA ranks (sensitivity scenario)
# -----------------------------------------------------------------------------
G <- V <- matrix(NA_real_, nrow(meta), length(CONTRASTS), dimnames = list(meta$gene, CONTRASTS))
for (ds in CONTRASTS) {
  d <- inv_results[[ds]]; i <- match(d$gene, meta$gene)
  G[i[!is.na(i)], ds] <- d$g[!is.na(i)]; V[i[!is.na(i)], ds] <- d$v[!is.na(i)]
}
norm_rank <- function(x, decreasing) {
  r <- rank(if (decreasing) -x else x, na.last = "keep", ties.method = "average")
  r / sum(!is.na(x))
}
R_up   <- apply(G, 2, norm_rank, decreasing = TRUE)
R_down <- apply(G, 2, norm_rank, decreasing = FALSE)

# F7: genes detected in malignant cells of enough single-cell tumours
n_tum <- rowSums(pb$pct >= PRI_SC_MIN_PCT)
SC_GENES <- names(n_tum)[n_tum >= PRI_SC_MIN_TUMOURS]
n_tum_lin <- rowSums(pb$lineage_pct >= PRI_SC_MIN_PCT)          # same rule in all tumour-lineage cells
SC_GENES_LIN <- names(n_tum_lin)[n_tum_lin >= PRI_SC_MIN_TUMOURS]
message("F7: ", length(SC_GENES), " genes detected in >= ", PRI_SC_MIN_PCT, "% of malignant cells in >= ",
        PRI_SC_MIN_TUMOURS, " of ", ncol(pb$pct), " tumours")

# -----------------------------------------------------------------------------
# Filters and ranking (written for any meta-analysis table, so the bootstrap reuses them)
# -----------------------------------------------------------------------------
sens_ok <- function(m, genes, g_main, mode = "nominal") {
  i <- match(genes, m$gene)
  sig <- if (mode == "nominal") m$p_re[i] < 0.05 else m$padj_re[i] < META_ADJP & abs(m$g_re[i]) >= META_G_MIN
  (!is.na(i) & sig & sign(m$g_re[i]) == sign(g_main)) %in% TRUE
}
apply_filters <- function(m, g97, sl, f4 = "nominal", use_f6 = TRUE, use_f7 = TRUE, sc_genes = SC_GENES) {
  # g97: named vector of GSE97508 g; sl: list of sensitivity meta-analyses
  m$F1 <- (m$padj_re < META_ADJP & abs(m$g_re) >= META_G_MIN) %in% TRUE
  m$F2 <- m$k >= PRI_MIN_CONTRASTS
  m$F3 <- m$n_agree >= PRI_MIN_AGREE
  m$F4 <- Reduce(`&`, lapply(sl, sens_ok, genes = m$gene, g_main = m$g_re, mode = f4))
  m$F5 <- !m$gene %in% LABEL_MARKERS
  x97  <- g97[m$gene]
  m$F6 <- if (use_f6) is.na(x97) | sign(x97) == sign(m$g_re) else TRUE
  m$F7 <- if (use_f7) m$gene %in% sc_genes else TRUE
  m$pass <- m$F1 & m$F2 & m$F3 & m$F4 & m$F5 & m$F6 & m$F7
  m
}
lower_bound <- function(m) abs(m$g_re) - Z975 * m$se_re
shortlist_of <- function(m, rank = "lower", n = PRI_N_PER_DIR) {
  p <- m[m$pass, ]
  score <- if (rank == "lower") -lower_bound(p) else
    ifelse(p$g_re > 0, rra_rho(R_up[p$gene, , drop = FALSE]), rra_rho(R_down[p$gene, , drop = FALSE]))
  p <- p[order(score, p$p_re), ]
  c(head(p$gene[p$g_re > 0], n), head(p$gene[p$g_re < 0], n))
}

g97_main <- setNames(G[, "GSE97508"], rownames(G))
main <- apply_filters(meta, g97_main, sens)
main$lower_bound <- lower_bound(main)
main$rra <- NA_real_
up_i <- which(main$F1 & main$g_re > 0); dn_i <- which(main$F1 & main$g_re < 0)
main$rra[up_i] <- rra_rho(R_up[main$gene[up_i], , drop = FALSE])
main$rra[dn_i] <- rra_rho(R_down[main$gene[dn_i], , drop = FALSE])
SHORT <- shortlist_of(main)
main$shortlist <- main$gene %in% SHORT
cum <- Reduce(`&`, main[paste0("F", 1:7)], accumulate = TRUE)
message("Filter counts: ", paste0("F1-F", 1:7, " ", vapply(cum, sum, numeric(1)), collapse = " | "))
message("Shortlist (lower CI bound): ", paste(SHORT, collapse = ", "))
stopifnot(length(SHORT) >= 2, all(main$pass[main$shortlist]))

# -----------------------------------------------------------------------------
# Sensitivity scenarios for the choices that change the shortlist
# -----------------------------------------------------------------------------
SCEN <- list(main = list(), F4_BH = list(f4 = "bh"), rank_RRA = list(rank = "rra"),
             top10 = list(n = 10), without_F6 = list(use_f6 = FALSE), without_F7 = list(use_f7 = FALSE),
             F7_tumour_lineage = list(sc_genes = SC_GENES_LIN))
scen_tbl <- dplyr::bind_rows(lapply(names(SCEN), function(nm) {
  a <- SCEN[[nm]]
  d <- do.call(apply_filters, c(list(meta, g97_main, sens), a[intersect(names(a), c("f4", "use_f6", "use_f7", "sc_genes"))]))
  s <- shortlist_of(d, rank = if (is.null(a$rank)) "lower" else a$rank, n = if (is.null(a$n)) PRI_N_PER_DIR else a$n)
  data.frame(scenario = nm, n_pass = sum(d$pass), n_pass_up = sum(d$pass & d$g_re > 0),
             n_pass_down = sum(d$pass & d$g_re < 0), main_shortlist_retained = length(intersect(s, SHORT)),
             shortlist = paste(s, collapse = ","))
}))
save_tsv(scen_tbl, file.path(DIR_PRI, "shortlist_scenarios.tsv"))
print(scen_tbl[, c("scenario", "n_pass", "main_shortlist_retained", "shortlist")])

# -----------------------------------------------------------------------------
# Bootstrap: tumours resampled within each contrast and label
# -----------------------------------------------------------------------------
normalized <- readRDS(file.path(DIR_BULK, "normalized_matrices.rds"))
rnaseq_mat <- readRDS(file.path(DIR_BULK, "rnaseq_log2tpm_GSE268136.rds"))
lab_tbl    <- read.delim(file.path(DIR_INV, "sample_labels.tsv"), stringsAsFactors = FALSE)
mat_of <- function(ds) if (ds == "GSE268136_ClusterC") rnaseq_mat else normalized[[ds]]
cohorts <- lapply(setNames(CONTRASTS, CONTRASTS), function(ds) {
  l <- lab_tbl[lab_tbl$dataset == ds & !is.na(lab_tbl$label), ]
  list(X = mat_of(ds)[inv_results[[ds]]$gene, l$sample_id, drop = FALSE], inv = l$label == "invasive")
})
# Hedges' g from the ordinary two-group t (same as step 04)
g_of <- function(X, inv) {
  n1 <- sum(inv); n2 <- sum(!inv)
  m1 <- rowMeans(X[, inv, drop = FALSE]); m2 <- rowMeans(X[, !inv, drop = FALSE])
  ss <- rowSums((X[, inv, drop = FALSE] - m1)^2) + rowSums((X[, !inv, drop = FALSE] - m2)^2)
  t  <- (m1 - m2) / (sqrt(ss / (n1 + n2 - 2)) * sqrt(1 / n1 + 1 / n2))
  t[!is.finite(t)] <- NA
  data.frame(gene = rownames(X), hedges_g(t, n1, n2))
}
# Self-check: on the original tumours g_of reproduces step 04's g for every contrast
# (genes whose group means are equal to machine precision are skipped: their t is
# rounding noise in either implementation)
for (ds in CONTRASTS) {
  chk <- g_of(cohorts[[ds]]$X, cohorts[[ds]]$inv)
  real <- abs(inv_results[[ds]]$log2FC) > 1e-8
  dif <- max(abs(chk$g - inv_results[[ds]]$g)[real], na.rm = TRUE)
  if (dif > 1e-6) stop("Bootstrap g does not reproduce step 04 for ", ds, " (max diff ", signif(dif, 3), ")")
}
excl_sets <- lapply(SENS_EXCLUDES, function(x) strsplit(x, ",")[[1]])
message("\n=== Bootstrap (", N_BOOT, " resamples) ===")
boot <- lapply(seq_len(N_BOOT), function(b) {
  set.seed(GLOBAL_SEED + b)
  res <- lapply(cohorts, function(co) {
    idx <- c(sample(which(co$inv), replace = TRUE), sample(which(!co$inv), replace = TRUE))
    g_of(co$X[, idx, drop = FALSE], co$inv[idx])
  })
  m  <- meta_re(res)
  sl <- lapply(excl_sets, function(x) meta_re(res[setdiff(names(res), x)]))
  d  <- apply_filters(m, setNames(res$GSE97508$g, res$GSE97508$gene), sl)
  if (b %% 20 == 0) message("  ", b, " / ", N_BOOT)
  list(top5 = shortlist_of(d), top10 = shortlist_of(d, n = 10), pass = d$gene[d$pass])
})
tab <- function(x) table(unlist(lapply(boot, `[[`, x)))
t5 <- tab("top5"); t10 <- tab("top10"); tp <- tab("pass")
bp <- data.frame(gene = names(tp), p_pass = as.numeric(tp) / N_BOOT)
bp$p_top5  <- as.numeric(t5[bp$gene]) / N_BOOT;  bp$p_top5[is.na(bp$p_top5)] <- 0
bp$p_top10 <- as.numeric(t10[bp$gene]) / N_BOOT; bp$p_top10[is.na(bp$p_top10)] <- 0
bp$direction <- ifelse(meta$g_re[match(bp$gene, meta$gene)] > 0, "up", "down")
bp$main_shortlist <- bp$gene %in% SHORT
bp <- bp[order(-bp$p_top5, -bp$p_top10, -bp$p_pass), ]
CORE <- bp$gene[bp$p_top5 >= PRI_CORE_PROB]
bp$core <- bp$gene %in% CORE
save_tsv(bp, file.path(DIR_PRI, "bootstrap_shortlist_probability.tsv"))
message("Core genes (top ", PRI_N_PER_DIR, " in >= ", 100 * PRI_CORE_PROB, "% of resamples): ",
        paste(CORE, collapse = ", "))
print(head(bp, 25))

# -----------------------------------------------------------------------------
# Annotate candidates with the reported (non-filter) evidence
# -----------------------------------------------------------------------------
cand <- main[main$F1, ]
gc  <- read.delim(file.path(DIR_INT, "gene_invasion_immune_corr.tsv"), stringsAsFactors = FALSE)
adj <- read.delim(file.path(DIR_SUB, "subtype_adjusted_effects.tsv"), stringsAsFactors = FALSE)
stb <- read.delim(file.path(DIR_INV, "lasso_stability.tsv"), stringsAsFactors = FALSE)
scg <- read.delim(file.path(DIR_INV, "scrna_gse249995_candidates.tsv"), stringsAsFactors = FALSE)
adj_w <- tidyr::pivot_wider(adj[, c("gene", "signature", "beta_adjusted", "beta_adjusted_lo", "beta_adjusted_hi",
                                    "padj_adjusted")],
                            names_from = signature, values_from = c(beta_adjusted, beta_adjusted_lo,
                                                                    beta_adjusted_hi, padj_adjusted))
cand <- cand |>
  dplyr::left_join(gc[, intersect(names(gc), c("gene", "mean_rho_invasion", "meta_padj_invasion",
                                               "mean_rho__Macrophage_Monocyte", "mean_rho_partial_macrophage"))], by = "gene") |>
  dplyr::left_join(adj_w, by = "gene") |>
  dplyr::left_join(stb[, c("gene", "selection_freq")], by = "gene") |>
  dplyr::left_join(scg[, c("gene", "concordant", "concordant_lineage")] |>
                     dplyr::rename(scrna_gse249995_concordant = concordant,
                                   scrna_gse249995_concordant_lineage = concordant_lineage), by = "gene") |>
  dplyr::left_join(bp[, c("gene", "p_pass", "p_top5", "p_top10", "core")], by = "gene")
cand$sc_tumours_detected <- unname(n_tum[cand$gene])
gw <- as.data.frame(G[cand$gene, , drop = FALSE]); names(gw) <- paste0("g__", names(gw))
cand <- cbind(cand, gw) |> dplyr::arrange(dplyr::desc(shortlist), dplyr::desc(pass), dplyr::desc(lower_bound))
save_tsv(cand, file.path(DIR_PRI, "priority_all_candidates.tsv"))
save_tsv(cand[cand$shortlist | cand$core %in% TRUE, ], file.path(DIR_PRI, "priority_shortlist.tsv"))

# -----------------------------------------------------------------------------
# Figures: bootstrap probabilities, forest plots, single-cell pseudobulk
# -----------------------------------------------------------------------------
top_bp <- dplyr::bind_rows(lapply(c("up", "down"), function(dr) head(bp[bp$direction == dr, ], 15)))
top_bp$gene <- factor(top_bp$gene, levels = rev(unique(top_bp$gene)))
p_bp <- ggplot2::ggplot(top_bp, ggplot2::aes(p_top5, gene, fill = core)) +
  ggplot2::geom_col() + ggplot2::geom_point(ggplot2::aes(x = p_top10), shape = 124, size = 4) +
  ggplot2::geom_vline(xintercept = PRI_CORE_PROB, linetype = "dashed", colour = "grey40") +
  ggplot2::facet_wrap(~ direction, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c(`TRUE` = COL_INVASIVE, `FALSE` = "grey70"), name = "Core") +
  ggplot2::labs(x = paste0("Bootstrap probability of top ", PRI_N_PER_DIR, " (bar) and top 10 (tick)"), y = NULL) +
  PLOT_THEME + ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"))
save_plot(p_bp, file.path(DIR_FIG, "11_bootstrap_shortlist_probability"), width = 11, height = 6)

SHOW <- union(SHORT, CORE)
mrow <- match(SHOW, main$gene)
fp <- dplyr::bind_rows(lapply(seq_along(SHOW), function(j) { g <- SHOW[j]; r <- main[mrow[j], ]; dplyr::bind_rows(
  data.frame(gene = g, row = CONTRASTS, g = G[g, ], lo = G[g, ] - Z975 * sqrt(V[g, ]),
             hi = G[g, ] + Z975 * sqrt(V[g, ]), kind = "contrast"),
  data.frame(gene = g, row = c("Pooled (REML, Wald 95% CI)", "Pooled (HKSJ 95% CI)", "95% prediction interval"),
             g = r$g_re, lo = c(r$ci_lo, r$ci_lo_hksj, r$pi_lo), hi = c(r$ci_hi, r$ci_hi_hksj, r$pi_hi),
             kind = c("pooled", "hksj", "pi"))) }))
fp$gene <- factor(fp$gene, levels = SHOW)
fp$row <- factor(fp$row, levels = rev(c(CONTRASTS, "Pooled (REML, Wald 95% CI)", "Pooled (HKSJ 95% CI)",
                                        "95% prediction interval")))
lab <- setNames(paste0(SHOW, ifelse(SHOW %in% CORE, " [core]", ""), "  (I2 = ",
                       round(100 * main$I2[mrow]), "%)"), SHOW)
p_fp <- ggplot2::ggplot(fp, ggplot2::aes(g, row, colour = kind)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0.2, na.rm = TRUE) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_colour_manual(values = c(contrast = "grey30", pooled = COL_INVASIVE, hksj = "#8E44AD", pi = "#E67E22"),
                               guide = "none") +
  ggplot2::facet_wrap(~ gene, ncol = 5, labeller = ggplot2::as_labeller(lab)) +
  ggplot2::labs(x = "Hedges' g, aggressive vs reference", y = NULL) + PLOT_THEME
save_plot(p_fp, file.path(DIR_FIG, "11_forest_shortlist"), width = 16, height = 3.5 * ceiling(length(SHOW) / 5))

sg <- intersect(SHOW, rownames(pb$mean))
dot <- dplyr::bind_rows(lapply(c("Malignant cells", "All tumour-lineage cells"), function(set) {
  mm <- if (set == "Malignant cells") pb$mean else pb$lineage_mean
  pp <- if (set == "Malignant cells") pb$pct else pb$lineage_pct
  dplyr::bind_rows(lapply(sg, function(g) data.frame(cells = set, gene = g, tumour_id = colnames(mm),
                                                     mean = mm[g, ], pct = pp[g, ])))
}))
dot <- dplyr::left_join(dot, pb$tumours[, c("tumour_id", "dataset", "invasion_status")], by = "tumour_id")
dot$tumour_lab <- paste0(dot$tumour_id, ifelse(is.na(dot$invasion_status), "",
                         ifelse(dot$invasion_status == "invasive", " [extraocular]", " [intraocular]")))
dot$gene <- factor(dot$gene, levels = sg)
dot$cells <- factor(dot$cells, levels = c("Malignant cells", "All tumour-lineage cells"))
p_dot <- ggplot2::ggplot(dot, ggplot2::aes(gene, tumour_lab, size = pct, colour = mean)) +
  ggplot2::geom_point() +
  ggplot2::scale_colour_gradient(low = "grey85", high = COL_INVASIVE, name = "Mean\nexpression") +
  ggplot2::scale_size_area(max_size = 6, name = "% cells") +
  ggplot2::facet_grid(dataset ~ cells, scales = "free_y", space = "free_y") +
  ggplot2::labs(x = NULL, y = NULL, title = "Shortlisted and core genes per single-cell tumour") +
  PLOT_THEME + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "italic"))
save_plot(p_dot, file.path(DIR_FIG, "11_scrna_pseudobulk_shortlist"), width = 14, height = 8)

message("\n[11_priority_ranking.R] Complete.")
