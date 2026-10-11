# =============================================================================
# utils.R — helper functions used by the analysis scripts
# =============================================================================

#' Quietly load a package, install if missing
load_pkg <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# -----------------------------------------------------------------------------
# I/O helpers
# -----------------------------------------------------------------------------

#' Save an R object as .rds with a message
save_rds <- function(obj, path, label = NULL) {
  saveRDS(obj, path)
  lbl <- if (!is.null(label)) label else deparse(substitute(obj))
  message("  Saved: ", basename(path), "  [", lbl, "]")
}

#' Save a ggplot as PDF + PNG
save_plot <- function(p, path_no_ext, width = PLOT_WIDTH_WIDE,
                      height = PLOT_HEIGHT, dpi = PLOT_DPI) {
  ggplot2::ggsave(paste0(path_no_ext, ".pdf"), plot = p,
                  width = width, height = height, units = "in")
  ggplot2::ggsave(paste0(path_no_ext, ".png"), plot = p,
                  width = width, height = height, units = "in", dpi = dpi)
  message("  Saved: ", basename(path_no_ext), ".pdf / .png")
}

#' Write a data frame to TSV
save_tsv <- function(df, path) {
  readr::write_tsv(df, path)
  message("  Saved: ", basename(path))
}

#' Remove duplicate gene symbols by keeping the row with highest IQR
dedup_by_iqr <- function(expr_mat, gene_col = "gene_symbol") {
  if (is.data.frame(expr_mat)) {
    genes  <- expr_mat[[gene_col]]
    mat    <- as.matrix(expr_mat[, !names(expr_mat) %in% gene_col])
  } else {
    genes  <- rownames(expr_mat)
    mat    <- expr_mat
  }
  iqrs   <- apply(mat, 1, IQR, na.rm = TRUE)
  best   <- tapply(seq_along(genes), genes, function(idx) idx[which.max(iqrs[idx])])
  mat[unlist(best), , drop = FALSE]
}

#' Standard volcano plot for DEG results
plot_volcano <- function(deg_df,
                         log2fc_col = "log2FC",
                         adjp_col   = "adj.P.Val",
                         gene_col   = "gene",
                         title      = "Volcano Plot",
                         log2fc_cut = 1.0,
                         adjp_cut   = 0.05,
                         label_top  = 20) {
  load_pkg("ggplot2")
  load_pkg("ggrepel", bioc = FALSE)

  df <- deg_df
  df$direction <- dplyr::case_when(
    df[[log2fc_col]] >= log2fc_cut & df[[adjp_col]] <= adjp_cut ~ "Up",
    df[[log2fc_col]] <= -log2fc_cut & df[[adjp_col]] <= adjp_cut ~ "Down",
    TRUE ~ "NS"
  )
  df$neg_log10_p <- -log10(df[[adjp_col]] + 1e-300)

  label_df <- df |>
    dplyr::filter(direction != "NS") |>
    dplyr::arrange(dplyr::desc(abs(.data[[log2fc_col]]))) |>
    dplyr::slice_head(n = label_top)

  ggplot2::ggplot(df, ggplot2::aes(
    x     = .data[[log2fc_col]],
    y     = neg_log10_p,
    color = direction
  )) +
    ggplot2::geom_point(alpha = 0.6, size = 1.2) +
    ggrepel::geom_text_repel(
      data          = label_df,
      ggplot2::aes(label = .data[[gene_col]]),
      size          = 2.8,
      max.overlaps  = 30,
      segment.color = "grey50"
    ) +
    ggplot2::scale_color_manual(
      values = c("Up" = COL_UPREGULATED, "Down" = COL_DOWNREGULATED, "NS" = COL_NONSIG)
    ) +
    ggplot2::geom_vline(xintercept = c(-log2fc_cut, log2fc_cut),
                        linetype = "dashed", color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(adjp_cut),
                        linetype = "dashed", color = "grey40") +
    ggplot2::labs(
      title = title,
      x     = expression(log[2]~"Fold Change"),
      y     = expression(-log[10]~"(adj. p-value)"),
      color = "Direction"
    ) +
    PLOT_THEME
}

# -----------------------------------------------------------------------------
# Heatmap helper
# -----------------------------------------------------------------------------

#' Scaled expression heatmap (ComplexHeatmap)
plot_heatmap <- function(expr_mat, annotation_df = NULL, title = "",
                         show_rownames = TRUE, max_genes = 100) {
  load_pkg("ComplexHeatmap", bioc = TRUE)
  load_pkg("circlize", bioc = FALSE)

  if (nrow(expr_mat) > max_genes) {
    iqrs <- apply(expr_mat, 1, IQR)
    expr_mat <- expr_mat[order(iqrs, decreasing = TRUE)[seq_len(max_genes)], ]
  }
  scaled <- t(scale(t(expr_mat)))

  col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#377EB8", "white", "#E41A1C"))

  ha <- NULL
  if (!is.null(annotation_df)) {
    ha <- ComplexHeatmap::HeatmapAnnotation(df = annotation_df)
  }

  ComplexHeatmap::Heatmap(
    scaled,
    name              = "z-score",
    col               = col_fun,
    top_annotation    = ha,
    show_row_names    = show_rownames,
    show_column_names = ncol(expr_mat) <= 60,
    row_names_gp      = grid::gpar(fontsize = 7),
    column_title      = title,
    clustering_method_rows    = "ward.D2",
    clustering_method_columns = "ward.D2"
  )
}

# -----------------------------------------------------------------------------
# Signed, sample-size-weighted Stouffer meta-analysis of per-dataset limma tables.
# res_list: named list of data.frames with gene, log2FC, P.Value (raw, two-sided)
# n       : named numeric vector, samples per dataset (weights = sqrt(n))
# Raw p-values are converted to signed z (positive = higher in invasive), so a
# gene that is up in some datasets and down in others cancels instead of
# reaching significance (Fisher on two-sided p-values does not do this).
# -----------------------------------------------------------------------------
meta_stouffer <- function(res_list, n, labels = c("UP_invasive", "DOWN_invasive")) {
  stopifnot(all(names(res_list) %in% names(n)))
  genes <- unique(unlist(lapply(res_list, `[[`, "gene")))
  z <- fc <- matrix(NA_real_, length(genes), length(res_list),
                    dimnames = list(genes, names(res_list)))
  for (nm in names(res_list)) {
    d <- res_list[[nm]]
    i <- match(d$gene, genes)
    z[i, nm]  <- sign(d$log2FC) * stats::qnorm(pmax(d$P.Value, 1e-300) / 2, lower.tail = FALSE)
    fc[i, nm] <- d$log2FC
  }
  w      <- sqrt(n[colnames(z)])
  w_mat  <- sweep(!is.na(z), 2, w, `*`)           # weight only where the gene was measured
  z_meta <- rowSums(sweep(z, 2, w, `*`), na.rm = TRUE) / sqrt(rowSums(w_mat^2))
  p      <- 2 * stats::pnorm(-abs(z_meta))
  out <- data.frame(
    gene        = genes,
    mean_log2FC = rowMeans(fc, na.rm = TRUE),
    meta_z      = z_meta,
    direction   = ifelse(z_meta >= 0, labels[1], labels[2]),
    n_datasets  = rowSums(!is.na(z)),
    meta_p      = p,
    meta_padj   = stats::p.adjust(p, method = "BH"),
    stringsAsFactors = FALSE, row.names = NULL
  )
  out[order(out$meta_padj), ]
}

# Self-check: opposite effects cancel; concordant effects combine.
local({
  a <- data.frame(gene = c("g1", "g2"), log2FC = c(1, 1),  P.Value = c(1e-4, 1e-4))
  b <- data.frame(gene = c("g1", "g2"), log2FC = c(1, -1), P.Value = c(1e-4, 1e-4))
  m <- meta_stouffer(list(a = a, b = b), c(a = 10, b = 10))
  stopifnot(m$meta_p[m$gene == "g1"] < 1e-6, m$meta_p[m$gene == "g2"] > 0.99)
})

# -----------------------------------------------------------------------------
# Invasion labels used by steps 04, 05, 08 and 09 — one definition so the steps agree.
# Return "invasive" / "non-invasive" / NA per row of a GEO pheno table.
# -----------------------------------------------------------------------------
label_invasion_gse97508 <- function(p) dplyr::case_when(
  grepl("non.?invasi", p$title, ignore.case = TRUE) ~ "non-invasive",
  grepl("invasi",      p$title, ignore.case = TRUE) ~ "invasive",
  TRUE ~ NA_character_            # adult retina controls
)
label_invasion_gse110811 <- function(p) dplyr::case_when(
  grepl("severe", p$characteristics_ch1.1, ignore.case = TRUE) ~ "invasive",
  grepl("mild",   p$characteristics_ch1.1, ignore.case = TRUE) ~ "non-invasive",
  TRUE ~ NA_character_            # Moderate (13) and unannotated (normal retina, retinocytoma) excluded
)
label_invasion_gse59983 <- function(p) dplyr::case_when(
  # GEO "uhc-subgroup: group 2" (n = 4) is Kooi 2015's small Ward cluster that the authors
  # excluded as likely normal-retina contamination, so it is not used as a reference group.
  grepl("group 1", p$characteristics_ch1.1, ignore.case = TRUE) ~ "invasive",
  grepl("group 3", p$characteristics_ch1.1, ignore.case = TRUE) ~ "non-invasive",
  TRUE ~ NA_character_
)

# -----------------------------------------------------------------------------
# Standardised effect size from a two-group ordinary (unmoderated) t statistic.
# Hedges' g with small-sample correction J and its sampling variance.
# -----------------------------------------------------------------------------
hedges_g <- function(t, n1, n2) {
  J <- 1 - 3 / (4 * (n1 + n2) - 9)
  g <- J * t * sqrt(1 / n1 + 1 / n2)
  data.frame(g = g, v = 1 / n1 + 1 / n2 + g^2 / (2 * (n1 + n2)))
}

# -----------------------------------------------------------------------------
# REML estimate of tau^2 per row of G (effects) / V (variances), vectorised.
# Restricted log-likelihood is evaluated on a grid (0 and 120 log-spaced points up to
# U = max(1, 10 x var(g))), then refined by golden-section search around the best grid
# point; the best of the bracket ends and the refined point is returned.
# -----------------------------------------------------------------------------
reml_tau2 <- function(G, V, k) {
  ll <- function(t2) {
    W <- 1 / (V + t2); s1 <- rowSums(W, na.rm = TRUE); mu <- rowSums(W * G, na.rm = TRUE) / s1
    -0.5 * (rowSums(log(V + t2), na.rm = TRUE) + log(s1) + rowSums(W * (G - mu)^2, na.rm = TRUE))
  }
  vg <- apply(G, 1, stats::var, na.rm = TRUE); vg[!is.finite(vg)] <- 0
  U <- pmax(1, 10 * vg)
  grid <- c(0, exp(seq(log(1e-6), 0, length.out = 120)))
  L <- vapply(grid, function(f) ll(f * U), numeric(nrow(G)))
  L[!is.finite(L)] <- -Inf
  j <- max.col(L, ties.method = "first")
  lo <- grid[pmax(j - 1, 1)] * U; hi <- grid[pmin(j + 1, length(grid))] * U
  a <- lo; b <- hi; r <- (sqrt(5) - 1) / 2
  for (i in 1:80) {
    c1 <- b - r * (b - a); d1 <- a + r * (b - a); left <- ll(c1) > ll(d1)
    a <- ifelse(left, a, c1); b <- ifelse(left, d1, b)
  }
  cand <- cbind(lo, hi, (a + b) / 2)
  Lc <- cbind(ll(lo), ll(hi), ll((a + b) / 2)); Lc[!is.finite(Lc)] <- -Inf
  t2 <- cand[cbind(seq_len(nrow(G)), max.col(Lc, ties.method = "first"))]
  ifelse(k > 1, t2, 0)
}

# -----------------------------------------------------------------------------
# Random-effects meta-analysis, vectorised over genes.
# res_list: named list of data.frames with gene, g, v (one per contrast).
# tau^2 by REML, maximised per gene (reml_tau2: grid + golden section; the REML
# likelihood can have several stationary points, so fixed-point iteration is not used).
# Discovery inference: Wald z on the REML pooled g (p_re, padj_re).
# Also returned: Hartung-Knapp-Sidik-Jonkman (HKSJ) SE/p/CI with k-1 df and the
# 95% prediction interval (t with k-2 df, k >= 3). BH is applied only to genes
# measured in >= 2 contrasts (single-contrast genes cannot be meta-analysed).
# -----------------------------------------------------------------------------
meta_re <- function(res_list) {
  genes <- unique(unlist(lapply(res_list, `[[`, "gene")))
  G <- V <- matrix(NA_real_, length(genes), length(res_list),
                   dimnames = list(genes, names(res_list)))
  for (nm in names(res_list)) {
    d <- res_list[[nm]]; i <- match(d$gene, genes)
    G[i, nm] <- d$g; V[i, nm] <- d$v
  }
  G[!is.finite(G) | !is.finite(V)] <- NA; V[is.na(G)] <- NA
  k   <- rowSums(!is.na(G))
  w   <- 1 / V
  sw  <- rowSums(w, na.rm = TRUE)
  gfe <- rowSums(w * G, na.rm = TRUE) / sw
  Q   <- rowSums(w * (G - gfe)^2, na.rm = TRUE)
  tau2 <- reml_tau2(G, V, k)
  ws    <- 1 / (V + tau2); sws <- rowSums(ws, na.rm = TRUE)
  g_re  <- rowSums(ws * G, na.rm = TRUE) / sws
  se_re <- sqrt(1 / sws)
  p_re  <- 2 * stats::pnorm(-abs(g_re / se_re))
  q_h   <- ifelse(k > 1, rowSums(ws * (G - g_re)^2, na.rm = TRUE) / (k - 1), NA_real_)
  se_h  <- sqrt(q_h / sws)
  df_h  <- ifelse(k > 1, k - 1, NA_real_)
  p_h   <- 2 * stats::pt(-abs(g_re / se_h), df_h)
  t_h   <- stats::qt(0.975, df_h)
  t_pi  <- ifelse(k > 2, stats::qt(0.975, pmax(k - 2, 1)), NA_real_)
  multi <- k >= 2
  bh    <- function(p) { out <- rep(NA_real_, length(p)); out[multi] <- stats::p.adjust(p[multi], "BH"); out }
  out <- data.frame(
    gene = genes, g_re = g_re, se_re = se_re, z_re = g_re / se_re, p_re = p_re, padj_re = bh(p_re),
    ci_lo = g_re - stats::qnorm(0.975) * se_re, ci_hi = g_re + stats::qnorm(0.975) * se_re,
    se_hksj = se_h, p_hksj = p_h, padj_hksj = bh(p_h),
    ci_lo_hksj = g_re - t_h * se_h, ci_hi_hksj = g_re + t_h * se_h,
    pi_lo = g_re - t_pi * sqrt(tau2 + se_re^2), pi_hi = g_re + t_pi * sqrt(tau2 + se_re^2),
    tau2 = tau2, I2 = ifelse(k > 1, ifelse(Q > 0, pmax(0, (Q - (k - 1)) / Q), 0), NA_real_), Q = Q,
    k = k, n_agree = rowSums(sign(G) == sign(g_re), na.rm = TRUE),
    direction = ifelse(g_re >= 0, "UP_invasive", "DOWN_invasive"),
    stringsAsFactors = FALSE, row.names = NULL
  )
  out[order(out$p_re), ]
}

# Self-check: identical effects -> no heterogeneity; opposite effects -> no pooled effect,
# and for k = 2 with equal v the REML tau^2 has the closed form (g1 - g2)^2 / 2 - v.
local({
  a <- data.frame(gene = c("g1", "g2", "g3"), g = c(1, 1, 1),  v = c(0.1, 0.1, 0.1))
  b <- data.frame(gene = c("g1", "g2"),       g = c(1, -1),    v = c(0.1, 0.1))
  m <- meta_re(list(a = a, b = b))
  r1 <- m[m$gene == "g1", ]; r2 <- m[m$gene == "g2", ]; r3 <- m[m$gene == "g3", ]
  stopifnot(abs(r1$g_re - 1) < 1e-12, r1$tau2 == 0, r1$I2 == 0, r1$n_agree == 2,
            abs(r2$g_re) < 1e-12, abs(r2$tau2 - 1.9) < 1e-6, r2$p_re > 0.99,
            is.na(r3$padj_re), r3$k == 1)
})

# -----------------------------------------------------------------------------
# Robust rank aggregation (Kolde et al. 2012) on normalised ranks in (0, 1].
# rmat: genes x lists matrix of normalised ranks (NA = not measured).
# Returns the Bonferroni-corrected rho score per gene (smaller = consistently high).
# -----------------------------------------------------------------------------
rra_rho <- function(rmat) {
  apply(rmat, 1, function(r) {
    r <- sort(r[!is.na(r)]); n <- length(r)
    if (n == 0) return(NA_real_)
    min(1, n * min(stats::pbeta(r, seq_len(n), n - seq_len(n) + 1)))
  })
}
stopifnot(rra_rho(matrix(c(0.01, 0.02, 0.9, 0.95), 2, byrow = TRUE,
                         dimnames = list(c("top", "bottom"), NULL)))[["top"]] < 0.01)

# -----------------------------------------------------------------------------
# Signed phenotype score: mean z(up genes) - mean z(down genes), genes z-scored
# across the samples of one cohort. Leave-one-gene-out version for gene g:
# the gene is removed from its own set (exact, no re-scoring needed).
# -----------------------------------------------------------------------------
zscore_rows <- function(mat) t(scale(t(mat)))
score_updown <- function(Z, up, down) {
  up <- intersect(up, rownames(Z)); down <- intersect(down, rownames(Z))
  colMeans(Z[up, , drop = FALSE]) - colMeans(Z[down, , drop = FALSE])
}
score_logo <- function(Z, up, down, g) score_updown(Z, setdiff(up, g), setdiff(down, g))

# Partial Spearman correlation of x and y given z (ranks, residuals, df = n - 3)
partial_spearman <- function(x, y, z) {
  ok <- stats::complete.cases(x, y, z); n <- sum(ok)
  if (n < 6) return(c(rho = NA_real_, p = NA_real_))
  rx <- stats::resid(stats::lm(rank(x[ok]) ~ rank(z[ok])))
  ry <- stats::resid(stats::lm(rank(y[ok]) ~ rank(z[ok])))
  r  <- stats::cor(rx, ry)
  tt <- r * sqrt((n - 3) / (1 - r^2))
  c(rho = r, p = 2 * stats::pt(-abs(tt), df = n - 3))
}

# AUC with DeLong 95% CI; under complete separation (AUC 0 or 1) the DeLong variance is
# zero, so the interval is reported as not estimable (NA) instead of a degenerate 1.00-1.00.
auc_ci <- function(roc) {
  ci <- suppressWarnings(as.numeric(pROC::ci.auc(roc, method = "delong")))
  if (ci[2] %in% c(0, 1)) return(c(auc = ci[2], lo = NA_real_, hi = NA_real_))
  c(auc = ci[2], lo = max(0, ci[1]), hi = min(1, ci[3]))
}

# Sample annotation from a GEO series-matrix file without parsing its data table.
# Returns title, gsm, platform, source and one column per "key: value" characteristic.
read_series_pheno <- function(f) {
  L <- readLines(gzfile(f))
  field <- function(tag) lapply(grep(paste0("^!", tag, "\t"), L, value = TRUE),
                                function(l) gsub('^"|"$', "", strsplit(l, "\t")[[1]][-1]))
  ch <- field("Sample_characteristics_ch1")
  p <- data.frame(title = field("Sample_title")[[1]], gsm = field("Sample_geo_accession")[[1]],
                  platform = field("Sample_platform_id")[[1]], source = field("Sample_source_name_ch1")[[1]],
                  stringsAsFactors = FALSE)
  kv <- lapply(seq_len(nrow(p)), function(j) {
    x <- vapply(ch, `[`, character(1), j); x <- x[!is.na(x) & nzchar(x)]
    setNames(trimws(sub("^[^:]*:", "", x)), make.names(trimws(sub(":.*$", "", x))))
  })
  for (k in unique(unlist(lapply(kv, names)))) p[[k]] <- vapply(kv, function(v) if (k %in% names(v)) unname(v[[k]]) else NA_character_, "")
  p
}
