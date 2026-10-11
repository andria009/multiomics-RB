# =============================================================================
# Step 05 — Exploratory classifier (not used to rank or select genes)
#
# 1. Stability selection: LASSO (glmnet, lambda.1se) on GSE29683 + GSE59983 in
#    N_BOOT bootstrap resamples (stratified by cohort and label); selection
#    frequency per gene. Each cohort is z-scored gene-wise on its own, so no
#    label-informed batch correction is used.
# 2. Leave-one-cohort-out (LOCO) validation: for each held-out microarray cohort,
#    candidate features come from a meta-analysis WITHOUT that cohort, the model is
#    trained on the other labelled cohorts, and the held-out AUC is reported with a
#    DeLong 95% CI. GSE97508 (3 vs 3) is reported but cannot test anything.
# Output: results/02_invasion/lasso_stability.tsv, loco_auc.tsv; figures/05_*
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
if (length(INVASION_EXCLUDE) > 0) stop("Step 05 runs on the main results only. Unset INVASION_EXCLUDE.")

load_pkg("glmnet", bioc = FALSE)
load_pkg("pROC",   bioc = FALSE)
load_pkg("dplyr")
load_pkg("ggplot2")

set.seed(GLOBAL_SEED)

normalized  <- readRDS(file.path(DIR_BULK, "normalized_matrices.rds"))
inv_results <- readRDS(file.path(DIR_INV, "invasion_deg_per_dataset.rds"))
cand        <- readRDS(file.path(DIR_INV, "invasion_candidates.rds"))
lab_tbl     <- read.delim(file.path(DIR_INV, "sample_labels.tsv"), stringsAsFactors = FALSE)

COHORTS <- c("GSE29683", "GSE59983", "GSE110811", "GSE97508")
cohort_xy <- lapply(setNames(COHORTS, COHORTS), function(ds) {
  l <- lab_tbl[lab_tbl$dataset == ds & !is.na(lab_tbl$label), ]
  list(Z = zscore_rows(normalized[[ds]][, l$sample_id, drop = FALSE]),
       y = as.integer(l$label == "invasive"))
})

candidate_genes <- function(res_list) {
  m <- meta_re(res_list)
  m$gene[m$padj_re < META_ADJP & abs(m$g_re) >= META_G_MIN & m$k >= 2]
}
stack_xy <- function(dss, genes) {
  X <- do.call(rbind, lapply(dss, function(ds) t(cohort_xy[[ds]]$Z[genes, , drop = FALSE])))
  list(X = X, y = unlist(lapply(dss, function(ds) cohort_xy[[ds]]$y)),
       cohort = rep(dss, vapply(dss, function(ds) length(cohort_xy[[ds]]$y), integer(1))))
}
fit_lasso <- function(X, y) {
  nf <- max(3, min(ML_CV_FOLDS, min(table(y))))
  glmnet::cv.glmnet(X, y, family = "binomial", alpha = LASSO_ALPHA, nfolds = nf,
                    type.measure = "auc", standardize = TRUE)
}
complete_genes <- function(dss, genes) {
  ok <- Reduce(intersect, lapply(dss, function(ds) rownames(cohort_xy[[ds]]$Z)))
  g <- intersect(genes, ok)
  g[vapply(g, function(x) all(vapply(dss, function(ds) all(is.finite(cohort_xy[[ds]]$Z[x, ])), logical(1))), logical(1))]
}

# -----------------------------------------------------------------------------
# 1. Stability selection on the two large cohorts
# -----------------------------------------------------------------------------
TRAIN <- c("GSE29683", "GSE59983")
feat  <- complete_genes(TRAIN, cand$gene)
tr    <- stack_xy(TRAIN, feat)
message("Stability selection: ", nrow(tr$X), " tumours, ", length(feat), " candidate features, ",
        N_BOOT, " bootstraps")

strata <- paste(tr$cohort, tr$y)
sel <- matrix(0, length(feat), N_BOOT, dimnames = list(feat, NULL))
coefs <- sel
for (b in seq_len(N_BOOT)) {
  idx <- unlist(lapply(split(seq_along(strata), strata), function(i) sample(i, length(i), replace = TRUE)))
  cv  <- fit_lasso(tr$X[idx, ], tr$y[idx])
  cf  <- as.matrix(coef(cv, s = "lambda.1se"))[feat, 1]
  sel[, b] <- cf != 0; coefs[, b] <- cf
}
stab <- data.frame(gene = feat, selection_freq = rowMeans(sel),
                   mean_coef_when_selected = ifelse(rowSums(sel) > 0,
                     rowSums(coefs) / pmax(rowSums(sel), 1), NA_real_)) |>
  dplyr::arrange(dplyr::desc(selection_freq))
save_tsv(stab, file.path(DIR_INV, "lasso_stability.tsv"))
message("  Genes selected in >= 50% of bootstraps: ", sum(stab$selection_freq >= 0.5))

top30 <- head(stab, 30)
p_stab <- ggplot2::ggplot(top30, ggplot2::aes(reorder(gene, selection_freq), selection_freq,
                                              fill = mean_coef_when_selected > 0)) +
  ggplot2::geom_col() + ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = c(`TRUE` = COL_INVASIVE, `FALSE` = COL_NONINVASIVE),
                             labels = c(`TRUE` = "higher in aggressive group", `FALSE` = "lower"), name = NULL) +
  ggplot2::labs(title = "LASSO selection frequency (GSE29683 + GSE59983)",
                subtitle = paste0(N_BOOT, " stratified bootstrap resamples; top 30 genes"),
                x = NULL, y = "Selection frequency") + PLOT_THEME
save_plot(p_stab, file.path(DIR_FIG, "05_lasso_stability"), width = 7, height = 7)

# -----------------------------------------------------------------------------
# 2. Leave-one-cohort-out validation
# -----------------------------------------------------------------------------
LABELLED_TRAIN <- c("GSE29683", "GSE59983", "GSE110811")
loco <- dplyr::bind_rows(lapply(COHORTS, function(h) {
  train_ds <- setdiff(LABELLED_TRAIN, h)
  key_drop <- h                                    # contrast names match cohort IDs
  genes <- complete_genes(c(train_ds, h), candidate_genes(inv_results[setdiff(names(inv_results), key_drop)]))
  trn <- stack_xy(train_ds, genes)
  cv  <- fit_lasso(trn$X, trn$y)
  tst <- stack_xy(h, genes)
  pr  <- as.numeric(predict(cv, newx = tst$X, s = "lambda.1se", type = "link"))
  roc <- pROC::roc(tst$y, pr, levels = c(0, 1), direction = "<", quiet = TRUE)
  ci  <- unname(auc_ci(roc))[c(2, 1, 3)]          # lo, auc, hi (NA bounds if separation is complete)
  n_sel <- sum(as.matrix(coef(cv, s = "lambda.1se"))[-1, 1] != 0)
  message("  LOCO held-out ", h, ": AUC ", round(ci[2], 3), " (", round(ci[1], 2), "-", round(ci[3], 2),
          "); trained on ", paste(train_ds, collapse = "+"), ", ", length(genes), " features, ",
          n_sel, " selected")
  data.frame(held_out = h, trained_on = paste(train_ds, collapse = " + "),
             n_inv = sum(tst$y == 1), n_noninv = sum(tst$y == 0), n_features = length(genes),
             n_selected = n_sel, auc = ci[2], auc_lo = ci[1], auc_hi = ci[3])
}))
save_tsv(loco, file.path(DIR_INV, "loco_auc.tsv"))

p_loco <- ggplot2::ggplot(loco, ggplot2::aes(auc, reorder(held_out, auc))) +
  ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = auc_lo, xmax = auc_hi), height = 0.2) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(ggplot2::aes(label = paste0(n_inv, " vs ", n_noninv)), nudge_y = 0.3, size = 3) +
  ggplot2::coord_cartesian(xlim = c(0, 1)) +
  ggplot2::labs(title = "Leave-one-cohort-out LASSO", x = "Held-out AUC (DeLong 95% CI; none when separation is complete)", y = NULL,
                subtitle = "Features selected without the held-out cohort") + PLOT_THEME
save_plot(p_loco, file.path(DIR_FIG, "05_loco_auc"), width = 7, height = 4)

message("\n[05_invasion_signature.R] Complete.")
