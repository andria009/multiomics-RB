# =============================================================================
# Step 10 — Benchmark against published RB subtype / differentiation signatures
# Input : data/supplementary/41467_2021_25792_MOESM6_ESM.xlsx (Liu 2021 Supplementary
#         Data 3: subtype 2 vs 1, Institut Curie cohort) and data/supplementary/mmc10.xlsx
#         (Kooi 2015 Table S3: Ward cluster 1 vs 2 in GSE59983, with gene clusters),
#         results/01_bulk, results/02_invasion
# Signatures (defined from the published tables before they were scored here):
#   LIU2021_subtype2  up = logFC >= 1, BH < 0.05 (higher in subtype 2); down = logFC <= -1, BH < 0.05
#   KOOI2015_low_photoreceptorness  up = Kooi gene cluster 1 (M-phase/ribosome, high in
#     Ward cluster 1); down = gene cluster 3 (photoreceptor, defines "photoreceptorness").
#     Derived from GSE59983 itself, so it is circular for the GSE59983 contrast.
# Output: results/06_subtype/
#   subtype_scores.tsv            per-sample signature scores (labelled primary tumours)
#   subtype_label_auc.tsv         how well each signature separates each contrast's labels
#   subtype_adjusted_effects.tsv  per candidate gene: label effect after adjusting for
#                                 each signature score (ANCOVA; pooled beta with 95% CI,
#                                 unadjusted beta, Stouffer p across contrasts)
#   subtype_label_collinearity.tsv  label-score correlation and VIF per contrast
#   figures/10_subtype_vs_labels
# A gene's own membership in a signature is removed before adjusting (no part-whole).
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
if (length(INVASION_EXCLUDE) > 0) stop("Unset INVASION_EXCLUDE.")

load_pkg("dplyr")
load_pkg("tidyr")
load_pkg("ggplot2")
load_pkg("pROC", bioc = FALSE)

load_pkg("readxl", bioc = FALSE)
liu_f  <- here::here("data", "supplementary", "41467_2021_25792_MOESM6_ESM.xlsx")
kooi_f <- here::here("data", "supplementary", "mmc10.xlsx")
for (f in c(liu_f, kooi_f)) if (!file.exists(f)) stop("Missing ", f, " (see README.md, Data files).")
liu <- readxl::read_excel(liu_f, sheet = "genes differentially expressed", skip = 3)
liu <- data.frame(gene = liu[[1]], logFC = liu[[2]], padj = liu[[6]])
kooi <- readxl::read_excel(kooi_f)
kooi <- data.frame(gene = kooi[["Symbol"]], cluster = substr(as.character(kooi[["Gene cluster"]]), 1, 1))
k1 <- unique(na.omit(kooi$gene[kooi$cluster == "1"])); k3 <- unique(na.omit(kooi$gene[kooi$cluster == "3"]))
both <- intersect(k1, k3)
sigs <- rbind(
  data.frame(signature = "LIU2021_subtype2", gene = unique(liu$gene[liu$padj < 0.05 & liu$logFC >= 1]), set = "up"),
  data.frame(signature = "LIU2021_subtype2", gene = unique(liu$gene[liu$padj < 0.05 & liu$logFC <= -1]), set = "down"),
  data.frame(signature = "KOOI2015_low_photoreceptorness", gene = setdiff(k1, both), set = "up"),
  data.frame(signature = "KOOI2015_low_photoreceptorness", gene = setdiff(k3, both), set = "down"))
sigs$source <- ifelse(grepl("LIU", sigs$signature), "Liu 2021 Nat Commun, Supplementary Data 3",
                      "Kooi 2015 EBioMedicine, Table S3")
print(table(sigs$signature, sigs$set))
save_tsv(sigs, file.path(DIR_SUB, "signatures.tsv"))
SIGS <- unique(sigs$signature)

normalized  <- readRDS(file.path(DIR_BULK, "normalized_matrices.rds"))
rnaseq_mat  <- readRDS(file.path(DIR_BULK, "rnaseq_log2tpm_GSE268136.rds"))
inv_results <- readRDS(file.path(DIR_INV, "invasion_deg_per_dataset.rds"))
cand        <- readRDS(file.path(DIR_INV, "invasion_candidates.rds"))
lab_tbl     <- read.delim(file.path(DIR_INV, "sample_labels.tsv"), stringsAsFactors = FALSE)
n_eff       <- attr(inv_results, "n")
CONTRASTS   <- names(inv_results)
mat_of <- function(ds) if (ds == "GSE268136_ClusterC") rnaseq_mat else normalized[[ds]]

cohorts <- lapply(setNames(CONTRASTS, CONTRASTS), function(ds) {
  l <- lab_tbl[lab_tbl$dataset == ds & !is.na(lab_tbl$label), ]
  l <- l[l$sample_id %in% colnames(mat_of(ds)), ]
  Z <- zscore_rows(mat_of(ds)[, l$sample_id, drop = FALSE]); Z <- Z[stats::complete.cases(Z), , drop = FALSE]
  list(Z = Z, y = as.integer(l$label == "invasive"))
})
sig_sets <- lapply(setNames(SIGS, SIGS), function(s)
  list(up = sigs$gene[sigs$signature == s & sigs$set == "up"],
       down = sigs$gene[sigs$signature == s & sigs$set == "down"]))
sig_score <- function(Z, s, drop = character(0)) {
  up <- setdiff(intersect(sig_sets[[s]]$up, rownames(Z)), drop)
  dn <- setdiff(intersect(sig_sets[[s]]$down, rownames(Z)), drop)
  colMeans(Z[up, , drop = FALSE]) - if (length(dn)) colMeans(Z[dn, , drop = FALSE]) else 0
}

# -----------------------------------------------------------------------------
# 1. Signature scores and their separation of each contrast's labels
# -----------------------------------------------------------------------------
scores <- dplyr::bind_rows(lapply(CONTRASTS, function(ds) {
  co <- cohorts[[ds]]
  data.frame(dataset = ds, sample_id = colnames(co$Z), label = ifelse(co$y == 1, "invasive", "non-invasive"),
             sapply(SIGS, function(s) sig_score(co$Z, s)), check.names = FALSE)
}))
save_tsv(scores, file.path(DIR_SUB, "subtype_scores.tsv"))

auc_tbl <- dplyr::bind_rows(lapply(CONTRASTS, function(ds) dplyr::bind_rows(lapply(SIGS, function(s) {
  co <- cohorts[[ds]]
  cov <- c(up = length(intersect(sig_sets[[s]]$up, rownames(co$Z))),
           down = length(intersect(sig_sets[[s]]$down, rownames(co$Z))))
  r <- pROC::roc(co$y, sig_score(co$Z, s), levels = c(0, 1), direction = "<", quiet = TRUE)
  ci <- unname(auc_ci(r))[c(2, 1, 3)]
  data.frame(dataset = ds, signature = s, genes_up = cov[["up"]], genes_down = cov[["down"]],
             auc = ci[2], auc_lo = ci[1], auc_hi = ci[3])
}))))
save_tsv(auc_tbl, file.path(DIR_SUB, "subtype_label_auc.tsv"))
print(auc_tbl)

long <- tidyr::pivot_longer(scores, dplyr::all_of(SIGS), names_to = "signature", values_to = "score")
p_sub <- ggplot2::ggplot(long, ggplot2::aes(label, score, fill = label)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) + ggplot2::geom_jitter(width = 0.12, size = 0.8) +
  ggplot2::facet_grid(signature ~ dataset, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c(invasive = COL_INVASIVE, "non-invasive" = COL_NONINVASIVE)) +
  ggplot2::scale_x_discrete(labels = c(invasive = "aggr.", "non-invasive" = "ref.")) +
  ggplot2::labs(x = NULL, y = "Signature score", title = "Published RB subtype/differentiation scores by contrast label") +
  PLOT_THEME + ggplot2::theme(legend.position = "none")
save_plot(p_sub, file.path(DIR_FIG, "10_subtype_vs_labels"), width = 14, height = 6)

# -----------------------------------------------------------------------------
# 2. Candidate genes: label effect adjusted for each signature (ANCOVA per contrast:
#    gene_z ~ label + signature score). Per contrast the label coefficient is a mean
#    difference in SD units; it is pooled across contrasts by inverse variance (with 95%
#    CI) next to the unadjusted coefficient from the same samples. Stouffer (n_eff) on the
#    adjusted t gives the p-value. GSE97508 (3 vs 3) is not used: with six tumours the
#    label and the scores cannot be separated. Label-score collinearity is reported as VIF.
# -----------------------------------------------------------------------------
ADJ_SETS <- CONTRASTS[vapply(CONTRASTS, function(ds) length(cohorts[[ds]]$y) >= 15, logical(1))]
vif_tbl <- dplyr::bind_rows(lapply(ADJ_SETS, function(ds) dplyr::bind_rows(lapply(SIGS, function(s) {
  r <- stats::cor(cohorts[[ds]]$y, sig_score(cohorts[[ds]]$Z, s))
  data.frame(dataset = ds, signature = s, r_label_score = r, vif = 1 / (1 - r^2))
}))))
save_tsv(vif_tbl, file.path(DIR_SUB, "subtype_label_collinearity.tsv"))
print(vif_tbl)
pool_iv <- function(b, se) {
  w <- 1 / se^2; ok <- is.finite(w) & !is.na(b)
  est <- rowSums(b * w * ok, na.rm = TRUE) / rowSums(w * ok, na.rm = TRUE)
  sep <- sqrt(1 / rowSums(w * ok, na.rm = TRUE))
  cbind(est = est, lo = est - 1.96 * sep, hi = est + 1.96 * sep)
}
adj <- dplyr::bind_rows(lapply(SIGS, function(s) {
  per <- lapply(ADJ_SETS, function(ds) {
    co <- cohorts[[ds]]; base <- sig_score(co$Z, s)
    genes <- intersect(cand$gene, rownames(co$Z))
    members <- c(sig_sets[[s]]$up, sig_sets[[s]]$down)
    t(vapply(genes, function(g) {
      sc <- if (g %in% members) sig_score(co$Z, s, drop = g) else base
      f  <- summary(stats::lm(co$Z[g, ] ~ co$y + sc))$coefficients
      f0 <- summary(stats::lm(co$Z[g, ] ~ co$y))$coefficients
      c(t = f[2, 3], p = f[2, 4], b = f[2, 1], se = f[2, 2], b0 = f0[2, 1], se0 = f0[2, 2])
    }, numeric(6)))
  })
  names(per) <- ADJ_SETS
  genes <- unique(unlist(lapply(per, rownames)))
  col <- function(nm) sapply(ADJ_SETS, function(ds) per[[ds]][match(genes, rownames(per[[ds]])), nm])
  zmat <- sign(col("t")) * stats::qnorm(pmax(col("p"), 1e-300) / 2, lower.tail = FALSE)
  w <- sqrt(n_eff[ADJ_SETS])
  zc <- rowSums(sweep(zmat, 2, w, `*`), na.rm = TRUE) / sqrt(rowSums(sweep(!is.na(zmat), 2, w^2, `*`)))
  pa <- pool_iv(col("b"), col("se")); pu <- pool_iv(col("b0"), col("se0"))
  data.frame(gene = genes, signature = s, z_adjusted = zc, p_adjusted = 2 * stats::pnorm(-abs(zc)),
             beta_unadjusted = pu[, "est"], beta_unadjusted_lo = pu[, "lo"], beta_unadjusted_hi = pu[, "hi"],
             beta_adjusted = pa[, "est"], beta_adjusted_lo = pa[, "lo"], beta_adjusted_hi = pa[, "hi"],
             n_contrasts = rowSums(!is.na(zmat)))
})) |>
  dplyr::group_by(signature) |> dplyr::mutate(padj_adjusted = p.adjust(p_adjusted, "BH")) |> dplyr::ungroup()
save_tsv(adj, file.path(DIR_SUB, "subtype_adjusted_effects.tsv"))
for (s in SIGS) message("  Candidates still associated (BH < 0.05, same sign) after adjusting for ", s, ": ",
  sum(adj$signature == s & adj$padj_adjusted < 0.05 &
      sign(adj$z_adjusted) == sign(cand$g_re[match(adj$gene, cand$gene)]), na.rm = TRUE),
  " of ", sum(adj$signature == s))

message("\n[10_subtype_benchmark.R] Complete.")
