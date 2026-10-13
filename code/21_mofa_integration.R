# =============================================================================
# Step 21 — Multi-omics factor analysis (MOFA2) of expression, methylation and copy number
# Input : discovery  results/01_bulk/rnaseq_log2tpm_GSE268136.rds, limma_clusterC_vs_AB_GSE268136.tsv (step 03)
#                    results/08_methylation/meth_beta.rds, meth_samples.tsv, meth_phenotype_samples.tsv (13-14)
#         held out  results/10_liu_validation/liu_expression.rds, liu_beta_top20k.rds, liu_samples.tsv,
#                    liu_marker_values.tsv (step 18)
#         both      results/11_cnv/cnv_samples.tsv (step 20)
#         results/07_priority/priority_all_candidates.tsv (step 11)
# Output: results/12_mofa/
#   mofa_variance_explained.tsv  per cohort, factor and view: % variance explained
#   mofa_factors.tsv             factor values per tumour (oriented, see below)
#   mofa_factor_tests.tsv        factors vs phenotype score, Ryl cluster, platform, FGA (discovery) and
#                                vs subtype and high-risk histology (held out)
#   mofa_top_weights.tsv         top 25 features per view for the main factor of each cohort
#   mofa_shortlist_weights.tsv   RNA-view weights (and percentile) of the prioritised genes
#   mofa_stability.tsv           main factor recovered with other seeds (|r| with the seed-42 factor)
#   mofa_cross_cohort.tsv        agreement of the two cohorts' main factors (RNA-view gene weights)
#   mofa_rna_weights_shared.tsv  RNA-view weights of the main RNA factor in both cohorts, shared genes
#   mofa_main_factor_by_view.tsv, mofa_factor_correlation.tsv
#
# One MOFA model per cohort, unsupervised (no labels or histology are used to fit it). Views:
#   RNA          MOFA_N_GENES most variable genes (log2 expression), centred
#   Methylation  MOFA_N_CPGS most variable CpGs (M-values); in GSE267015 the array platform is
#                regressed out first (limma::removeBatchEffect)
#   CNV          39 arm-level log2 ratios + MYCN and RB1 focal values (step 20)
# Tumours missing a view are kept (MOFA handles missing views). The main factor of each cohort is
# the factor most correlated with the expression phenotype score; it is oriented so that higher
# values mean a more aggressive (de-differentiated) tumour. MOFA runs its Python back end through
# basilisk, which sets up its own Python environment on first use.
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
if (length(INVASION_EXCLUDE) > 0) stop("Unset INVASION_EXCLUDE.")

load_pkg("MOFA2", bioc = TRUE)
load_pkg("limma", bioc = TRUE)
load_pkg("dplyr")
load_pkg("pROC", bioc = FALSE)

clip <- function(b) pmin(pmax(b, METH_BETA_CLIP), 1 - METH_BETA_CLIP)
to_M <- function(b) log2(clip(b) / (1 - clip(b)))
top_var <- function(M, n) { v <- apply(M, 1, stats::var, na.rm = TRUE); M[head(order(v, decreasing = TRUE), n), , drop = FALSE] }
centre <- function(M) M - rowMeans(M, na.rm = TRUE)
pad <- function(M, ids) {   # features x all tumours, NA columns for tumours missing the view
  out <- matrix(NA_real_, nrow(M), length(ids), dimnames = list(rownames(M), ids))
  k <- intersect(colnames(M), ids); out[, k] <- M[, k]; out
}

fit_mofa <- function(views, seed) {
  mo <- MOFA2::create_mofa(views)
  dopt <- MOFA2::get_default_data_options(mo); dopt$scale_views <- TRUE
  mopt <- MOFA2::get_default_model_options(mo); mopt$num_factors <- MOFA_N_FACTORS
  topt <- MOFA2::get_default_training_options(mo); topt$seed <- seed; topt$convergence_mode <- "medium"
  topt$verbose <- FALSE
  mo <- MOFA2::prepare_mofa(mo, data_options = dopt, model_options = mopt, training_options = topt)
  MOFA2::run_mofa(mo, outfile = file.path(tempdir(), paste0("mofa_", seed, "_", as.integer(Sys.time()), ".hdf5")),
                  use_basilisk = TRUE, save_data = FALSE)
}

# -----------------------------------------------------------------------------
# 1. Views for each cohort
# -----------------------------------------------------------------------------
cnv <- read.delim(file.path(DIR_CNV, "cnv_samples.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
cnv_cols <- c(grep("^[0-9]+[pq]$", names(cnv), value = TRUE), "MYCN", "RB1")
# copy-number view: QC-passed arrays only (step 20)
cnv_view <- function(coh) { d <- cnv[cnv$cohort == coh & cnv$qc_pass, ]; m <- t(as.matrix(d[, cnv_cols])); colnames(m) <- d$tumour; centre(m) }

# discovery (GSE267015 + GSE268136): tumour ids = methylation titles (RB_xx)
smp  <- read.delim(file.path(DIR_METH, "meth_samples.tsv"), stringsAsFactors = FALSE); smp <- smp[smp$used, ]
ps   <- read.delim(file.path(DIR_METH, "meth_phenotype_samples.tsv"), stringsAsFactors = FALSE)
rna  <- readRDS(file.path(DIR_BULK, "rnaseq_log2tpm_GSE268136.rds"))
lim  <- read.delim(file.path(DIR_BULK, "limma_clusterC_vs_AB_GSE268136.tsv"), stringsAsFactors = FALSE)
beta <- readRDS(file.path(DIR_METH, "meth_beta.rds"))
E_d <- rna[intersect(lim$gene, rownames(rna)), smp$rnaseq_sample[smp$paired_rna], drop = FALSE]
colnames(E_d) <- smp$title[smp$paired_rna]
M_d <- limma::removeBatchEffect(to_M(beta[, smp$title]), batch = smp$platform)
ids_d <- smp$title
views_d <- list(RNA = pad(centre(top_var(E_d, MOFA_N_GENES)), ids_d), Methylation = pad(centre(top_var(M_d, MOFA_N_CPGS)), ids_d),
                CNV = pad(cnv_view("GSE267015"), ids_d))

# held out (Liu): tumour ids = RBxx
X_l <- readRDS(file.path(DIR_VAL, "liu_expression.rds"))
B_l <- readRDS(file.path(DIR_VAL, "liu_beta_top20k.rds"))
lab <- read.delim(file.path(DIR_VAL, "liu_samples.tsv"), stringsAsFactors = FALSE)
ids_l <- sort(unique(c(colnames(X_l), colnames(B_l), cnv$tumour[cnv$cohort == "GSE58783"])))
views_l <- list(RNA = pad(centre(top_var(X_l, MOFA_N_GENES)), ids_l), Methylation = pad(centre(top_var(to_M(B_l), MOFA_N_CPGS)), ids_l),
                CNV = pad(cnv_view("GSE58783"), ids_l))
for (nm in c("d", "l")) {
  v <- get(paste0("views_", nm))
  message(c(d = "GSE267015/GSE268136", l = "Liu")[nm], ": ", length(get(paste0("ids_", nm))), " tumours; ",
          paste(names(v), vapply(v, function(m) sum(colSums(!is.na(m)) > 0), 0), "with data", collapse = ", "))
}

# -----------------------------------------------------------------------------
# 2. Fit, pick and orient the main factor
# -----------------------------------------------------------------------------
pheno_d <- setNames(ps$score, ps$title)
mv <- read.delim(file.path(DIR_VAL, "liu_marker_values.tsv"), stringsAsFactors = FALSE)
mv <- mv[mv$marker == "phenotype score", ]; pheno_l <- setNames(mv$value, mv$tumour)

analyse <- function(views, pheno, cohort) {
  m <- fit_mofa(views, GLOBAL_SEED)
  Z <- MOFA2::get_factors(m)[[1]]
  r2 <- MOFA2::get_variance_explained(m)$r2_per_factor[[1]]
  keep <- rownames(r2)[apply(r2, 1, max) >= 2]                      # factors explaining >= 2% of any view
  k <- intersect(names(pheno), rownames(Z))
  rho <- vapply(colnames(Z), function(f) stats::cor(Z[k, f], pheno[k], method = "spearman"), 0)
  main <- names(which.max(abs(rho[keep])))
  # every factor oriented so that higher = more aggressive (positive correlation with the phenotype score)
  Z <- sweep(Z, 2, sign(rho), `*`)
  # main factor per view: among factors explaining >= 0.5% of that view, the one most correlated with the
  # phenotype score (the axis can be carried by different, correlated factors in different views)
  main_view <- vapply(colnames(r2), function(v) { f <- rownames(r2)[r2[, v] >= 0.5]; f[which.max(abs(rho[f]))] }, "")
  W <- lapply(setNames(names(main_view), names(main_view)), function(v)
    sign(rho[[main_view[[v]]]]) * MOFA2::get_weights(m, views = v, factors = main_view[[v]])[[1]][, 1])
  list(model = m, Z = Z[, keep, drop = FALSE], r2 = r2, main = main, main_view = main_view, W = W, rho = rho, cohort = cohort)
}
res <- list(d = analyse(views_d, pheno_d, "GSE267015"), l = analyse(views_l, pheno_l, "GSE58783"))

ve <- dplyr::bind_rows(lapply(res, function(r) data.frame(cohort = r$cohort, factor = rownames(r$r2), r$r2,
                                                          main = rownames(r$r2) == r$main, check.names = FALSE)))
save_tsv(ve, file.path(DIR_MOFA, "mofa_variance_explained.tsv")); print(ve)
fz <- dplyr::bind_rows(lapply(res, function(r) data.frame(cohort = r$cohort, tumour = rownames(r$Z), r$Z, check.names = FALSE)))
fz$main_factor <- mapply(function(co, i) { r <- res[[c(GSE267015 = "d", GSE58783 = "l")[[co]]]]; r$Z[i, r$main] }, fz$cohort, fz$tumour)
save_tsv(fz, file.path(DIR_MOFA, "mofa_factors.tsv"))
mv_tab <- dplyr::bind_rows(lapply(res, function(r) data.frame(cohort = r$cohort, view = names(r$main_view), factor = unname(r$main_view),
  r2_view = mapply(function(f, v) r$r2[f, v], r$main_view, names(r$main_view)), rho_phenotype = unname(abs(r$rho[r$main_view])))))
save_tsv(mv_tab, file.path(DIR_MOFA, "mofa_main_factor_by_view.tsv")); print(mv_tab)
fc <- dplyr::bind_rows(lapply(res, function(r) { C <- stats::cor(r$Z, method = "spearman")
  data.frame(cohort = r$cohort, factor1 = rownames(C)[row(C)], factor2 = colnames(C)[col(C)], rho = c(C))[c(upper.tri(C)), ] }))
save_tsv(fc, file.path(DIR_MOFA, "mofa_factor_correlation.tsv"))

# -----------------------------------------------------------------------------
# 3. Factor associations
# -----------------------------------------------------------------------------
auc_row <- function(x, y, label, cohort, factor) {
  ok <- !is.na(x) & !is.na(y); if (sum(y[ok]) < 3 || sum(1 - y[ok]) < 3) return(NULL)
  r <- pROC::roc(y[ok], x[ok], levels = c(0, 1), direction = "<", quiet = TRUE); ci <- auc_ci(r)
  data.frame(cohort = cohort, factor = factor, test = label, n = sum(ok), n_pos = sum(y[ok]), estimate = ci[["auc"]],
             lo = ci[["lo"]], hi = ci[["hi"]], p = stats::wilcox.test(x[ok][y[ok] == 1], x[ok][y[ok] == 0])$p.value)
}
rho_row <- function(x, y, label, cohort, factor) {
  ok <- !is.na(x) & !is.na(y); ct <- stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  data.frame(cohort = cohort, factor = factor, test = label, n = sum(ok), n_pos = NA, estimate = unname(ct$estimate),
             lo = NA, hi = NA, p = ct$p.value)
}
cnv_d <- cnv[cnv$cohort == "GSE267015" & cnv$qc_pass, ]; cnv_l <- cnv[cnv$cohort == "GSE58783" & cnv$qc_pass, ]
tests <- dplyr::bind_rows(
  lapply(colnames(res$d$Z), function(f) { z <- res$d$Z[ids_d, f]; s <- smp[match(ids_d, smp$title), ]
    dplyr::bind_rows(rho_row(z, pheno_d[ids_d], "Spearman vs phenotype score", "GSE267015", f),
                     rho_row(z, cnv_d$fga[match(ids_d, cnv_d$tumour)], "Spearman vs fraction of genome altered", "GSE267015", f),
                     auc_row(z, as.integer(s$cluster_ryl_et_al %in% c("B", "C")), "AUC, Ryl cluster B/C vs A", "GSE267015", f),
                     auc_row(z, as.integer(s$platform == "EPIC"), "AUC, EPIC vs 450K (residual batch)", "GSE267015", f)) }),
  lapply(colnames(res$l$Z), function(f) { z <- res$l$Z[ids_l, f]; d <- lab[match(ids_l, lab$tumour), ]
    dplyr::bind_rows(rho_row(z, pheno_l[ids_l], "Spearman vs phenotype score", "GSE58783", f),
                     rho_row(z, cnv_l$fga[match(ids_l, cnv_l$tumour)], "Spearman vs fraction of genome altered", "GSE58783", f),
                     auc_row(z, as.integer(d$subtype == 2), "AUC, Liu subtype 2", "GSE58783", f),
                     auc_row(z, as.integer(d$high_risk_broad), "AUC, high risk (broad)", "GSE58783", f),
                     auc_row(z, as.integer(d$high_risk_strict), "AUC, high risk (strict)", "GSE58783", f),
                     auc_row(z, as.integer(d$high_risk_broad_in_subtype2), "AUC, high risk (broad) within subtype 2", "GSE58783", f),
                     auc_row(z, as.integer(d$high_risk_broad_in_subtype1), "AUC, high risk (broad) within subtype 1", "GSE58783", f)) }))
tests <- tests |> dplyr::group_by(cohort, test) |> dplyr::mutate(padj = p.adjust(p, "BH")) |> dplyr::ungroup()
tests$main <- tests$factor == ifelse(tests$cohort == "GSE267015", res$d$main, res$l$main)
save_tsv(tests, file.path(DIR_MOFA, "mofa_factor_tests.tsv")); print(as.data.frame(tests[tests$main, ]))

# -----------------------------------------------------------------------------
# 4. Weights of the main factors
# -----------------------------------------------------------------------------
tw <- dplyr::bind_rows(lapply(res, function(r) dplyr::bind_rows(lapply(names(r$W), function(v) {
  w <- r$W[[v]]; o <- head(order(-abs(w)), 25)
  data.frame(cohort = r$cohort, factor = r$main_view[[v]], view = v, feature = sub(paste0("_", v, "$"), "", names(w)[o]), weight = w[o])
}))))
save_tsv(tw, file.path(DIR_MOFA, "mofa_top_weights.tsv"))
pri <- read.delim(file.path(DIR_PRI, "priority_all_candidates.tsv"), stringsAsFactors = FALSE)
genes <- pri$gene[pri$shortlist %in% TRUE | pri$core %in% TRUE]
genes <- unique(c(genes, "MYLK", "DNM3"))
sw <- dplyr::bind_rows(lapply(res, function(r) { w <- r$W$RNA; names(w) <- sub("_RNA$", "", names(w))
  pct <- 100 * rank(abs(w)) / length(w)
  data.frame(cohort = r$cohort, factor = r$main_view[["RNA"]], gene = genes, in_rna_view = genes %in% names(w), weight = unname(w[genes]),
             abs_weight_percentile = unname(pct[genes]), direction_meta = sign(pri$g_re[match(genes, pri$gene)])) }))
save_tsv(sw, file.path(DIR_MOFA, "mofa_shortlist_weights.tsv")); print(sw)

# cross-cohort: RNA-view gene weights of the two main factors, over genes in both views
wd <- res$d$W$RNA; names(wd) <- sub("_RNA$", "", names(wd)); wl <- res$l$W$RNA; names(wl) <- sub("_RNA$", "", names(wl))
g <- intersect(names(wd), names(wl))
cc <- if (length(g) >= 20) stats::cor.test(wd[g], wl[g], method = "spearman", exact = FALSE) else NULL
cross <- data.frame(shared_genes = length(g), rho = if (is.null(cc)) NA else unname(cc$estimate), p = if (is.null(cc)) NA else cc$p.value)
save_tsv(cross, file.path(DIR_MOFA, "mofa_cross_cohort.tsv")); print(cross)
save_tsv(data.frame(gene = g, weight_discovery = unname(wd[g]), weight_liu = unname(wl[g])),
         file.path(DIR_MOFA, "mofa_rna_weights_shared.tsv"))

# -----------------------------------------------------------------------------
# 5. Stability of the main factor across seeds
# -----------------------------------------------------------------------------
stab <- dplyr::bind_rows(lapply(c("d", "l"), function(nm) {
  r <- res[[nm]]; v <- get(paste0("views_", nm)); z0 <- r$Z[, r$main]
  dplyr::bind_rows(lapply(MOFA_SEEDS, function(s) {
    Z <- MOFA2::get_factors(fit_mofa(v, s))[[1]][names(z0), , drop = FALSE]
    data.frame(cohort = r$cohort, seed = s, max_abs_r = max(abs(stats::cor(Z, z0, use = "pairwise.complete.obs"))))
  }))
}))
save_tsv(stab, file.path(DIR_MOFA, "mofa_stability.tsv")); print(stab)
message("\n[21_mofa_integration.R] Complete.")
