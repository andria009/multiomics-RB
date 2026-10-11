# =============================================================================
# Step 16 — A small CpG panel marking the aggressive (de-differentiated) state,
#           for methylation assays on FFPE DNA (e.g. pyrosequencing)
# Input : results/08_methylation/ (steps 13-14)
# Output: results/08_methylation/
#   cpg_regions.tsv          candidate regions: CpGs, genes, delta beta per SD (95% CI), lower bound
#   cpg_panel.tsv            the panel: top CPG_PANEL_N regions per direction, with bootstrap probability
#   cpg_panel_bootstrap.tsv  every candidate region's probability of entering the panel
#   cpg_panel_cv.tsv         cross-validated performance (selection repeated inside each training fold)
#   cpg_panel_scores.tsv     panel score of every tumour (paired and methylation-only)
#   figures/16_cpg_panel
#
# Rules (fixed before the results were seen):
#   candidate CpGs  DMPs from step 14 (BH < 0.05, |delta beta| >= METH_DMP_DB_MIN per SD) with a
#                   beta range (10th-90th percentile) >= CPG_MIN_RANGE, so they are measurable
#   regions         neighbouring candidate CpGs of the same direction within CPG_MERGE_BP bp
#                   and correlated (r >= CPG_MERGE_R) are merged; region beta = mean of its CpGs
#   ranking         lower 95% bound of |delta beta per SD| (region beta ~ score + platform)
#   panel           top CPG_PANEL_N regions per direction; panel score = mean beta of the
#                   hyper-methylated regions minus mean beta of the hypo-methylated regions
#   stability       N_BOOT resamples of tumours (within platform); ranking recomputed
#   performance     5-fold cross-validation repeated 10 times with DMP calling, filtering,
#                   merging and ranking redone inside each training fold
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
if (length(INVASION_EXCLUDE) > 0) stop("Unset INVASION_EXCLUDE.")

load_pkg("limma", bioc = TRUE)
load_pkg("dplyr")
load_pkg("ggplot2")
set.seed(GLOBAL_SEED)

beta <- readRDS(file.path(DIR_METH, "meth_beta.rds"))
ann  <- readRDS(file.path(DIR_METH, "meth_probe_annotation.rds"))
smp  <- read.delim(file.path(DIR_METH, "meth_samples.tsv"), stringsAsFactors = FALSE)
ps   <- read.delim(file.path(DIR_METH, "meth_phenotype_samples.tsv"), stringsAsFactors = FALSE)
dmp  <- data.table::fread(file.path(DIR_METH, "dmp_phenotype.tsv.gz"), data.table = FALSE)
clip <- function(b) pmin(pmax(b, METH_BETA_CLIP), 1 - METH_BETA_CLIP)
to_M <- function(b) log2(clip(b) / (1 - clip(b)))
B <- beta[, ps$title, drop = FALSE]

# slope of each row of Y on the score, adjusted for platform (vectorised OLS)
slopes <- function(Y, score, platform) {
  X <- cbind(1, score, as.integer(platform == "EPIC"))
  if (length(unique(platform)) < 2) X <- X[, 1:2]
  XtXi <- solve(crossprod(X)); Bh <- Y %*% X %*% XtXi
  res <- Y - Bh %*% t(X); s2 <- rowSums(res^2) / (ncol(Y) - ncol(X))
  data.frame(slope = Bh[, 2], se = sqrt(s2 * XtXi[2, 2]))
}

# candidate CpGs -> regions, given a DMP table and the tumours used
make_regions <- function(d, Bm) {
  q <- apply(Bm[d$cpg, , drop = FALSE], 1, stats::quantile, c(0.1, 0.9))
  d <- d[(q[2, ] - q[1, ]) >= CPG_MIN_RANGE, ]
  if (!nrow(d)) return(NULL)
  d <- d[order(d$chr, d$pos), ]
  reg <- integer(nrow(d)); r <- 0L
  for (i in seq_len(nrow(d))) {
    new <- i == 1 || d$chr[i] != d$chr[i - 1] || d$pos[i] - d$pos[i - 1] > CPG_MERGE_BP ||
      sign(d$delta_beta[i]) != sign(d$delta_beta[i - 1]) ||
      stats::cor(Bm[d$cpg[i], ], Bm[d$cpg[i - 1], ]) < CPG_MERGE_R
    if (new) r <- r + 1L
    reg[i] <- r
  }
  d$region <- paste0("R", reg)
  d
}
region_matrix <- function(d, Bm) {
  idx <- split(d$cpg, d$region)
  t(vapply(idx, function(cg) colMeans(Bm[cg, , drop = FALSE]), numeric(ncol(Bm))))
}
rank_panel <- function(st, n = CPG_PANEL_N) {
  st$lower <- abs(st$slope) - 1.96 * st$se
  c(head(rownames(st)[st$slope > 0][order(-st$lower[st$slope > 0])], n),
    head(rownames(st)[st$slope < 0][order(-st$lower[st$slope < 0])], n))
}
panel_score <- function(R, panel, slope_sign) {
  hy <- intersect(panel[slope_sign[panel] > 0], rownames(R)); ho <- intersect(panel[slope_sign[panel] < 0], rownames(R))
  (if (length(hy)) colMeans(R[hy, , drop = FALSE]) else 0) - (if (length(ho)) colMeans(R[ho, , drop = FALSE]) else 0)
}

# -----------------------------------------------------------------------------
# 1. Candidate regions and the panel (all paired tumours)
# -----------------------------------------------------------------------------
cand <- make_regions(dmp[dmp$dmp, c("cpg", "chr", "pos", "delta_beta", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island")], B)
if (is.null(cand)) stop("No candidate CpGs pass the range filter.")
R <- region_matrix(cand, B)
st <- slopes(R, ps$score_sd, ps$platform); rownames(st) <- rownames(R)
st$lower <- abs(st$slope) - 1.96 * st$se
regions <- cand |> dplyr::group_by(region) |>
  dplyr::summarise(chr = chr[1], start = min(pos), end = max(pos), n_cpg = dplyr::n(), cpgs = paste(cpg, collapse = ";"),
                   genes = paste(unique(unlist(strsplit(UCSC_RefGene_Name, ";"))), collapse = ";"),
                   context = paste(unique(Relation_to_Island), collapse = ";"), .groups = "drop") |>
  dplyr::mutate(delta_beta_per_sd = st[region, "slope"], ci_lo = st[region, "slope"] - 1.96 * st[region, "se"],
                ci_hi = st[region, "slope"] + 1.96 * st[region, "se"], lower_bound = st[region, "lower"],
                mean_beta = rowMeans(R[region, , drop = FALSE]))
PANEL <- rank_panel(st)
message("Candidate CpGs: ", nrow(cand), " in ", nrow(R), " regions (", sum(st$slope > 0), " hyper, ",
        sum(st$slope < 0), " hypo in aggressive tumours)")

# -----------------------------------------------------------------------------
# 2. Bootstrap stability of the panel (regions fixed, tumours resampled within platform)
# -----------------------------------------------------------------------------
hits <- setNames(numeric(nrow(R)), rownames(R))
for (b in seq_len(N_BOOT)) {
  idx <- unlist(lapply(split(seq_len(ncol(R)), ps$platform), function(i) i[sample.int(length(i), replace = TRUE)]))
  sb <- slopes(R[, idx, drop = FALSE], ps$score_sd[idx], ps$platform[idx]); rownames(sb) <- rownames(R)
  pb <- rank_panel(sb); hits[pb] <- hits[pb] + 1
}
regions$p_panel <- unname(hits[regions$region] / N_BOOT)
regions$in_panel <- regions$region %in% PANEL
regions <- regions[order(-regions$in_panel, -regions$lower_bound), ]
save_tsv(regions, file.path(DIR_METH, "cpg_regions.tsv"))
panel <- regions[regions$in_panel, ]
save_tsv(panel, file.path(DIR_METH, "cpg_panel.tsv"))
save_tsv(regions[, c("region", "genes", "delta_beta_per_sd", "lower_bound", "p_panel")],
         file.path(DIR_METH, "cpg_panel_bootstrap.tsv"))
print(panel[, c("region", "chr", "start", "end", "n_cpg", "genes", "delta_beta_per_sd", "ci_lo", "ci_hi", "p_panel")])

# -----------------------------------------------------------------------------
# 3. Cross-validated performance (whole selection repeated inside training folds)
# -----------------------------------------------------------------------------
cv <- dplyr::bind_rows(lapply(1:10, function(rep) {
  fold <- sample(rep_len(1:5, ncol(B)))
  pred <- rep(NA_real_, ncol(B))
  for (k in 1:5) {
    tr <- fold != k; te <- fold == k
    ptr <- ps[tr, ]; ptr$score_sd <- as.numeric(scale(ptr$score))
    des <- model.matrix(~ score_sd + platform, data = ptr)
    fM <- limma::eBayes(limma::lmFit(to_M(B[, tr]), des))
    tt <- limma::topTable(fM, coef = "score_sd", number = Inf, sort.by = "none")
    db <- limma::lmFit(B[, tr], des)$coefficients[, "score_sd"]
    d <- data.frame(cpg = rownames(B), delta_beta = db, padj = tt$adj.P.Val)
    d <- d[d$padj < METH_DMP_ADJP & abs(d$delta_beta) >= METH_DMP_DB_MIN, ]
    d <- dplyr::left_join(d, ann[, c("cpg", "chr", "pos")], by = "cpg")
    cr <- make_regions(d, B[, tr, drop = FALSE]); if (is.null(cr)) next
    Rtr <- region_matrix(cr, B[, tr, drop = FALSE])
    s <- slopes(Rtr, ptr$score_sd, ptr$platform); rownames(s) <- rownames(Rtr)
    pn <- rank_panel(s)
    Rte <- region_matrix(cr[cr$region %in% pn, ], B[, te, drop = FALSE])
    pred[te] <- panel_score(Rte, pn, setNames(sign(s$slope), rownames(s)))
  }
  data.frame(repeat_id = rep, title = ps$title, platform = ps$platform, score = ps$score, cv_panel_score = pred)
}))
perf <- cv |> dplyr::group_by(repeat_id) |>
  dplyr::summarise(spearman_all = stats::cor(score, cv_panel_score, method = "spearman", use = "complete.obs"),
                   spearman_450K = stats::cor(score[platform == "450K"], cv_panel_score[platform == "450K"],
                                              method = "spearman", use = "complete.obs"),
                   spearman_EPIC = stats::cor(score[platform == "EPIC"], cv_panel_score[platform == "EPIC"],
                                              method = "spearman", use = "complete.obs"), .groups = "drop")
save_tsv(perf, file.path(DIR_METH, "cpg_panel_cv.tsv"))
message("Cross-validated Spearman rho, held-out panel score vs phenotype score: median ",
        round(stats::median(perf$spearman_all), 2), " (range ", paste(round(range(perf$spearman_all), 2), collapse = " to "), ")")

# -----------------------------------------------------------------------------
# 4. Panel score for every tumour (methylation-only tumours have no phenotype score)
# -----------------------------------------------------------------------------
sgn <- setNames(sign(st$slope), rownames(st))
Rall <- region_matrix(cand[cand$region %in% PANEL, ], beta)
sc_all <- data.frame(title = colnames(beta), panel_score = panel_score(Rall, PANEL, sgn)) |>
  dplyr::left_join(smp[, c("title", "platform", "cluster_ryl_et_al", "paired_rna")], by = "title") |>
  dplyr::left_join(ps[, c("title", "score")], by = "title")
save_tsv(sc_all, file.path(DIR_METH, "cpg_panel_scores.tsv"))
print(sc_all |> dplyr::group_by(cluster_ryl_et_al, has_score = !is.na(score)) |>
        dplyr::summarise(n = dplyr::n(), median_panel = stats::median(panel_score), .groups = "drop"))

p <- ggplot2::ggplot(sc_all[!is.na(sc_all$score), ], ggplot2::aes(score, panel_score, colour = cluster_ryl_et_al, shape = platform)) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::labs(x = "Expression-based phenotype score", y = "CpG panel score (in-sample)",
                title = "CpG panel score vs phenotype score (cross-validated performance in cpg_panel_cv.tsv)") + PLOT_THEME
save_plot(p, file.path(DIR_FIG, "16_cpg_panel"), width = 7.5, height = 5.5)
message("\n[16_cpg_marker_panel.R] Complete.")
