# =============================================================================
# Step 18 — Held-out validation in the Liu/Curie cohort, including histology
# Input : data/raw/liu/GSE58780_series_matrix.txt.gz            (expression, HG-U133 Plus 2.0,
#                                                                 Brainarray Entrez-gene CDF)
#         data/raw/liu/GSE58783_series_matrix.txt.gz            (450K sample annotation)
#         data/raw/liu/GSE58783_processed_data_i_methylation_retino_all_samples.csv.gz (betas)
#         data/supplementary/41467_2021_25792_MOESM4_ESM.xlsx   (Liu et al. 2021 Table S1:
#                                                                 histology, subtypes, 9 CpGs)
#         results/ from steps 04, 11, 13, 16, 17
# Output: results/10_liu_validation/
#   liu_samples.tsv              tumours, histology, high-risk labels, Liu subtype, data available
#   liu_expression_tests.tsv     genes and phenotype score vs high-risk histology and vs subtype
#   liu_methylation_tests.tsv    CpG panel, Liu 9-CpG predictor and promoter methylation vs the same
#   liu9_in_gse267015.tsv        Liu 9-CpG predictor applied to the discovery methylation cohort
#   liu_marker_values.tsv        per-tumour standardised marker values with histology and subtype
#   liu_power.tsv                power of the histology tests (Hanley-McNeil) for true AUCs 0.65-0.80
#   figures/18_liu_histology
#
# This cohort was not used for any selection (genes, filters, CpG panel). Definitions fixed
# before it was analysed: high risk (strict) = postlaminar optic nerve invasion, or deep
# choroidal invasion, or scleral invasion; broad = strict or "choroid_extended" (config.R).
# Tests: Wilcoxon and AUC (DeLong; not estimable under complete separation); logistic
# regression of high risk on the marker (per SD) with and without adjustment for Liu subtype.
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
if (length(INVASION_EXCLUDE) > 0) stop("Unset INVASION_EXCLUDE.")

load_pkg("GEOquery", bioc = TRUE)
load_pkg("readxl", bioc = FALSE)
load_pkg("dplyr")
load_pkg("ggplot2")
load_pkg("pROC", bioc = FALSE)
load_pkg("org.Hs.eg.db", bioc = TRUE)

f_expr <- file.path(DIR_RAW_LIU, "GSE58780_series_matrix.txt.gz")
f_mmat <- file.path(DIR_RAW_LIU, "GSE58783_series_matrix.txt.gz")
f_beta <- file.path(DIR_RAW_LIU, "GSE58783_processed_data_i_methylation_retino_all_samples.csv.gz")
f_supp <- here::here("data", "supplementary", LIU_SUPP_FILE)
for (f in c(f_expr, f_mmat, f_beta, f_supp)) if (!file.exists(f)) stop("Missing ", f, " (run data/download_raw_data.sh)")
avg_by <- function(M, g) sapply(split(seq_along(g), g), function(i) rowMeans(M[, i, drop = FALSE], na.rm = TRUE))
# "RB4_REP1 (methylation)" -> "RB4"; "RBsjd7 (expression)" -> "RBsjd7"; fetal retinas (RET...) are dropped later
tumour_id <- function(title) sub("_REP[0-9]+$", "", sub("\\s*\\(.*$", "", trimws(title)))
is_tumour <- function(id) grepl("^RB[A-Za-z]*[0-9]+$", id)

# -----------------------------------------------------------------------------
# 1. Histology and subtypes (Liu et al. Table S1)
# -----------------------------------------------------------------------------
pt <- readxl::read_excel(f_supp, sheet = "patient information", skip = 3, na = c("", "NA"))   # missing histology is coded "NA"
pt <- as.data.frame(pt[!is.na(pt$`sample ID`) & is_tumour(pt$`sample ID`), ])
# final subtype for all 102 tumours (omics or 9-CpG pyrosequencing); 3 = unclassified -> NA
om <- as.data.frame(readxl::read_excel(f_supp, sheet = "overview of subtype assign.", skip = 3))
om <- om[!is.na(om$`sample ID`), ]
om$subtype <- ifelse(om$`final subtype classification` %in% c(1, 2), om$`final subtype classification`, NA)
lab <- data.frame(tumour = pt$`sample ID`, optic_nerve = pt$`Optic nerve invasion`,
                  choroid_sclera = pt$`Choroid and sclera invasion`, chemotherapy = pt$Chemotherapy,
                  subtype = om$subtype[match(pt$`sample ID`, om$`sample ID`)])
known <- !is.na(lab$optic_nerve) | !is.na(lab$choroid_sclera)
hr <- function(ch_set) ifelse(!known, NA, (lab$optic_nerve %in% HR_STRICT_ON) | (lab$choroid_sclera %in% ch_set))
lab$high_risk_strict <- hr(HR_STRICT_CH)
lab$high_risk_broad  <- hr(HR_BROAD_CH)
lab$untreated <- lab$chemotherapy %in% "no"
# Revision: broad definition primary; sensitivity outcomes as copies of the labels set to NA outside
# the stratum (untreated eyes; within each subtype)
for (v in c("high_risk_broad", "high_risk_strict")) {
  lab[[paste0(v, "_untreated")]] <- ifelse(lab$untreated, lab[[v]], NA)
  for (k in 1:2) lab[[paste0(v, "_in_subtype", k)]] <- ifelse(lab$subtype %in% k, lab[[v]], NA)
}
HR_OUTCOMES <- c("high_risk_broad", "high_risk_strict", "high_risk_broad_untreated", "high_risk_strict_untreated",
                 "high_risk_broad_in_subtype2", "high_risk_broad_in_subtype1",
                 "high_risk_strict_in_subtype2", "high_risk_strict_in_subtype1")
message("Treated before enucleation: ", sum(!lab$untreated), " (high risk broad ",
        sum(!lab$untreated & lab$high_risk_broad %in% TRUE), ", strict ", sum(!lab$untreated & lab$high_risk_strict %in% TRUE), ")")
message("Liu tumours with histology: ", sum(known), "; high risk strict ", sum(lab$high_risk_strict, na.rm = TRUE),
        ", broad ", sum(lab$high_risk_broad, na.rm = TRUE))
print(table(strict = lab$high_risk_strict, subtype = lab$subtype, useNA = "ifany"))
# Do high-risk features differ between Liu subtypes? (Liu et al. reported no significant difference)
sh <- dplyr::bind_rows(lapply(c("high_risk_strict", "high_risk_broad"), function(v) {
  tb <- table(subtype = lab$subtype, hr = lab[[v]]); ft <- stats::fisher.test(tb)
  data.frame(definition = v, subtype1_hr = tb["1", "TRUE"], subtype1_n = sum(tb["1", ]),
             subtype2_hr = tb["2", "TRUE"], subtype2_n = sum(tb["2", ]), odds_ratio = unname(ft$estimate),
             or_lo = ft$conf.int[1], or_hi = ft$conf.int[2], p = ft$p.value)
}))
save_tsv(sh, file.path(DIR_VAL, "liu_subtype_vs_histology.tsv")); print(sh)

# -----------------------------------------------------------------------------
# 2. Expression (GSE58780): replicate arrays averaged, Entrez IDs -> symbols
# -----------------------------------------------------------------------------
ge <- GEOquery::getGEO(filename = f_expr, getGPL = FALSE)
X <- Biobase::exprs(ge); pe <- Biobase::pData(ge)
if (max(X, na.rm = TRUE) > 50) X <- log2(X + 1)
tid <- tumour_id(pe$title); ok <- is_tumour(tid)
X <- X[, ok, drop = FALSE]; tid <- tid[ok]
X <- avg_by(X, tid)
sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, sub("_at$", "", rownames(X)),
                                              column = "SYMBOL", keytype = "ENTREZID"))
keep <- !is.na(sym) & !duplicated(sym)
X <- X[keep, , drop = FALSE]; rownames(X) <- sym[keep]
message("GSE58780: ", nrow(X), " genes x ", ncol(X), " tumours")

# -----------------------------------------------------------------------------
# 3. Methylation (GSE58783): processed betas
# The file's columns "SAMPLE k" follow the order of the samples in the GEO series (k-th sample),
# NOT the "geo dataset serie sample id" characteristic (which uses another numbering). The
# mapping is verified below against the Illumina beta values Liu et al. list for their 9 CpGs.
# -----------------------------------------------------------------------------
pm <- read_series_pheno(f_mmat)
map <- data.frame(sample = paste("SAMPLE", seq_len(nrow(pm))), tumour = tumour_id(pm$title))
hdr <- names(data.table::fread(f_beta, nrows = 0))
message("GSE58783 processed file: ", length(hdr), " columns; first: ", paste(head(hdr, 4), collapse = " | "))
# columns: "SAMPLE n_AVG_Beta" and "SAMPLE n_Detection_PValue"; betas with detection p > 0.01 set to NA
bcols <- grep("_AVG_Beta$", hdr, value = TRUE); pcols <- sub("_AVG_Beta$", "_Detection_PValue", bcols)
stopifnot(length(bcols) > 0, all(pcols %in% hdr))
bt <- data.table::fread(f_beta, select = c(hdr[1], bcols, pcols), data.table = FALSE)
Bv <- as.matrix(bt[, bcols]); Pv <- as.matrix(bt[, pcols]); rownames(Bv) <- bt[[1]]; rm(bt)
Bv[Pv > 0.01] <- NA; rm(Pv)
message("  probes read: ", nrow(Bv), " (a complete 450K file has about 485,000; a much smaller number means a truncated download)")
stopifnot(nrow(Bv) > 400000)
stopifnot(setequal(sub("_AVG_Beta$", "", colnames(Bv)), map$sample))
btid <- map$tumour[match(sub("_AVG_Beta$", "", colnames(Bv)), map$sample)]
# verification: Liu Table S1 Illumina betas of the 9 CpGs must match the mapped columns
x9 <- as.data.frame(readxl::read_excel(f_supp, sheet = "9-CpGs data and subtype assig.", col_names = FALSE))
r_cg <- which(apply(x9, 1, function(r) any(grepl("^cg[0-9]+$", r))))[1]      # row naming the CpGs
r_hd <- which(x9[[1]] %in% "sample ID")[1]                                     # header row; data follow
cg_at <- which(grepl("^cg[0-9]+$", unlist(x9[r_cg, ]))); cg9 <- unlist(x9[r_cg, cg_at])
v9 <- x9[-(1:r_hd), c(1, 2, cg_at + 1)]; names(v9) <- c("tumour", "series", cg9)   # Illumina value next to pyro
v9 <- v9[v9$series %in% "initial series", ]
v9[cg9] <- lapply(v9[cg9], as.numeric); v9 <- v9[stats::complete.cases(v9[cg9]), ]
d9 <- vapply(seq_len(nrow(v9)), function(i) {
  cols <- which(btid == v9$tumour[i]); if (!length(cols)) return(NA_real_)
  min(colSums((Bv[cg9, cols, drop = FALSE] - unlist(v9[i, cg9]))^2, na.rm = TRUE))
}, numeric(1))
message("  sample mapping check: ", sum(d9 < 0.01, na.rm = TRUE), " of ", nrow(v9),
        " tumours match Liu's published 9-CpG Illumina values")
stopifnot(all(d9 < 0.01, na.rm = TRUE), sum(!is.na(d9)) >= 10)
okb <- !is.na(btid) & is_tumour(btid)
Bv <- Bv[, okb, drop = FALSE]; btid <- btid[okb]
Bv <- avg_by(Bv, btid)
message("GSE58783: ", nrow(Bv), " CpGs x ", ncol(Bv), " tumours")

# matrices for the multi-omics factor analysis (step 21): expression and the most variable CpGs
save_rds(X, file.path(DIR_VAL, "liu_expression.rds"), "Liu expression (genes x tumours)")
# Round-2 revision: the MOFA view uses the same probes as the discovery cohort (autosomal, SNP- and
# cross-hybridisation-filtered; step 13). The first revision used all probes, so sex-chromosome CpGs
# entered the Liu methylation view; that version is kept for the comparison in step 21.
top20k <- function(B) { v <- apply(B, 1, stats::var, na.rm = TRUE); B[head(order(v, decreasing = TRUE), 20000), , drop = FALSE] }
disc_cpg <- rownames(readRDS(file.path(DIR_METH, "meth_beta.rds")))
save_rds(top20k(Bv[rownames(Bv) %in% disc_cpg, , drop = FALSE]), file.path(DIR_VAL, "liu_beta_top20k.rds"),
         "Liu methylation, 20,000 most variable CpGs (discovery probe set)")
save_rds(top20k(Bv), file.path(DIR_VAL, "liu_beta_top20k_allprobes.rds"),
         "Liu methylation, 20,000 most variable CpGs (all probes, as in revision 1)")
lab$has_expression <- lab$tumour %in% colnames(X)
lab$has_methylation <- lab$tumour %in% colnames(Bv)
save_tsv(lab, file.path(DIR_VAL, "liu_samples.tsv"))

# -----------------------------------------------------------------------------
# 4. Tests
# -----------------------------------------------------------------------------
VALUES <- list()   # per-tumour standardised marker values (for the within-subtype figure)
test_marker <- function(x, ids, name, layer) {
  d <- lab[match(ids, lab$tumour), ]; d$x <- as.numeric(scale(x))
  VALUES[[paste(layer, name)]] <<- data.frame(tumour = ids, marker = name, layer = layer, value = d$x)
  dplyr::bind_rows(lapply(c(HR_OUTCOMES, "subtype2"), function(outcome) {
    y <- if (outcome == "subtype2") as.integer(d$subtype == 2) else as.integer(d[[outcome]])
    ok <- !is.na(y) & !is.na(d$x)
    if (sum(y[ok]) < 3 || sum(1 - y[ok]) < 3) return(NULL)
    r <- pROC::roc(y[ok], d$x[ok], levels = c(0, 1), direction = "<", quiet = TRUE)
    ci <- auc_ci(r)
    out <- data.frame(marker = name, layer = layer, outcome = outcome, n = sum(ok), n_pos = sum(y[ok]),
                      auc = ci[["auc"]], auc_lo = ci[["lo"]], auc_hi = ci[["hi"]],
                      wilcox_p = stats::wilcox.test(d$x[ok][y[ok] == 1], d$x[ok][y[ok] == 0])$p.value)
    if (outcome != "subtype2") {
      fit <- function(f) tryCatch(summary(stats::glm(f, family = binomial, data = data.frame(y = y[ok], x = d$x[ok], s = factor(d$subtype[ok]))))$coefficients["x", ],
                                  error = function(e) rep(NA_real_, 4))
      u <- fit(y ~ x); a <- if (length(unique(d$subtype[ok])) > 1) fit(y ~ x + s) else rep(NA_real_, 4)
      out$or_per_sd <- exp(u[1]); out$or_p <- u[4]
      out$or_adj_subtype <- exp(a[1]); out$or_adj_lo <- exp(a[1] - 1.96 * a[2]); out$or_adj_hi <- exp(a[1] + 1.96 * a[2]); out$or_adj_p <- a[4]
    }
    out
  }))
}

pri  <- read.delim(file.path(DIR_PRI, "priority_all_candidates.tsv"), stringsAsFactors = FALSE)
mo   <- read.delim(file.path(DIR_MO, "multiomic_shortlist.tsv"), stringsAsFactors = FALSE)
cand <- read.delim(file.path(DIR_INV, "invasion_candidates.tsv"), stringsAsFactors = FALSE)
genes <- unique(c(pri$gene[pri$shortlist %in% TRUE | pri$core %in% TRUE], mo$gene))

# Expression: phenotype score (all 3,298 candidates; none selected with this cohort) and genes
Z <- zscore_rows(X); Z <- Z[stats::complete.cases(Z), , drop = FALSE]
up <- intersect(cand$gene[cand$g_re > 0], rownames(Z)); dn <- intersect(cand$gene[cand$g_re < 0], rownames(Z))
score <- score_updown(Z, up, dn)
expr_tests <- dplyr::bind_rows(
  test_marker(score, colnames(Z), "phenotype score", "expression"),
  dplyr::bind_rows(lapply(intersect(genes, rownames(Z)), function(g) {
    s <- if (pri$g_re[match(g, pri$gene)] < 0) -1 else 1    # oriented so that higher = more aggressive
    test_marker(s * Z[g, ], colnames(Z), g, "expression")
  })))
# BH within outcome and layer; the held-out subtype AUC is shown next to every histology result (an
# odds ratio adjusted for subtype is not interpretable when the marker nearly separates subtype)
add_cols <- function(tt) {
  tt <- tt |> dplyr::group_by(outcome) |> dplyr::mutate(padj = p.adjust(wilcox_p, "BH")) |> dplyr::ungroup()
  sa <- tt[tt$outcome == "subtype2", c("marker", "auc")]
  tt$subtype_auc <- sa$auc[match(tt$marker, sa$marker)]
  tt
}
expr_tests <- add_cols(expr_tests)
save_tsv(expr_tests, file.path(DIR_VAL, "liu_expression_tests.tsv"))
print(as.data.frame(expr_tests[expr_tests$marker == "phenotype score" | expr_tests$outcome == "high_risk_strict", ]))

# Methylation: CpG panel score (step 16), Liu 9-CpG predictor, promoter methylation of genes
panel <- read.delim(file.path(DIR_METH, "cpg_panel.tsv"), stringsAsFactors = FALSE)
reg_score <- function(B) {
  v <- vapply(seq_len(nrow(panel)), function(i) {
    cg <- intersect(strsplit(panel$cpgs[i], ";")[[1]], rownames(B))
    if (!length(cg)) return(rep(NA_real_, ncol(B)))
    colMeans(B[cg, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(B)))
  v <- matrix(v, nrow = ncol(B))
  hy <- panel$delta_beta_per_sd > 0
  rowMeans(v[, hy, drop = FALSE], na.rm = TRUE) - rowMeans(v[, !hy, drop = FALSE], na.rm = TRUE)
}
l9 <- as.data.frame(readxl::read_excel(f_supp, sheet = "9-CpGs characteristics", skip = 3))
l9 <- l9[grepl("^cg[0-9]+$", l9$`CpG position`), ]
liu9_score <- function(B) {    # distance to subtype-1 centroid minus distance to subtype-2 centroid
  cg <- intersect(l9$`CpG position`, rownames(B)); c1 <- l9$`centroid subtype 1`[match(cg, l9$`CpG position`)]
  c2 <- l9$`centroid subtype 2`[match(cg, l9$`CpG position`)]
  b <- B[cg, , drop = FALSE]
  sqrt(colSums((b - c1)^2, na.rm = TRUE)) - sqrt(colSums((b - c2)^2, na.rm = TRUE))
}
message("Liu 9 CpGs present: GSE58783 ", sum(l9$`CpG position` %in% rownames(Bv)), " of ", nrow(l9))
pg <- readRDS(file.path(DIR_METH, "meth_probe_gene.rds"))
meth_tests <- dplyr::bind_rows(
  test_marker(reg_score(Bv), colnames(Bv), "CpG panel score", "methylation"),
  test_marker(liu9_score(Bv), colnames(Bv), "Liu 9-CpG predictor", "methylation"),
  dplyr::bind_rows(lapply(intersect(genes, pg$gene[pg$promoter]), function(g) {
    cg <- intersect(pg$cpg[pg$promoter & pg$gene == g], rownames(Bv)); if (!length(cg)) return(NULL)
    s <- if (pri$g_re[match(g, pri$gene)] < 0) 1 else -1   # promoter methylation oriented: higher = more aggressive
    test_marker(s * colMeans(Bv[cg, , drop = FALSE], na.rm = TRUE), colnames(Bv), paste0(g, " promoter"), "methylation")
  })))
meth_tests <- add_cols(meth_tests)
save_tsv(meth_tests, file.path(DIR_VAL, "liu_methylation_tests.tsv"))
print(as.data.frame(meth_tests[meth_tests$layer == "methylation" & !grepl("promoter", meth_tests$marker), ]))

# Paired tumours: does the CpG panel track the expression phenotype score in an independent cohort?
both <- intersect(colnames(Z), colnames(Bv))
if (length(both) >= 10) message("Liu cohort, paired tumours (n = ", length(both), "): Spearman rho CpG panel vs expression phenotype score = ",
                                round(stats::cor(reg_score(Bv[, both, drop = FALSE]), score[both], method = "spearman"), 2))

# Liu 9-CpG predictor in the discovery methylation cohort (GSE267015)
b267 <- readRDS(file.path(DIR_METH, "meth_beta.rds"))
ps <- read.delim(file.path(DIR_METH, "meth_phenotype_samples.tsv"), stringsAsFactors = FALSE)
smp <- read.delim(file.path(DIR_METH, "meth_samples.tsv"), stringsAsFactors = FALSE)
l9_267 <- data.frame(title = colnames(b267), liu9 = liu9_score(b267)) |>
  dplyr::left_join(smp[, c("title", "platform", "cluster_ryl_et_al")], by = "title") |>
  dplyr::left_join(ps[, c("title", "score")], by = "title")
save_tsv(l9_267, file.path(DIR_VAL, "liu9_in_gse267015.tsv"))
message("GSE267015: Liu 9 CpGs present after filtering: ", sum(l9$`CpG position` %in% rownames(b267)), " of ", nrow(l9),
        "; Spearman rho with phenotype score = ", round(stats::cor(l9_267$liu9, l9_267$score, method = "spearman", use = "complete.obs"), 2))

# Figure: phenotype score and CpG panel score by histology
fd <- dplyr::bind_rows(
  data.frame(tumour = colnames(Z), value = as.numeric(scale(score)), marker = "Expression phenotype score"),
  data.frame(tumour = colnames(Bv), value = as.numeric(scale(reg_score(Bv))), marker = "CpG panel score")) |>
  dplyr::left_join(lab[, c("tumour", "high_risk_broad", "subtype")], by = "tumour") |>
  dplyr::filter(!is.na(high_risk_broad))
p <- ggplot2::ggplot(fd, ggplot2::aes(ifelse(high_risk_broad, "high risk", "other"), value, colour = factor(subtype))) +
  ggplot2::geom_boxplot(outlier.shape = NA, colour = "grey40") + ggplot2::geom_jitter(width = 0.15, size = 1.5) +
  ggplot2::facet_wrap(~ marker) + ggplot2::labs(x = NULL, y = "z-score", colour = "Liu subtype",
                                                title = "Liu/Curie cohort (held out): markers by histological risk (broad)") + PLOT_THEME
save_plot(p, file.path(DIR_FIG, "18_liu_histology"), width = 8, height = 5)
# -----------------------------------------------------------------------------
# 5. Per-tumour marker values and power of the held-out histology tests (revision)
# -----------------------------------------------------------------------------
vals <- dplyr::bind_rows(VALUES) |>
  dplyr::left_join(lab[, c("tumour", "subtype", "high_risk_broad", "high_risk_strict", "untreated")], by = "tumour")
save_tsv(vals, file.path(DIR_VAL, "liu_marker_values.tsv"))
# Hanley-McNeil SE of the AUC under the alternative; normal approximation; Bonferroni over the
# number of markers tested in that layer
hm_power <- function(A, n1, n0, alpha) {
  q1 <- A / (2 - A); q2 <- 2 * A^2 / (1 + A)
  se <- sqrt((A * (1 - A) + (n1 - 1) * (q1 - A^2) + (n0 - 1) * (q2 - A^2)) / (n1 * n0))
  stats::pnorm((A - 0.5) / se - stats::qnorm(1 - alpha / 2))
}
stopifnot(abs(hm_power(0.70, 12, 32, 0.05) - 0.56) < 0.01)   # reviewer's value for 12 vs 32
alltests <- dplyr::bind_rows(expr_tests, meth_tests)
pw <- alltests |> dplyr::filter(outcome %in% c("high_risk_broad", "high_risk_strict")) |>
  dplyr::group_by(layer, outcome) |>
  dplyr::summarise(n_markers = dplyr::n(), n_pos = max(n_pos), n = max(n), .groups = "drop") |>
  dplyr::mutate(n_neg = n - n_pos)
pw <- dplyr::bind_rows(lapply(HR_POWER_AUC, function(A) dplyr::mutate(pw, true_auc = A,
  power_alpha05 = hm_power(A, n_pos, n_neg, 0.05), power_bonferroni = hm_power(A, n_pos, n_neg, 0.05 / n_markers))))
save_tsv(pw, file.path(DIR_VAL, "liu_power.tsv")); print(as.data.frame(pw))

message("\n[18_liu_validation.R] Complete.")
