# =============================================================================
# Step 24 — Within-subtype pattern in the held-out cohort (round-2 revision; EXPLORATORY, post hoc)
# Input : results/10_liu_validation/liu_marker_values.tsv, liu_samples.tsv, liu_expression.rds (step 18)
#         results/11_cnv/cnv_samples.tsv (step 20), results/12_mofa/mofa_factors.tsv (step 21)
# Output: results/10_liu_validation/
#   liu_within_subtype_pattern.tsv  per marker: AUC within subtype 1 and within subtype 2 (broad, strict),
#                                   and whether they point in opposite directions
#   liu_boundary_tests.tsv          per marker: "intermediate position" tests: AUC of closeness to the
#                                   midpoint between the subtype medians, and the quadratic-term p of
#                                   logistic(high risk ~ z + z^2)
#   liu_s2_composition.tsv          within subtype 2: composition proxies vs high risk, and the inverse
#                                   marker associations adjusted for stromal and RPE/melanocyte proxies
# All marker values are standardised and oriented so that higher = more aggressive (step 18).
# These tests were suggested by the reviewer after the results were seen; they are hypothesis-generating.
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
source(here::here("code", "R", "gene_lists.R"))
if (length(INVASION_EXCLUDE) > 0) stop("Unset INVASION_EXCLUDE.")
load_pkg("dplyr")
load_pkg("pROC", bioc = FALSE)

lab  <- read.delim(file.path(DIR_VAL, "liu_samples.tsv"), stringsAsFactors = FALSE)
vals <- read.delim(file.path(DIR_VAL, "liu_marker_values.tsv"), stringsAsFactors = FALSE)
cnv  <- read.delim(file.path(DIR_CNV, "cnv_samples.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
cnv  <- cnv[cnv$cohort == "GSE58783" & cnv$qc_pass, ]
fz   <- read.delim(file.path(DIR_MOFA, "mofa_factors.tsv"), stringsAsFactors = FALSE)
fz   <- fz[fz$cohort == "GSE58783", ]
extra <- dplyr::bind_rows(
  data.frame(tumour = cnv$tumour, marker = "fraction of genome altered", layer = "copy number", value = as.numeric(scale(cnv$fga))),
  data.frame(tumour = fz$tumour, marker = "MOFA joint factor", layer = "joint", value = as.numeric(scale(fz$main_factor))))
extra <- dplyr::left_join(extra, lab[, c("tumour", "subtype", "high_risk_broad", "high_risk_strict", "untreated")], by = "tumour")
vals <- dplyr::bind_rows(vals, extra)
vals <- vals[!is.na(vals$subtype) & vals$untreated %in% TRUE, ]
OUT <- c(broad = "high_risk_broad", strict = "high_risk_strict")

auc_of <- function(x, y) {
  ok <- !is.na(x) & !is.na(y); if (sum(y[ok]) < 3 || sum(1 - y[ok]) < 3) return(c(auc = NA, lo = NA, hi = NA, n_pos = sum(y[ok]), n = sum(ok)))
  c(auc_ci(pROC::roc(y[ok], x[ok], levels = c(0, 1), direction = "<", quiet = TRUE)), n_pos = sum(y[ok]), n = sum(ok))
}
stopifnot(abs(auc_of(c(1, 2, 3, 4, 5, 6), c(0, 0, 0, 1, 1, 1))[["auc"]] - 1) < 1e-9)

# -----------------------------------------------------------------------------
# 1. Within-subtype AUCs side by side
# -----------------------------------------------------------------------------
pat <- dplyr::bind_rows(lapply(split(vals, paste(vals$layer, vals$marker)), function(d) dplyr::bind_rows(lapply(names(OUT), function(o) {
  y <- as.integer(d[[OUT[[o]]]])
  a1 <- auc_of(d$value[d$subtype == 1], y[d$subtype == 1]); a2 <- auc_of(d$value[d$subtype == 2], y[d$subtype == 2])
  data.frame(marker = d$marker[1], layer = d$layer[1], definition = o,
             auc_s1 = a1[["auc"]], lo_s1 = a1[["lo"]], hi_s1 = a1[["hi"]], n_pos_s1 = a1[["n_pos"]], n_s1 = a1[["n"]],
             auc_s2 = a2[["auc"]], lo_s2 = a2[["lo"]], hi_s2 = a2[["hi"]], n_pos_s2 = a2[["n_pos"]], n_s2 = a2[["n"]])
}))))
pat$opposite <- (pat$auc_s1 > 0.5 & pat$auc_s2 < 0.5) | (pat$auc_s1 < 0.5 & pat$auc_s2 > 0.5)
pat$s1_up_s2_down <- pat$auc_s1 > 0.5 & pat$auc_s2 < 0.5
save_tsv(pat, file.path(DIR_VAL, "liu_within_subtype_pattern.tsv"))
print(pat |> dplyr::group_by(definition) |>
        dplyr::summarise(markers = sum(!is.na(auc_s1) & !is.na(auc_s2)), s1_above_s2_below = sum(s1_up_s2_down, na.rm = TRUE),
                         s1_below_s2_above = sum(auc_s1 < 0.5 & auc_s2 > 0.5, na.rm = TRUE), .groups = "drop"))

# -----------------------------------------------------------------------------
# 2. Intermediate-position ("subtype boundary") tests
#    closeness = -|z - midpoint|, midpoint = mean of the two subtype medians of that marker;
#    AUC > 0.5 means high-risk tumours lie nearer the boundary. Quadratic logistic term p, and the
#    same with subtype as a covariate.
# -----------------------------------------------------------------------------
bnd <- dplyr::bind_rows(lapply(split(vals, paste(vals$layer, vals$marker)), function(d) dplyr::bind_rows(lapply(names(OUT), function(o) {
  mid <- mean(c(stats::median(d$value[d$subtype == 1], na.rm = TRUE), stats::median(d$value[d$subtype == 2], na.rm = TRUE)))
  y <- as.integer(d[[OUT[[o]]]]); ok <- !is.na(y) & !is.na(d$value)
  a <- auc_of(-abs(d$value - mid), y)
  qp <- function(f) tryCatch(summary(stats::glm(f, family = binomial, data = data.frame(y = y[ok], z = d$value[ok], s = factor(d$subtype[ok]))))$coefficients["I(z^2)", 4],
                             error = function(e) NA_real_)
  data.frame(marker = d$marker[1], layer = d$layer[1], definition = o, n = a[["n"]], n_pos = a[["n_pos"]],
             auc_closeness = a[["auc"]], lo = a[["lo"]], hi = a[["hi"]],
             p_quadratic = qp(y ~ z + I(z^2)), p_quadratic_adj_subtype = qp(y ~ z + I(z^2) + s))
}))))
bnd <- bnd |> dplyr::group_by(definition) |> dplyr::mutate(padj_quadratic = p.adjust(p_quadratic, "BH")) |> dplyr::ungroup()
save_tsv(bnd, file.path(DIR_VAL, "liu_boundary_tests.tsv"))
print(as.data.frame(bnd[bnd$marker %in% c("phenotype score", "TFF1", "EBF3", "SSTR2", "MOFA joint factor", "Liu 9-CpG predictor"), ]), digits = 2)

# -----------------------------------------------------------------------------
# 3. Within subtype 2: do high-risk eyes contain more non-tumour tissue? (composition dilution)
# -----------------------------------------------------------------------------
X <- readRDS(file.path(DIR_VAL, "liu_expression.rds"))
px <- proxy_scores(X, COMPOSITION_PROXIES[c("stroma", "endothelial", "immune", "rpe_melanocyte", "rod", "muller")])
s2 <- lab[lab$subtype %in% 2 & lab$untreated %in% TRUE & lab$tumour %in% px$sample, ]
px2 <- px[match(s2$tumour, px$sample), ]
comp <- dplyr::bind_rows(lapply(names(OUT), function(o) { y <- as.integer(s2[[OUT[[o]]]])
  dplyr::bind_rows(lapply(setdiff(names(px2), c("sample", "sex_score")), function(v) {
    a <- auc_of(px2[[v]], y); ok <- !is.na(y)
    data.frame(test = "proxy vs high risk", definition = o, variable = v, n = a[["n"]], n_pos = a[["n_pos"]], auc = a[["auc"]],
               lo = a[["lo"]], hi = a[["hi"]], p = tryCatch(stats::wilcox.test(px2[[v]][ok & y == 1], px2[[v]][ok & y == 0])$p.value, error = function(e) NA_real_))
  })) }))
mk <- vals[vals$subtype == 2 & vals$layer == "expression" & vals$marker %in% c("TFF1", "EBF3", "phenotype score"), ]
adj <- dplyr::bind_rows(lapply(names(OUT), function(o) dplyr::bind_rows(lapply(unique(mk$marker), function(g) {
  d <- mk[mk$marker == g, ]; d <- cbind(d, px[match(d$tumour, px$sample), c("stroma", "rpe_melanocyte")])
  d$y <- as.integer(d[[OUT[[o]]]]); d <- d[!is.na(d$y), ]
  co <- function(f) tryCatch(summary(stats::glm(f, family = binomial, data = d))$coefficients["value", ], error = function(e) rep(NA_real_, 4))
  u <- co(y ~ value); a <- co(y ~ value + stroma + rpe_melanocyte)
  data.frame(test = "marker OR per SD, unadjusted vs adjusted for stroma and RPE/melanocyte", definition = o, variable = g,
             n = nrow(d), n_pos = sum(d$y), or_unadj = exp(u[1]), p_unadj = u[4], or_adj = exp(a[1]), p_adj = a[4])
}))))
s2c <- dplyr::bind_rows(comp, adj)
save_tsv(s2c, file.path(DIR_VAL, "liu_s2_composition.tsv")); print(as.data.frame(s2c), digits = 2)

message("\n[24_within_subtype_exploratory.R] Complete.")
