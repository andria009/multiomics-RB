# =============================================================================
# Step 23 — Methylation sensitivity analyses requested in the second review round
# Input : results/08_methylation/ meth_beta.rds, meth_phenotype_samples.tsv, dmp_phenotype.tsv.gz,
#         meth_probe_annotation.rds, cpg_panel.tsv (steps 13, 14, 16)
#         data/raw/methylation/idat/ (IDATs of GSE267015, extracted by step 20)
# Output: results/08_methylation/
#   meth_r2_composition_samples.tsv  per tumour: EpiDISH fractions (epithelial-like/tumour, fibroblast,
#                                    immune), surrogate variables, RNA-based covariates of step 14
#   meth_r2_models.tsv               per model: DMP count, main DMPs reproduced, and the change in effect
#                                    size of the main DMPs (slope through the origin of adjusted on main
#                                    delta beta; median ratio of |delta beta|)
#   meth_r2_panel_cpg.tsv            the 10 panel CpGs under every model
#   meth_r2_epic_check.tsv           EPIC artefact check by chromosome group (see section 4)
# Models (all: M-value ~ phenotype score (per SD) + platform [+ covariates], limma, as step 14):
#   rna_composition   + ESTIMATE StromalScore, purity, Muller-glia score (step 14; from the RNA-seq)
#   epidish           + fibroblast and immune fractions estimated from the methylation data themselves
#                       (EpiDISH, RPC, centEpiFibIC.m reference; Teschendorff et al. 2017)
#   sva               + surrogate variables (reference-free; sva on the 50,000 most variable CpGs,
#                       number by num.sv "be", at most 10), protecting the phenotype score
#   noob_single_pipeline  beta values recomputed from the IDATs with one pipeline (minfi noob, each
#                       platform separately), same model as the main analysis
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
if (length(INVASION_EXCLUDE) > 0) stop("Unset INVASION_EXCLUDE.")

load_pkg("limma", bioc = TRUE)
load_pkg("dplyr")
load_pkg("minfi", bioc = TRUE)
load_pkg("EpiDISH", bioc = TRUE)
load_pkg("sva", bioc = TRUE)
for (p in c("IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylationEPICmanifest")) load_pkg(p, bioc = TRUE)

clip <- function(b) pmin(pmax(b, METH_BETA_CLIP), 1 - METH_BETA_CLIP)
to_M <- function(b) log2(clip(b) / (1 - clip(b)))
beta <- readRDS(file.path(DIR_METH, "meth_beta.rds"))
ann  <- readRDS(file.path(DIR_METH, "meth_probe_annotation.rds"))
ps   <- read.delim(file.path(DIR_METH, "meth_phenotype_samples.tsv"), stringsAsFactors = FALSE)
dmp  <- as.data.frame(data.table::fread(file.path(DIR_METH, "dmp_phenotype.tsv.gz")))
panel <- read.delim(file.path(DIR_METH, "cpg_panel.tsv"), stringsAsFactors = FALSE)
panel_cpg <- unique(unlist(strsplit(panel$cpgs[panel$in_panel %in% TRUE], "[;, ]+")))
b_ps <- beta[, ps$title, drop = FALSE]
main_dmp <- dmp$cpg[dmp$dmp]
stopifnot(length(main_dmp) > 100, length(panel_cpg) >= 5)

fit_cpg <- function(b, design, coef = "score_sd") {    # as step 14
  fM <- limma::eBayes(limma::lmFit(to_M(b), design)); fB <- limma::lmFit(b, design)
  tt <- limma::topTable(fM, coef = coef, number = Inf, sort.by = "none")
  data.frame(cpg = rownames(b), delta_beta = fB$coefficients[, coef], t = tt$t, p = tt$P.Value, padj = tt$adj.P.Val)
}
# change in effect size of the main DMPs: slope through the origin (adjusted on main) and median ratio
compare <- function(r, model, n_tumours) {
  i <- match(main_dmp, r$cpg); x <- dmp$delta_beta[match(main_dmp, dmp$cpg)]; y <- r$delta_beta[i]; ok <- !is.na(y)
  is_dmp <- r$padj < METH_DMP_ADJP & abs(r$delta_beta) >= METH_DMP_DB_MIN
  data.frame(model = model, n_tumours = n_tumours, n_cpg = nrow(r), n_dmp = sum(is_dmp, na.rm = TRUE),
             main_dmp_tested = sum(ok),
             pct_main_dmp_reproduced = 100 * mean(sign(y[ok]) == sign(x[ok]) & r$p[i][ok] < 0.05),
             effect_slope_vs_main = sum(x[ok] * y[ok]) / sum(x[ok]^2),
             median_abs_ratio_vs_main = stats::median(abs(y[ok]) / abs(x[ok])),
             r_delta_beta_with_main = stats::cor(x[ok], y[ok]))
}
stopifnot(abs(sum(c(1, 2) * c(2, 4)) / sum(c(1, 2)^2) - 2) < 1e-12)   # slope helper: y = 2x gives 2

# -----------------------------------------------------------------------------
# 1. Composition estimated from the methylation data
# -----------------------------------------------------------------------------
data("centEpiFibIC.m", package = "EpiDISH", envir = environment())
ref <- centEpiFibIC.m[rownames(centEpiFibIC.m) %in% rownames(b_ps), , drop = FALSE]
message("EpiDISH reference CpGs present: ", nrow(ref), " of ", nrow(centEpiFibIC.m))
stopifnot(nrow(ref) >= 300)
frac <- EpiDISH::epidish(beta.m = b_ps, ref.m = ref, method = "RPC")$estF
ps$epi_fraction <- frac[ps$title, "Epi"]; ps$fib_fraction <- frac[ps$title, "Fib"]; ps$ic_fraction <- frac[ps$title, "IC"]

mod  <- model.matrix(~ score_sd + platform, data = ps); mod0 <- model.matrix(~ platform, data = ps)
Mall <- to_M(b_ps); v <- apply(Mall, 1, stats::var); Mtop <- Mall[head(order(v, decreasing = TRUE), 50000), ]
n_sv <- min(10L, sva::num.sv(Mtop, mod, method = "be"))
message("Surrogate variables: ", n_sv)
sv <- if (n_sv > 0) sva::sva(Mtop, mod, mod0, n.sv = n_sv)$sv else NULL
if (!is.null(sv)) { sv <- as.matrix(sv); colnames(sv) <- paste0("SV", seq_len(ncol(sv))); ps <- cbind(ps, sv) }

cs <- ps[, c("title", "platform", "cluster_ryl_et_al", "score", "stroma", "purity", "muller", "epi_fraction", "fib_fraction", "ic_fraction",
             grep("^SV[0-9]+$", names(ps), value = TRUE))]
save_tsv(cs, file.path(DIR_METH, "meth_r2_composition_samples.tsv"))
cv <- c("stroma", "purity", "muller", "epi_fraction", "fib_fraction", "ic_fraction", grep("^SV[0-9]+$", names(ps), value = TRUE))
print(data.frame(covariate = cv, rho_with_phenotype_score = vapply(cv, function(k) stats::cor(ps[[k]], ps$score, method = "spearman"), 0)), digits = 2)

# -----------------------------------------------------------------------------
# 2. Models on the deposited beta values
# -----------------------------------------------------------------------------
svf <- if (is.null(sv)) "" else paste("+", paste(colnames(sv), collapse = " + "))
forms <- list(main = ~ score_sd + platform,
              rna_composition = ~ score_sd + platform + scale(stroma) + scale(purity) + scale(muller),
              epidish = ~ score_sd + platform + fib_fraction + ic_fraction)
if (!is.null(sv)) forms$sva <- stats::as.formula(paste("~ score_sd + platform", svf))
fits <- lapply(forms, function(f) fit_cpg(b_ps, model.matrix(f, data = ps)))
models <- dplyr::bind_rows(lapply(names(fits), function(nm) compare(fits[[nm]], nm, nrow(ps))))

# -----------------------------------------------------------------------------
# 3. Single-pipeline beta values from the IDATs (noob, per platform)
# -----------------------------------------------------------------------------
idat_dir <- file.path(DIR_RAW_METH, "idat")
f <- list.files(idat_dir, "_Grn\\.idat(\\.gz)?$", full.names = TRUE)
if (!length(f)) stop("No IDATs in ", idat_dir, " (run step 20 first, which extracts GSE267015_RAW.tar)")
idat <- data.frame(basename = sub("_Grn\\.idat(\\.gz)?$", "", f), gsm = sub("^(GSM[0-9]+)_.*$", "\\1", basename(f)))
ps$basename <- idat$basename[match(ps$gsm, idat$gsm)]
stopifnot(!anyNA(ps$basename))
bn <- lapply(split(ps, ps$platform), function(d) {
  rg <- minfi::read.metharray(d$basename, force = TRUE)
  b <- minfi::getBeta(minfi::preprocessNoob(rg)); colnames(b) <- d$title; b })
common <- Reduce(intersect, c(lapply(bn, rownames), list(rownames(b_ps))))
b_noob <- do.call(cbind, lapply(bn, function(b) b[common, , drop = FALSE]))[, ps$title]
b_noob <- b_noob[stats::complete.cases(b_noob), , drop = FALSE]
message("noob beta: ", nrow(b_noob), " CpGs (of ", nrow(b_ps), " in the main analysis)")
r_noob <- fit_cpg(b_noob, model.matrix(~ score_sd + platform, data = ps))
fits$noob_single_pipeline <- r_noob
models <- dplyr::bind_rows(models, compare(r_noob, "noob_single_pipeline", nrow(ps)))
k <- intersect(rownames(b_noob), rownames(b_ps))
models$r_beta_with_deposited <- ifelse(models$model == "noob_single_pipeline",
  stats::median(vapply(ps$title, function(s) stats::cor(b_noob[k, s], b_ps[k, s]), 0)), NA)
save_tsv(models, file.path(DIR_METH, "meth_r2_models.tsv")); print(models, digits = 3)

pc <- dplyr::bind_rows(lapply(names(fits), function(nm) { r <- fits[[nm]]; i <- match(panel_cpg, r$cpg)
  data.frame(model = nm, cpg = panel_cpg, delta_beta = r$delta_beta[i], padj = r$padj[i]) }))
save_tsv(pc, file.path(DIR_METH, "meth_r2_panel_cpg.tsv"))

# -----------------------------------------------------------------------------
# 4. EPIC artefact check. The EPIC IDATs gave apparent losses of chromosomes 3-8 and gains of 17-22 in
#    copy number (step 20). If a related artefact reached the beta values, the CpGs on those chromosomes
#    would (a) be over-represented among DMPs, (b) show more platform x phenotype interaction, or
#    (c) agree less between the 450K-only and EPIC-only fits.
# -----------------------------------------------------------------------------
grp <- ifelse(ann$chr[match(rownames(b_ps), ann$cpg)] %in% METH_EPIC_ARTEFACT_CHR, "chr3-8,17-22", "other")
fi <- function(hit, label) {
  t2 <- table(factor(grp, c("chr3-8,17-22", "other")), factor(hit, c(TRUE, FALSE)))
  ft <- stats::fisher.test(t2)
  data.frame(check = label, pct_of_tested_on_artefact_chr = 100 * mean(grp == "chr3-8,17-22"),
             n_hits = sum(hit), pct_of_hits_on_artefact_chr = 100 * t2[1, 1] / max(1, sum(t2[, 1])),
             odds_ratio = unname(ft$estimate), p = ft$p.value)
}
r_main <- fits$main[match(rownames(b_ps), fits$main$cpg), ]
is_dmp <- r_main$padj < METH_DMP_ADJP & abs(r_main$delta_beta) >= METH_DMP_DB_MIN
fint <- limma::eBayes(limma::lmFit(to_M(b_ps), model.matrix(~ score_sd * platform, data = ps)))
icoef <- grep("score_sd:platform", colnames(fint$coefficients), value = TRUE)
tint <- limma::topTable(fint, coef = icoef, number = Inf, sort.by = "none")
p450 <- ps$platform == "450K"; pE <- ps$platform == "EPIC"
sd_in <- function(x) as.numeric(scale(x))
t450 <- limma::eBayes(limma::lmFit(to_M(b_ps[, p450]), model.matrix(~ sd_in(ps$score[p450]))))$t[, 2]
tEPIC <- limma::eBayes(limma::lmFit(to_M(b_ps[, pE]), model.matrix(~ sd_in(ps$score[pE]))))$t[, 2]
agree <- function(g) stats::cor(t450[grp == g], tEPIC[grp == g])
epic <- dplyr::bind_rows(
  fi(is_dmp & r_main$delta_beta < 0, "DMPs, lower methylation in aggressive tumours"),
  fi(is_dmp & r_main$delta_beta > 0, "DMPs, higher methylation in aggressive tumours"),
  fi(tint$adj.P.Val < 0.05, "platform x phenotype interaction, BH < 0.05"),
  fi(tint$P.Value < 0.001, "platform x phenotype interaction, p < 0.001"))
epic <- dplyr::bind_rows(epic, data.frame(check = c("r of t, 450K-only vs EPIC-only fits: chr3-8,17-22", "r of t, 450K-only vs EPIC-only fits: other chromosomes"),
                                          r = c(agree("chr3-8,17-22"), agree("other"))))
save_tsv(epic, file.path(DIR_METH, "meth_r2_epic_check.tsv")); print(epic, digits = 3)

message("\n[23_methylation_sensitivity_r2.R] Complete.")
