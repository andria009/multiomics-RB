# =============================================================================
# Step 20 — Copy number from the methylation arrays (conumee)
# Input : data/raw/methylation/GSE267015_RAW.tar   (IDATs; discovery tumours, 450K and EPIC)
#         data/raw/liu/GSE58783_RAW.tar            (IDATs; held-out Liu cohort, 450K)
#         data/raw/liu/GSE58783_Raw_data_i_methylation_retino_all_samples.csv.gz (fallback if the
#                                                   RAW.tar holds no IDATs: GenomeStudio signals)
#         results/08_methylation/meth_samples.tsv, meth_phenotype_samples.tsv (steps 13-14)
#         results/10_liu_validation/liu_samples.tsv (step 18)
# Output: results/11_cnv/
#   cnv_samples.tsv            one row per tumour: cohort, platform, QC pass, fit noise, fraction of genome
#                              altered (FGA), MYCN and RB1 focal values, arm-level values
#   cnv_segments.tsv.gz        all segments (CBS) per array
#   cnv_arm_frequency.tsv      per cohort and arm: % tumours with gain / loss
#   cnv_tests_discovery.tsv    GSE267015: copy number vs phenotype score and Ryl cluster;
#                              MYCN call vs the submitters' MYCN annotation
#   cnv_tests_liu.tsv          Liu cohort: copy number vs subtype and high-risk histology
#   cache/                     per-array conumee results (re-used on later runs)
#
# Method: Meth + Unmeth intensities (minfi preprocessIllumina) of each array are normalised against
# copy-neutral 450K reference arrays (CopyNeutralIMA, 51 arrays; shared probes for EPIC arrays) by
# conumee CNV.fit, then binned and segmented (CBS). Arm values are length-weighted means of
# segment log2 ratios; gains / losses use CNV_ARM_THR. Replicate arrays of a tumour are averaged.
# The copy-number values are relative to tumour purity: low tumour content shrinks them.
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
if (length(INVASION_EXCLUDE) > 0) stop("Unset INVASION_EXCLUDE.")

load_pkg("minfi", bioc = TRUE)
load_pkg("conumee", bioc = TRUE)
load_pkg("CopyNeutralIMA", bioc = TRUE)
for (p in c("IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylationEPICmanifest",
            "IlluminaHumanMethylationEPICanno.ilm10b2.hg19")) load_pkg(p, bioc = TRUE)
load_pkg("dplyr")
load_pkg("pROC", bioc = FALSE)
set.seed(GLOBAL_SEED)

cache_dir <- file.path(DIR_CNV, "cache"); dir.create(cache_dir, showWarnings = FALSE)
data(exclude_regions, package = "conumee"); data(detail_regions, package = "conumee")

# hg19 centromere positions (Mb, UCSC cytoBand acen midpoints) used to split arms
CEN <- c(chr1 = 125.0, chr2 = 93.3, chr3 = 91.0, chr4 = 50.4, chr5 = 48.4, chr6 = 61.0, chr7 = 59.9, chr8 = 45.6,
         chr9 = 49.0, chr10 = 40.2, chr11 = 53.7, chr12 = 35.8, chr13 = 17.9, chr14 = 17.6, chr15 = 19.0,
         chr16 = 36.6, chr17 = 24.0, chr18 = 17.2, chr19 = 26.5, chr20 = 27.5, chr21 = 13.2, chr22 = 14.7) * 1e6
ACRO <- c("chr13", "chr14", "chr15", "chr21", "chr22")   # p arms not on the arrays
FOCAL <- data.frame(name = c("MYCN", "RB1"), chr = c("chr2", "chr13"), start = c(15940550, 48877883), end = c(15947004, 49056026))

untar_idats <- function(tar, dir) {
  if (!file.exists(tar)) stop("Missing ", tar, " (run data/download_raw_data.sh)")
  if (!length(list.files(dir, "_Grn\\.idat(\\.gz)?$"))) { dir.create(dir, showWarnings = FALSE, recursive = TRUE); utils::untar(tar, exdir = dir) }
  f <- list.files(dir, "_Grn\\.idat(\\.gz)?$", full.names = TRUE)
  data.frame(basename = sub("_Grn\\.idat(\\.gz)?$", "", f), gsm = sub("^(GSM[0-9]+)_.*$", "\\1", basename(f)))
}

# conumee for one array type: returns per-array segments, bins and fit noise (cached per array)
# All arrays (450K and EPIC) use the 51 450K reference arrays and the 450K annotation, restricted to
# the probes shared with the query. This gives every tumour the same bins and reference; the 13 EPIC
# reference arrays gave much noisier EPIC profiles in a first run (median noise 0.74 vs 0.29 for 450K).
ms_ctrl <- minfi::preprocessIllumina(CopyNeutralIMA::getCopyNeutralRGSet("IlluminaHumanMethylation450k"))
ctrl <- conumee::CNV.load(ms_ctrl)
anno_450k <- conumee::CNV.create_anno(array_type = "450k", exclude_regions = exclude_regions, detail_regions = detail_regions)
run_conumee <- function(mset, cohort) {
  anno <- anno_450k
  # EPIC arrays of different manifest versions (and the references) do not carry every annotated probe;
  # conumee needs intensities for all annotated probes, so the annotation is restricted to shared probes
  shared <- intersect(rownames(mset), rownames(ms_ctrl))
  n0 <- length(anno@probes); anno@probes <- anno@probes[names(anno@probes) %in% shared]
  message("  ", cohort, " (", ncol(mset), " arrays): ", length(anno@probes), " of ", n0, " annotated probes present in query and reference arrays")
  q <- conumee::CNV.load(mset)
  lapply(colnames(mset), function(s) {
    cf <- file.path(cache_dir, paste0(cohort, "_", s, "_ref450k.rds"))
    if (file.exists(cf)) return(readRDS(cf))
    x <- conumee::CNV.segment(conumee::CNV.bin(conumee::CNV.fit(q[s], ctrl, anno)))
    out <- list(id = s, noise = x@fit$noise, seg = conumee::CNV.write(x, what = "segments"),
                bins = conumee::CNV.write(x, what = "bins"))
    saveRDS(out, cf); out
  })
}

# per-array summaries: arm values, FGA, focal values
summarise_array <- function(r) {
  seg <- r$seg; names(seg)[names(seg) == "chrom"] <- "chr"
  seg$len <- seg$loc.end - seg$loc.start
  arms <- unlist(lapply(names(CEN), function(ch) {
    s <- seg[seg$chr == ch, ]
    one <- function(lo, hi) {
      st <- pmax(s$loc.start, lo); en <- pmin(s$loc.end, hi); w <- pmax(0, en - st)
      if (sum(w) == 0) NA_real_ else sum(w * s$seg.mean) / sum(w)
    }
    a <- c(q = one(CEN[[ch]], Inf)); if (!ch %in% ACRO) a <- c(p = one(0, CEN[[ch]]), a)
    setNames(a, paste0(sub("chr", "", ch), names(a)))
  }))
  b <- r$bins; vcol <- names(b)[ncol(b)]
  foc <- vapply(seq_len(nrow(FOCAL)), function(i) {
    k <- b$Chromosome == FOCAL$chr[i] & b$End >= FOCAL$start[i] & b$Start <= FOCAL$end[i]
    if (any(k)) mean(b[[vcol]][k]) else NA_real_
  }, numeric(1))
  fga <- sum(seg$len[abs(seg$seg.mean) >= CNV_SEG_THR]) / sum(seg$len)
  c(noise = r$noise, fga = fga, setNames(foc, FOCAL$name), arms)
}

# -----------------------------------------------------------------------------
# 1. Discovery tumours (GSE267015): IDATs of the tumours used in steps 13-17
# -----------------------------------------------------------------------------
smp <- read.delim(file.path(DIR_METH, "meth_samples.tsv"), stringsAsFactors = FALSE)
smp <- smp[smp$used, ]
idat_d <- untar_idats(file.path(DIR_RAW_METH, "GSE267015_RAW.tar"), file.path(DIR_RAW_METH, "idat"))
smp$basename <- idat_d$basename[match(smp$gsm, idat_d$gsm)]
if (anyNA(smp$basename)) stop("IDATs missing for: ", paste(smp$title[is.na(smp$basename)], collapse = ", "))
summaries <- function(rr) data.frame(gsm = vapply(rr, `[[`, "", "id"), do.call(rbind, lapply(rr, summarise_array)),
                                    check.names = FALSE)
res_d <- dplyr::bind_rows(lapply(split(smp, smp$platform), function(d) {
  rg <- minfi::read.metharray(d$basename, force = TRUE)
  ms <- minfi::preprocessIllumina(rg); colnames(ms) <- d$gsm
  summaries(run_conumee(ms, "GSE267015"))
}))
cnv_d <- dplyr::left_join(smp[, c("title", "gsm", "platform", "cluster_ryl_et_al", "mycn_amplification", "rb1_status")], res_d, by = "gsm")
cnv_d$cohort <- "GSE267015"; cnv_d$tumour <- cnv_d$title
# QC (fixed after inspecting the first run): the EPIC arrays of GSE267015 show a systematic bias against
# the 450K and the EPIC references alike (losses of chromosomes 3-8, gains of 17-22 in most tumours,
# e.g. 6q "loss" in 71%), so copy number is used from 450K arrays only; arrays with fit noise above
# CNV_MAX_NOISE are also excluded
cnv_d$qc_pass <- cnv_d$platform == "450K" & cnv_d$noise <= CNV_MAX_NOISE
message("GSE267015: copy number for ", nrow(cnv_d), " tumours; median fit noise ", round(stats::median(cnv_d$noise), 2))
print(cnv_d |> dplyr::group_by(platform) |> dplyr::summarise(n = dplyr::n(), median_noise = stats::median(noise),
      median_fga = stats::median(fga), pct_1q_gain = 100 * mean(`1q` >= CNV_ARM_THR), pct_16q_loss = 100 * mean(`16q` <= -CNV_ARM_THR)))

# -----------------------------------------------------------------------------
# 2. Held-out Liu tumours (GSE58783): IDATs, or the raw GenomeStudio signal table as fallback
# -----------------------------------------------------------------------------
tumour_id <- function(title) sub("_REP[0-9]+$", "", sub("\\s*\\(.*$", "", trimws(title)))
is_tumour <- function(id) grepl("^RB[A-Za-z]*[0-9]+$", id)
pm <- read_series_pheno(file.path(DIR_RAW_LIU, "GSE58783_series_matrix.txt.gz"))
pm$tumour <- tumour_id(pm$title)
pm <- pm[is_tumour(pm$tumour), ]
idat_l <- tryCatch(untar_idats(file.path(DIR_RAW_LIU, "GSE58783_RAW.tar"), file.path(DIR_RAW_LIU, "idat")),
                   error = function(e) { message("  ", conditionMessage(e)); data.frame(basename = character(), gsm = character()) })
if (nrow(idat_l)) {
  pm$basename <- idat_l$basename[match(pm$gsm, idat_l$gsm)]
  pm <- pm[!is.na(pm$basename), ]
  rg <- minfi::read.metharray(pm$basename, force = TRUE)
  ms <- minfi::preprocessIllumina(rg); colnames(ms) <- pm$gsm
  message("Liu: ", ncol(ms), " tumour arrays read from IDATs")
} else {
  # shortcut: GenomeStudio signal columns, mapped to GEO samples by position as in step 18;
  # upgrade if GEO adds IDATs for every sample
  f_raw <- file.path(DIR_RAW_LIU, "GSE58783_Raw_data_i_methylation_retino_all_samples.csv.gz")
  hdr <- names(data.table::fread(f_raw, nrows = 0))
  ucol <- grep("(Signal_A|Unmethylated)", hdr, value = TRUE); mcol <- grep("(Signal_B|[^n]Methylated)", hdr, value = TRUE)
  stopifnot(length(ucol) > 0, length(ucol) == length(mcol))
  raw <- data.table::fread(f_raw, select = c(hdr[1], ucol, mcol), data.table = FALSE)
  all_pm <- read_series_pheno(file.path(DIR_RAW_LIU, "GSE58783_series_matrix.txt.gz"))
  k <- as.integer(sub("^SAMPLE ([0-9]+).*$", "\\1", ucol)); gsm <- all_pm$gsm[k]
  keep <- gsm %in% pm$gsm
  ms <- minfi::MethylSet(Meth = as.matrix(raw[, mcol[keep]]), Unmeth = as.matrix(raw[, ucol[keep]]),
                         annotation = c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"))
  rownames(ms) <- raw[[1]]; colnames(ms) <- gsm[keep]
  pm <- pm[match(colnames(ms), pm$gsm), ]
  message("Liu: ", ncol(ms), " tumour arrays read from the raw signal table (no IDATs found)")
}
rr_l <- run_conumee(ms, "GSE58783")
sum_l <- summaries(rr_l)
sum_l$tumour <- pm$tumour[match(sum_l$gsm, pm$gsm)]
message("Liu: ", sum(sum_l$noise > CNV_MAX_NOISE), " of ", nrow(sum_l), " arrays excluded (fit noise > ", CNV_MAX_NOISE, ")")
sum_l <- sum_l[sum_l$noise <= CNV_MAX_NOISE, ]
num <- setdiff(names(sum_l), c("gsm", "tumour"))
cnv_l <- sum_l |> dplyr::group_by(tumour) |> dplyr::summarise(dplyr::across(dplyr::all_of(num), mean), n_arrays = dplyr::n(), .groups = "drop")
cnv_l$cohort <- "GSE58783"; cnv_l$platform <- "450K"; cnv_l$qc_pass <- TRUE
message("Liu: copy number for ", nrow(cnv_l), " tumours")

# segments of both cohorts (from the cache)
seg_all <- dplyr::bind_rows(
  lapply(list.files(cache_dir, "_ref450k\\.rds$", full.names = TRUE), function(f) { r <- readRDS(f)
    cbind(cohort = sub("_.*$", "", basename(f)), sample = r$id, r$seg) }))
data.table::fwrite(seg_all, file.path(DIR_CNV, "cnv_segments.tsv.gz"), sep = "\t")

cnv <- dplyr::bind_rows(cnv_d, cnv_l)
arm_cols <- grep("^[0-9]+[pq]$", names(cnv), value = TRUE)
save_tsv(cnv, file.path(DIR_CNV, "cnv_samples.tsv"))

# self-check: the arm splitter returns 39 arms (22 q + 17 p)
stopifnot(length(arm_cols) == 39)

# -----------------------------------------------------------------------------
# 3. Arm frequencies per cohort
# -----------------------------------------------------------------------------
cnv$group <- ifelse(cnv$qc_pass, cnv$cohort, paste(cnv$cohort, cnv$platform, "excluded"))
freq <- dplyr::bind_rows(lapply(split(cnv, cnv$group), function(d) data.frame(
  cohort = d$group[1], arm = arm_cols, n = nrow(d),
  pct_gain = 100 * colMeans(d[, arm_cols] >= CNV_ARM_THR, na.rm = TRUE),
  pct_loss = 100 * colMeans(d[, arm_cols] <= -CNV_ARM_THR, na.rm = TRUE))))
save_tsv(freq, file.path(DIR_CNV, "cnv_arm_frequency.tsv"))
print(freq[freq$arm %in% CNV_KEY_ARMS, ])

# -----------------------------------------------------------------------------
# 4. Discovery: copy number vs phenotype score, Ryl cluster and the MYCN annotation
# -----------------------------------------------------------------------------
ps <- read.delim(file.path(DIR_METH, "meth_phenotype_samples.tsv"), stringsAsFactors = FALSE)
cnv_d$score <- ps$score[match(cnv_d$title, ps$title)]
vars <- c("fga", CNV_KEY_ARMS, "MYCN", "RB1")
td <- dplyr::bind_rows(lapply(vars, function(v) {
  all_ok <- !is.na(cnv_d$score) & !is.na(cnv_d[[v]])
  rho_all <- unname(stats::cor(cnv_d[[v]][all_ok], cnv_d$score[all_ok], method = "spearman"))   # transparency: all arrays
  ok <- all_ok & cnv_d$qc_pass
  ct <- stats::cor.test(cnv_d[[v]][ok], cnv_d$score[ok], method = "spearman", exact = FALSE)
  bc <- cnv_d$cluster_ryl_et_al %in% c("B", "C")
  wt <- stats::wilcox.test(cnv_d[[v]][bc & cnv_d$qc_pass], cnv_d[[v]][!bc & cnv_d$qc_pass])
  data.frame(variable = v, n = sum(ok), rho_phenotype = unname(ct$estimate), p_phenotype = ct$p.value,
             rho_all_arrays_incl_EPIC = rho_all, n_cluster_BC = sum(bc & cnv_d$qc_pass),
             median_cluster_A = stats::median(cnv_d[[v]][!bc & cnv_d$qc_pass], na.rm = TRUE),
             median_cluster_BC = stats::median(cnv_d[[v]][bc & cnv_d$qc_pass], na.rm = TRUE), p_cluster = wt$p.value)
}))
td$padj_phenotype <- p.adjust(td$p_phenotype, "BH")
mycn_call <- cnv_d$MYCN >= CNV_AMP_THR; mycn_ann <- cnv_d$mycn_amplification == "MYCN_amplification"
td_mycn <- data.frame(variable = "MYCN call vs annotation (all arrays)", n = nrow(cnv_d),
                      both = sum(mycn_call & mycn_ann, na.rm = TRUE), call_only = sum(mycn_call & !mycn_ann, na.rm = TRUE),
                      annotation_only = sum(!mycn_call & mycn_ann, na.rm = TRUE))
save_tsv(dplyr::bind_rows(td, td_mycn), file.path(DIR_CNV, "cnv_tests_discovery.tsv")); print(td); print(td_mycn)

# -----------------------------------------------------------------------------
# 5. Held-out: copy number vs subtype and high-risk histology (same tests as step 18)
# -----------------------------------------------------------------------------
lab <- read.delim(file.path(DIR_VAL, "liu_samples.tsv"), stringsAsFactors = FALSE)
outc <- c("subtype2", "high_risk_broad", "high_risk_strict", "high_risk_broad_in_subtype2", "high_risk_broad_in_subtype1")
tl <- dplyr::bind_rows(lapply(vars, function(v) dplyr::bind_rows(lapply(outc, function(o) {
  d <- lab[match(cnv_l$tumour, lab$tumour), ]
  y <- if (o == "subtype2") as.integer(d$subtype == 2) else as.integer(d[[o]])
  x <- cnv_l[[v]]; if (v %in% c("16q", "RB1")) x <- -x          # oriented: higher = expected in aggressive
  ok <- !is.na(y) & !is.na(x)
  if (sum(y[ok]) < 3 || sum(1 - y[ok]) < 3) return(NULL)
  r <- pROC::roc(y[ok], x[ok], levels = c(0, 1), direction = "<", quiet = TRUE); ci <- auc_ci(r)
  data.frame(variable = v, outcome = o, n = sum(ok), n_pos = sum(y[ok]), auc = ci[["auc"]], auc_lo = ci[["lo"]],
             auc_hi = ci[["hi"]], wilcox_p = stats::wilcox.test(x[ok][y[ok] == 1], x[ok][y[ok] == 0])$p.value)
}))))
tl <- tl |> dplyr::group_by(outcome) |> dplyr::mutate(padj = p.adjust(wilcox_p, "BH")) |> dplyr::ungroup()
save_tsv(tl, file.path(DIR_CNV, "cnv_tests_liu.tsv")); print(as.data.frame(tl))
message("\n[20_cnv_from_methylation.R] Complete.")
