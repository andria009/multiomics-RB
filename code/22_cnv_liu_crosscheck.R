# =============================================================================
# Step 22 — Check of the methylation-array copy number against Liu et al.'s array-CGH
# =============================================================================
# Copy number inferred from the 450K arrays (step 20) is compared, in the same Liu tumours, with the
# independent BAC array-CGH profiles deposited by Liu et al. (GSE58784; GLAD-smoothed log2 ratios).
# Input : results/11_cnv/cnv_samples.tsv (step 20)
#         data/raw/liu/GSE58784-GPL18861_series_matrix.txt.gz, GSE58784-GPL29807_series_matrix.txt.gz
#         data/raw/liu/GPL18861_family.soft.gz, GPL29807_family.soft.gz (BAC positions, hg19)
# Output: results/11_cnv/cnv_liu_bac_samples.tsv   per-tumour BAC arm values, FGA, MYCN, RB1
#         results/11_cnv/cnv_liu_crosscheck.tsv    per variable: Spearman rho, call agreement, kappa
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
load_pkg("GEOquery", bioc = TRUE)

# same arm split and focal regions as step 20
CEN <- c(chr1 = 125.0, chr2 = 93.3, chr3 = 91.0, chr4 = 50.4, chr5 = 48.4, chr6 = 61.0, chr7 = 59.9, chr8 = 45.6,
         chr9 = 49.0, chr10 = 40.2, chr11 = 53.7, chr12 = 35.8, chr13 = 17.9, chr14 = 17.6, chr15 = 19.0,
         chr16 = 36.6, chr17 = 24.0, chr18 = 17.2, chr19 = 26.5, chr20 = 27.5, chr21 = 13.2, chr22 = 14.7) * 1e6
ACRO <- c("chr13", "chr14", "chr15", "chr21", "chr22")
FOCAL <- data.frame(name = c("MYCN", "RB1"), chr = c("chr2", "chr13"), start = c(15940550, 48877883), end = c(15947004, 49056026))

kappa_cat <- function(a, b) {
  po <- mean(a == b); lv <- union(a, b)
  pe <- sum(vapply(lv, function(k) mean(a == k) * mean(b == k), 0))
  (po - pe) / (1 - pe)
}
stopifnot(kappa_cat(c(1, 0, -1, 0), c(1, 0, -1, 0)) == 1, abs(kappa_cat(c(1, 1, 0, 0), c(1, 0, 1, 0))) < 1e-12)

read_gpl <- function(f) {
  L <- readLines(gzfile(f))
  i <- grep("^!platform_table_begin", L); j <- grep("^!platform_table_end", L)
  g <- read.delim(text = L[(i + 1):(j - 1)], stringsAsFactors = FALSE, check.names = FALSE)
  names(g) <- tolower(names(g))
  data.frame(id = as.character(g$id), chr = paste0("chr", g$chromosome), start = as.numeric(g$start), end = as.numeric(g$end))
}

bac_summary <- function(platform) {
  es <- GEOquery::getGEO(filename = file.path(DIR_RAW_LIU, paste0("GSE58784-", platform, "_series_matrix.txt.gz")), getGPL = FALSE)
  gpl <- read_gpl(file.path(DIR_RAW_LIU, paste0(platform, "_family.soft.gz")))
  X <- Biobase::exprs(es); pos <- gpl[match(rownames(X), gpl$id), ]
  keep <- !is.na(pos$chr) & pos$chr %in% names(CEN); X <- X[keep, , drop = FALSE]; pos <- pos[keep, ]
  mid <- (pos$start + pos$end) / 2
  tumour <- sub("\\s*\\(.*$", "", trimws(Biobase::pData(es)$title))
  out <- lapply(seq_len(ncol(X)), function(s) {
    v <- X[, s]; ok <- !is.na(v)
    arms <- unlist(lapply(names(CEN), function(ch) {
      a <- c(q = mean(v[ok & pos$chr == ch & mid > CEN[[ch]]]))
      if (!ch %in% ACRO) a <- c(p = mean(v[ok & pos$chr == ch & mid <= CEN[[ch]]]), a)
      setNames(a, paste0(sub("chr", "", ch), names(a)))
    }))
    foc <- vapply(seq_len(nrow(FOCAL)), function(i) {
      k <- ok & pos$chr == FOCAL$chr[i] & pos$end >= FOCAL$start[i] & pos$start <= FOCAL$end[i]
      if (!any(k)) k <- ok & pos$chr == FOCAL$chr[i] & abs(mid - (FOCAL$start[i] + FOCAL$end[i]) / 2) < 1e6
      if (any(k)) mean(v[k]) else NA_real_
    }, 0)
    # BAC clones are roughly evenly spaced, so the fraction of clones beyond the threshold approximates FGA
    c(n_bac = sum(ok), fga = mean(abs(v[ok]) >= CNV_SEG_THR), setNames(foc, FOCAL$name), arms)
  })
  data.frame(tumour = tumour, platform = platform, do.call(rbind, out), check.names = FALSE)
}

bac <- rbind(bac_summary("GPL18861"), bac_summary("GPL29807"))
bac <- bac[!duplicated(bac$tumour), ]
save_tsv(bac, file.path(DIR_CNV, "cnv_liu_bac_samples.tsv"))

ours <- read.delim(file.path(DIR_CNV, "cnv_samples.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
ours <- ours[ours$cohort == "GSE58783", ]
common <- intersect(ours$tumour, bac$tumour)
A <- ours[match(common, ours$tumour), ]; B <- bac[match(common, bac$tumour), ]
arm_cols <- grep("^[0-9]+[pq]$", names(bac), value = TRUE)
stopifnot(length(arm_cols) == 39, all(arm_cols %in% names(ours)))
message("Liu tumours with both 450K-derived and BAC-CGH copy number: ", length(common))

call3 <- function(x) ifelse(x >= CNV_ARM_THR, 1, ifelse(x <= -CNV_ARM_THR, -1, 0))
res <- do.call(rbind, lapply(c("fga", "MYCN", "RB1", arm_cols), function(v) {
  ok <- !is.na(A[[v]]) & !is.na(B[[v]])
  ct <- suppressWarnings(stats::cor.test(A[[v]][ok], B[[v]][ok], method = "spearman", exact = FALSE))
  is_arm <- v %in% arm_cols
  data.frame(variable = v, n = sum(ok), rho = unname(ct$estimate), p = ct$p.value,
             pct_agree = if (is_arm) 100 * mean(call3(A[[v]][ok]) == call3(B[[v]][ok])) else NA_real_,
             kappa = if (is_arm) kappa_cat(call3(A[[v]][ok]), call3(B[[v]][ok])) else NA_real_,
             pct_altered_450k = if (is_arm) 100 * mean(call3(A[[v]][ok]) != 0) else NA_real_,
             pct_altered_bac = if (is_arm) 100 * mean(call3(B[[v]][ok]) != 0) else NA_real_)
}))
# per-tumour agreement of the 39-arm profiles
prof <- vapply(seq_along(common), function(i)
  suppressWarnings(stats::cor(unlist(A[i, arm_cols]), unlist(B[i, arm_cols]), method = "spearman", use = "complete.obs")), 0)
res <- rbind(res, data.frame(variable = "per-tumour 39-arm profile (median rho)", n = length(common), rho = stats::median(prof, na.rm = TRUE),
                             p = NA, pct_agree = NA, kappa = NA, pct_altered_450k = NA, pct_altered_bac = NA))
save_tsv(res, file.path(DIR_CNV, "cnv_liu_crosscheck.tsv"))
print(res[res$variable %in% c("fga", "MYCN", "RB1", CNV_KEY_ARMS) | grepl("profile", res$variable), ], digits = 2)
message("Median rho over the 39 arms: ", round(stats::median(res$rho[res$variable %in% arm_cols]), 2))
message("\n[22_cnv_liu_crosscheck.R] Complete.")
