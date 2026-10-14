# =============================================================================
# Step 06b — Independent check of the malignant-cell calls with CopyKAT
# Input : results/03_scrna/integrated_seurat.rds (step 06)
# Output: results/03_scrna/copykat_calls.tsv          per cell: step-06 call vs CopyKAT
#         results/03_scrna/cnv_caller_agreement.tsv   per tumour (and overall): % agreement,
#                                                     Cohen's kappa, % malignant by each caller
#         results/03_scrna/cnv_caller_by_celltype.tsv % called malignant per cell type
#         results/03_scrna/copykat_seed_summary.tsv   per seed and tumour: separated or not, agreement, kappa
#   The first three files are from the seed-42 run; CopyKAT is repeated with COPYKAT_SEEDS (round-2
#   revision) because repeated runs with an identical seed gave different results (unseeded parallel steps).
# CopyKAT (Gao et al. 2021) segments expression along the genome (Bayesian segmentation)
# and classifies cells as aneuploid / diploid. It is run per tumour on raw counts with
# known normal cells = immune/stromal cells pooled over all tumours of the same data set
# (same chemistry and batch; at most COPYKAT_MAX_REF, random). Several tumours have fewer
# than 10 such cells of their own, too few for CopyKAT to anchor its diploid cluster.
# At most COPYKAT_MAX_CELLS non-reference cells per tumour are used (random, seed 42).
# When CopyKAT's diploid and aneuploid consensus profiles correlate >= 0.6 it labels every
# cell diploid; such tumours are reported as "no separation" and left out of agreement.
# Install once: remotes::install_github("navinlabcode/copykat")
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
if (length(INVASION_EXCLUDE) > 0) stop("Unset INVASION_EXCLUDE.")

load_pkg("Seurat", bioc = FALSE)
load_pkg("Matrix", bioc = FALSE)
load_pkg("dplyr")
if (!requireNamespace("copykat", quietly = TRUE))
  stop('CopyKAT is not installed: remotes::install_github("navinlabcode/copykat")')
# Attach (not just ::) because copykat's gene annotation looks up its bundled data
# (full.anno) on the search path
suppressPackageStartupMessages(library(copykat))

COPYKAT_MAX_CELLS <- 1000L
COPYKAT_MAX_REF   <- 300L
COPYKAT_CORES     <- max(1L, parallel::detectCores() - 1L)
CNV_REFERENCE     <- c("Microglia_TAM", "T_NK", "B_plasma", "Endothelial", "Fibroblast")  # as step 06

# Cohen's kappa for two logical vectors
kappa2 <- function(a, b) {
  po <- mean(a == b); pe <- mean(a) * mean(b) + mean(!a) * mean(!b)
  if (pe == 1) return(NA_real_)
  (po - pe) / (1 - pe)
}
stopifnot(kappa2(c(TRUE, FALSE, TRUE, FALSE), c(TRUE, FALSE, TRUE, FALSE)) == 1,
          kappa2(c(TRUE, TRUE, FALSE, FALSE), c(TRUE, FALSE, TRUE, FALSE)) == 0)

seu <- readRDS(file.path(DIR_SC, "integrated_seurat.rds"))
counts <- SeuratObject::GetAssayData(seu, layer = "counts")
md <- seu@meta.data
md$cell_id <- colnames(seu)
md$ref <- md$cluster_type %in% CNV_REFERENCE

run_copykat <- function(tid, seed) {
  m <- md[md$tumour_id == tid, ]
  set.seed(seed)
  non_ref <- m$cell_id[!m$ref]
  if (length(non_ref) > COPYKAT_MAX_CELLS) non_ref <- sample(non_ref, COPYKAT_MAX_CELLS)
  ref_pool <- md$cell_id[md$ref & md$dataset == m$dataset[1]]
  if (length(ref_pool) > COPYKAT_MAX_REF) ref_pool <- sample(ref_pool, COPYKAT_MAX_REF)
  cells <- c(ref_pool, non_ref)
  raw <- as.matrix(counts[, cells, drop = FALSE])
  safe <- paste0("c", seq_along(cells))              # CopyKAT rewrites some characters in names
  colnames(raw) <- safe
  norm <- if (length(ref_pool) >= 5) safe[seq_along(ref_pool)] else ""
  wd <- file.path(tempdir(), paste0("copykat_", make.names(tid), "_", seed))
  dir.create(wd, showWarnings = FALSE); old <- setwd(wd); on.exit(setwd(old))
  message("  ", tid, ": ", length(non_ref), " non-reference + ", length(ref_pool), " reference cells (",
          sum(ref_pool %in% m$cell_id), " from this tumour)",
          if (identical(norm, "")) " (reference inferred by CopyKAT)" else "")
  ck <- tryCatch(
    copykat::copykat(rawmat = raw, id.type = "S", ngene.chr = 5, win.size = 25, KS.cut = 0.1,
                     sam.name = make.names(tid), distance = "euclidean", norm.cell.names = norm,
                     output.seg = FALSE, plot.genes = FALSE, genome = "hg20", n.cores = COPYKAT_CORES),
    error = function(e) { message("    CopyKAT failed for ", tid, ": ", conditionMessage(e)); NULL })
  pred <- if (is.null(ck)) rep(NA_character_, length(cells)) else
    ck$prediction$copykat.pred[match(safe, ck$prediction$cell.names)]
  data.frame(cell_id = cells, tumour_id = tid, copykat_pred = pred, own_cell = cells %in% m$cell_id)
}

analyse_seed <- function(seed) {
  message("=== CopyKAT per tumour, seed ", seed, " ===")
  ck <- dplyr::bind_rows(lapply(unique(md$tumour_id), run_copykat, seed = seed)) |>
    dplyr::filter(own_cell) |> dplyr::select(-own_cell) |>
    dplyr::left_join(md[, c("cell_id", "dataset", "cluster_type", "ref", "malignant", "cnv_residual_r")], by = "cell_id")
  if (all(is.na(ck$copykat_pred))) stop("CopyKAT produced no calls for any tumour (see messages above).")

  # Agreement on non-reference cells with a definite CopyKAT call, in tumours CopyKAT separated
  status <- ck |> dplyr::filter(!ref) |> dplyr::group_by(tumour_id) |>
    dplyr::summarise(copykat_status = dplyr::case_when(
      all(is.na(copykat_pred)) ~ "failed",
      !any(copykat_pred %in% "aneuploid") ~ "no separation (all diploid)",
      TRUE ~ "separated"), .groups = "drop")
  ev <- ck[!ck$ref & ck$copykat_pred %in% c("aneuploid", "diploid") &
           ck$tumour_id %in% status$tumour_id[status$copykat_status == "separated"], ]
  ev$copykat_malignant <- ev$copykat_pred == "aneuploid"
  agree_of <- function(d) data.frame(n_cells = nrow(d),
    pct_malignant_step06 = 100 * mean(d$malignant), pct_aneuploid_copykat = 100 * mean(d$copykat_malignant),
    pct_agreement = 100 * mean(d$malignant == d$copykat_malignant), kappa = kappa2(d$malignant, d$copykat_malignant))
  agr <- dplyr::bind_rows(
    dplyr::bind_rows(lapply(split(ev, ev$tumour_id), function(d) cbind(tumour_id = d$tumour_id[1], agree_of(d)))),
    cbind(tumour_id = "ALL", agree_of(ev)))
  agr <- dplyr::full_join(status, agr, by = "tumour_id")
  agr <- dplyr::left_join(agr, ck |> dplyr::group_by(tumour_id) |>
                            dplyr::summarise(n_not_called = sum(!ref & !copykat_pred %in% c("aneuploid", "diploid")),
                                             .groups = "drop"), by = "tumour_id")

  by_ct <- ev |> dplyr::group_by(cluster_type) |>
    dplyr::summarise(n = dplyr::n(), pct_malignant_step06 = 100 * mean(malignant),
                     pct_aneuploid_copykat = 100 * mean(copykat_malignant), .groups = "drop")
  list(ck = ck, agr = agr, by_ct = by_ct)
}

runs <- lapply(setNames(COPYKAT_SEEDS, COPYKAT_SEEDS), analyse_seed)
r42 <- runs[[as.character(GLOBAL_SEED)]]
save_tsv(r42$ck, file.path(DIR_SC, "copykat_calls.tsv"))
save_tsv(r42$agr, file.path(DIR_SC, "cnv_caller_agreement.tsv")); print(r42$agr)
save_tsv(r42$by_ct, file.path(DIR_SC, "cnv_caller_by_celltype.tsv")); print(r42$by_ct)
seed_sum <- dplyr::bind_rows(lapply(names(runs), function(sd) cbind(seed = as.integer(sd), runs[[sd]]$agr)))
save_tsv(seed_sum, file.path(DIR_SC, "copykat_seed_summary.tsv"))
ov <- seed_sum[seed_sum$tumour_id == "ALL", ]
ov$n_separated <- vapply(ov$seed, function(sd) sum(seed_sum$seed == sd & seed_sum$copykat_status %in% "separated"), 0L)
print(ov[, c("seed", "n_separated", "pct_agreement", "kappa")])
message("Across ", length(COPYKAT_SEEDS), " seeds: ", min(ov$n_separated), "-", max(ov$n_separated),
        " of ", length(unique(md$tumour_id)), " tumours separated; overall kappa ", round(min(ov$kappa), 2), "-", round(max(ov$kappa), 2))

message("\n[06b_copykat_validation.R] Complete. Overall kappa (seed ", GLOBAL_SEED, "): ", round(r42$agr$kappa[r42$agr$tumour_id == "ALL"], 2))
