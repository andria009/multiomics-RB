# =============================================================================
# Step 07 — Per-tumour pseudobulk within malignant cells (descriptive)
# Input : results/03_scrna/integrated_seurat.rds (step 06),
#         results/02_invasion/invasion_candidates.tsv (step 04)
# Output: results/03_scrna/pseudobulk_malignant.rds (mean expression and % expressing,
#         genes x tumours, in malignant cells and in all tumour-lineage cells),
#         pseudobulk_tumours.tsv,
#         results/02_invasion/scrna_gse249995_candidates.tsv
# Tumours with < 50 cells in a set are not summarised. In GSE249995 (2 extraocular
# vs 2 intraocular tumours) a gene is "concordant" when both extraocular tumours lie
# beyond both intraocular tumours in the bulk direction; a random gene meets this with
# probability 1/6, so it is reported as a description, not as evidence.
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
source(here::here("code", "R", "gene_lists.R"))
if (length(INVASION_EXCLUDE) > 0) stop("Unset INVASION_EXCLUDE.")

load_pkg("Seurat", bioc = FALSE)
load_pkg("Matrix", bioc = FALSE)
load_pkg("dplyr")

seu_all <- readRDS(file.path(DIR_SC, "integrated_seurat.rds"))
X_all <- SeuratObject::GetAssayData(seu_all, layer = "data")

# Per-tumour mean expression and % expressing over a set of cells (tumours with >= 50 cells)
pseudobulk <- function(cells) {
  tid <- seu_all$tumour_id[cells]
  n <- table(tid); tumours <- names(n)[n >= 50]
  grp <- factor(tid, levels = tumours); keep <- !is.na(grp)
  M <- Matrix::sparse.model.matrix(~ 0 + grp, data.frame(grp = grp[keep]))
  M <- M %*% Matrix::Diagonal(x = 1 / Matrix::colSums(M))       # column means
  colnames(M) <- tumours
  Xs <- X_all[, which(cells)[keep], drop = FALSE]
  list(mean = as.matrix(Xs %*% M), pct = as.matrix((Xs > 0) %*% M) * 100, n = n)
}
# Malignant cells (step 06 call) and, less dependent on that call, all tumour-lineage
# (neuronal/photoreceptor-like) cells whether called malignant or not
mal <- pseudobulk(seu_all$malignant)
lin <- pseudobulk(seu_all$cluster_type %in% TUMOUR_LINEAGE)
message("Malignant cells per tumour:"); print(mal$n)
message("Tumour-lineage cells per tumour:"); print(lin$n)
mean_expr <- mal$mean; pct_expr <- mal$pct
tumours <- union(colnames(mal$mean), colnames(lin$mean))

tum <- data.frame(tumour_id = tumours, n_malignant = as.integer(mal$n[tumours]),
                  n_tumour_lineage = as.integer(lin$n[tumours]),
                  dataset = seu_all$dataset[match(tumours, seu_all$tumour_id)],
                  invasion_status = seu_all$invasion_status[match(tumours, seu_all$tumour_id)])
save_rds(list(mean = mean_expr, pct = pct_expr, tumours = tum,
              lineage_mean = lin$mean, lineage_pct = lin$pct),
         file.path(DIR_SC, "pseudobulk_malignant.rds"), "per-tumour pseudobulk: malignant and tumour-lineage cells")
save_tsv(tum, file.path(DIR_SC, "pseudobulk_tumours.tsv"))

# GSE249995 description for every bulk candidate gene
cand <- read.delim(file.path(DIR_INV, "invasion_candidates.tsv"), stringsAsFactors = FALSE)
ex  <- intersect(tum$tumour_id[tum$invasion_status %in% "invasive"], colnames(mean_expr))
inn <- intersect(tum$tumour_id[tum$invasion_status %in% "non-invasive"], colnames(mean_expr))
stopifnot(length(ex) == 2, length(inn) == 2)
g <- intersect(cand$gene, rownames(mean_expr))
up <- cand$g_re[match(g, cand$gene)] > 0
mx <- mean_expr[g, , drop = FALSE]
sc <- data.frame(gene = g, bulk_direction = ifelse(up, "up", "down"),
                 mean_ex1 = mx[, ex[1]], mean_ex2 = mx[, ex[2]],
                 mean_in1 = mx[, inn[1]], mean_in2 = mx[, inn[2]],
                 pct_ex1 = pct_expr[g, ex[1]], pct_ex2 = pct_expr[g, ex[2]],
                 pct_in1 = pct_expr[g, inn[1]], pct_in2 = pct_expr[g, inn[2]])
sc$concordant <- ifelse(up, pmin(sc$mean_ex1, sc$mean_ex2) > pmax(sc$mean_in1, sc$mean_in2),
                            pmax(sc$mean_ex1, sc$mean_ex2) < pmin(sc$mean_in1, sc$mean_in2))
lx <- lin$mean[g, , drop = FALSE]   # same gene rows as the malignant pseudobulk
sc$concordant_lineage <- ifelse(up, pmin(lx[, ex[1]], lx[, ex[2]]) > pmax(lx[, inn[1]], lx[, inn[2]]),
                                    pmax(lx[, ex[1]], lx[, ex[2]]) < pmin(lx[, inn[1]], lx[, inn[2]]))
attr(sc, "tumours") <- c(ex, inn)
save_tsv(sc, file.path(DIR_INV, "scrna_gse249995_candidates.tsv"))
message("Concordant in GSE249995 tumour-lineage cells: ", sum(sc$concordant_lineage, na.rm = TRUE), " of ", nrow(sc))
message("Candidates concordant in GSE249995 malignant cells: ", sum(sc$concordant, na.rm = TRUE),
        " of ", nrow(sc), " (", round(100 * mean(sc$concordant, na.rm = TRUE)), "%; chance = 17%)")
message("\n[07_scrna_validation.R] Complete.")
