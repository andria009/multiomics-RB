# =============================================================================
# Step 06 — scRNA-seq: retinoblastoma samples only, integration, malignant-cell calls
# Input : data/raw/scrna/<GSE>_RAW.tar — GSE168434 (Wu 2022; 7 tumours, RB01-03 in
#         technical duplicate), GSE249995 (Zhang 2025; 2 intraocular + 2 extraocular),
#         GSE166173 (Collin 2021; fetal retina, normal tissue and 2 RBs — only the
#         2 RB samples are used); extracted in place on first run
# Output: results/03_scrna/integrated_seurat.rds, cell_metadata.tsv,
#         cnv_arm_profiles_by_tumour.tsv, cnv_arm_reference_profiles.tsv,
#         cell_composition_by_tumour.tsv; figures/06_*
# Method: QC; RB samples only; <= 20,000 cells per dataset; Harmony on dataset and
#         tumour; Louvain clusters; non-malignant cell types from marker module
#         scores (layers joined first); malignant cells from inferred arm-level copy
#         number (expression averaged over chromosome arms relative to immune and
#         stromal reference cells, in the spirit of inferCNV)
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
source(here::here("code", "R", "gene_lists.R"))

load_pkg("Seurat",   bioc = FALSE)
load_pkg("harmony",  bioc = FALSE)
load_pkg("Matrix",   bioc = FALSE)
load_pkg("dplyr")
load_pkg("tidyr")
load_pkg("ggplot2")

set.seed(GLOBAL_SEED)

# =============================================================================
# SECTION 1: Download and extract scRNA-seq RAW.tar from GEO
# Each dataset is distributed as a single RAW.tar containing per-sample 10x MEX dirs
# =============================================================================

#' Download and extract GEO RAW.tar supplementary file
#' Returns path to extraction directory, or NULL if failed
download_geo_raw_tar <- function(gse_id, raw_dir) {
  extract_dir <- file.path(raw_dir, gse_id)

  # Already extracted if directory has subdirs with matrix files
  if (dir.exists(extract_dir)) {
    mex_files <- list.files(extract_dir, pattern = "matrix\\.mtx(\\.gz)?$",
                            recursive = TRUE, full.names = FALSE)
    if (length(mex_files) > 0) {
      message("  [CACHE] ", gse_id, " already extracted (", length(mex_files), " samples)")
      return(extract_dir)
    }
  }

  # Build GEO FTP URL
  prefix <- substr(gse_id, 1, nchar(gse_id) - 3)   # e.g. GSE168nnn
  tar_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/", prefix, "nnn/",
    gse_id, "/suppl/", gse_id, "_RAW.tar"
  )
  tar_file <- file.path(raw_dir, paste0(gse_id, "_RAW.tar"))

  if (!file.exists(tar_file)) {
    message("  Downloading ", gse_id, "_RAW.tar (this may take several minutes) ...")
    options(timeout = 3600)
    dl_result <- tryCatch(
      download.file(tar_url, tar_file, mode = "wb", quiet = FALSE),
      error = function(e) { message("  Download failed: ", e$message); return(-1) }
    )
    if (dl_result != 0 || !file.exists(tar_file)) {
      message("  SKIP: Could not download ", gse_id)
      return(NULL)
    }
  }

  message("  Extracting ", tar_file, " ...")
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  untar(tar_file, exdir = extract_dir)

  # Check for nested .tar.gz per-sample archives and extract them
  inner_tars <- list.files(extract_dir, pattern = "\\.tar\\.gz$|\\.tgz$",
                            full.names = TRUE, recursive = FALSE)
  if (length(inner_tars) > 0) {
    message("  Extracting ", length(inner_tars), " per-sample archives ...")
    for (tg in inner_tars) {
      sample_name <- sub("\\.tar\\.gz$|\\.tgz$", "", basename(tg))
      sample_dir  <- file.path(extract_dir, sample_name)
      dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
      untar(tg, exdir = sample_dir)
    }
  }

  # Some datasets store per-sample files at top level with sample prefix in filename
  # Reorganize flat files (barcodes/features/matrix with sample prefix) into subdirs
  flat_matrices <- list.files(extract_dir, pattern = "matrix\\.mtx\\.gz$",
                               full.names = TRUE, recursive = FALSE)
  if (length(flat_matrices) > 0) {
    message("  Reorganizing ", length(flat_matrices), " flat sample files into subdirs ...")
    for (mat_f in flat_matrices) {
      sample_prefix <- sub("_matrix\\.mtx\\.gz$", "", basename(mat_f))
      sdir <- file.path(extract_dir, sample_prefix)
      dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
      # Move/copy matrix, barcodes, features files
      for (suffix in c("_matrix.mtx.gz", "_barcodes.tsv.gz", "_features.tsv.gz",
                        "_barcodes.tsv", "_features.tsv", "_genes.tsv")) {
        src <- file.path(extract_dir, paste0(sample_prefix, suffix))
        if (file.exists(src)) {
          # Rename to standard 10x names
          std_name <- sub(paste0("^", sample_prefix, "_"), "", paste0(sample_prefix, suffix))
          std_name <- sub("^_", "", std_name)
          std_name <- gsub("^genes\\.tsv", "features.tsv", std_name)
          file.copy(src, file.path(sdir, std_name), overwrite = TRUE)
        }
      }
    }
  }

  extract_dir
}

#' Build Seurat object from extracted GEO RAW directory
#' Reads all per-sample 10x MEX subdirectories and merges them
create_seurat_from_raw <- function(gse_id, raw_dir, dataset_label) {
  cache_path <- file.path(DIR_SC, paste0(gse_id, "_seurat.rds"))

  if (file.exists(cache_path)) {
    message("  [CACHE] ", gse_id, " Seurat object")
    return(readRDS(cache_path))
  }

  extract_dir <- download_geo_raw_tar(gse_id, raw_dir)
  if (is.null(extract_dir)) {
    message("  SKIP: ", gse_id, " — no data available")
    return(NULL)
  }

  # Find per-sample directories (must contain matrix + barcodes + features)
  all_subdirs <- list.dirs(extract_dir, recursive = TRUE, full.names = TRUE)
  sample_dirs <- Filter(function(d) {
    has_matrix   <- any(file.exists(file.path(d, c("matrix.mtx.gz", "matrix.mtx"))))
    has_barcodes <- any(file.exists(file.path(d, c("barcodes.tsv.gz", "barcodes.tsv"))))
    has_features <- any(file.exists(file.path(d, c("features.tsv.gz", "features.tsv",
                                                     "genes.tsv.gz",   "genes.tsv"))))
    has_matrix && has_barcodes && has_features
  }, all_subdirs)

  if (length(sample_dirs) == 0) {
    message("  WARNING: No valid 10x MEX directories found in ", extract_dir)
    message("  Contents: ", paste(head(list.files(extract_dir), 10), collapse=", "))
    return(NULL)
  }

  message("  Found ", length(sample_dirs), " sample directories for ", gse_id)

  # Read each sample
  seu_list <- lapply(seq_along(sample_dirs), function(i) {
    sdir   <- sample_dirs[[i]]
    slabel <- basename(sdir)
    message("    Reading sample ", i, "/", length(sample_dirs), ": ", slabel)
    tryCatch({
      counts <- Seurat::Read10X(data.dir = sdir)
      # Handle gene expression only (ignore ADT etc. if multi-modal)
      if (is.list(counts)) counts <- counts[["Gene Expression"]]
      seu_s  <- Seurat::CreateSeuratObject(
        counts       = counts,
        project      = slabel,
        min.cells    = SCRNA_MIN_CELLS,
        min.features = SCRNA_MIN_FEATURES
      )
      seu_s$sample_id <- slabel
      seu_s$dataset   <- dataset_label
      seu_s$gse_id    <- gse_id
      seu_s[["pct_mt"]] <- Seurat::PercentageFeatureSet(seu_s, pattern = "^MT-")
      seu_s
    }, error = function(e) {
      message("    ERROR reading ", slabel, ": ", e$message)
      NULL
    })
  })

  seu_list <- Filter(Negate(is.null), seu_list)
  if (length(seu_list) == 0) {
    message("  ERROR: All samples failed for ", gse_id)
    return(NULL)
  }

  # Merge all samples
  if (length(seu_list) == 1) {
    seu <- seu_list[[1]]
  } else {
    seu <- merge(
      x            = seu_list[[1]],
      y            = seu_list[-1],
      add.cell.ids = sapply(seu_list, function(s) unique(s$sample_id)[1]),
      merge.data   = FALSE
    )
  }

  message("    Raw merged: ", ncol(seu), " cells x ", nrow(seu), " genes")
  saveRDS(seu, cache_path)
  message("  Saved: ", cache_path)
  seu
}

message("=== Loading scRNA-seq datasets ===")
seu_wu    <- create_seurat_from_raw("GSE168434", DIR_RAW_SCRNA, "Wu2022")
seu_zhang <- create_seurat_from_raw("GSE249995", DIR_RAW_SCRNA, "Zhang2025")
seu_col   <- create_seurat_from_raw("GSE166173", DIR_RAW_SCRNA, "Collin2021")
seu_col$dataset <- "Collin2021"   # older caches carry a wrong dataset label

# Retinoblastoma samples only; one tumour ID per tumour (technical replicates share it)
seu_col <- subset(seu_col, cells = colnames(seu_col)[grepl("Retinoblastoma", seu_col$sample_id)])
seu_wu$tumour_id    <- sub("_rep[0-9]+.*|_filtered.*", "", seu_wu$sample_id)
seu_zhang$tumour_id <- seu_zhang$sample_id
seu_col$tumour_id   <- seu_col$sample_id
for (s in list(seu_wu, seu_zhang, seu_col)) print(table(s$tumour_id))

# GSE249995 site: _in = intraocular, _ex = extraocular (Zhang 2025)
seu_zhang$invasion_status <- dplyr::case_when(
  grepl("_in[0-9]", seu_zhang$sample_id) ~ "non-invasive",
  grepl("_ex[0-9]", seu_zhang$sample_id) ~ "invasive",
  TRUE ~ NA_character_)

# -----------------------------------------------------------------------------
# QC, per-dataset preprocessing and downsampling
# -----------------------------------------------------------------------------
qc <- function(seu, label) {
  before <- ncol(seu)
  seu <- subset(seu, subset = nFeature_RNA > SCRNA_MIN_FEATURES &
                  nFeature_RNA < SCRNA_MAX_FEATURES & pct_mt < SCRNA_MAX_MT_PCT)
  message("  QC ", label, ": ", ncol(seu), "/", before, " cells")
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = SCRNA_N_HVG, verbose = FALSE)
  if (ncol(seu) > SCRNA_MAX_CELLS) {
    set.seed(GLOBAL_SEED)
    seu <- seu[, sample(colnames(seu), SCRNA_MAX_CELLS)]
  }
  seu
}
SCRNA_MAX_CELLS <- 20000L
objs <- list(Wu2022 = qc(seu_wu, "Wu2022"), Zhang2025 = qc(seu_zhang, "Zhang2025"),
             Collin2021 = qc(seu_col, "Collin2021"))
rm(seu_wu, seu_zhang, seu_col); invisible(gc())

# -----------------------------------------------------------------------------
# Merge and integrate (Harmony on dataset and tumour)
# -----------------------------------------------------------------------------
shared_hvg <- Seurat::SelectIntegrationFeatures(object.list = unname(objs), nfeatures = SCRNA_N_HVG)
seu_all <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs), merge.data = FALSE)
rm(objs); invisible(gc())
seu_all <- Seurat::NormalizeData(seu_all, verbose = FALSE)
Seurat::VariableFeatures(seu_all) <- shared_hvg
seu_all <- Seurat::ScaleData(seu_all, features = shared_hvg, verbose = FALSE)
seu_all <- Seurat::RunPCA(seu_all, features = shared_hvg, npcs = SCRNA_N_DIMS, verbose = FALSE)
seu_all <- harmony::RunHarmony(seu_all, group.by.vars = c("dataset", "tumour_id"),
                               reduction.use = "pca", verbose = FALSE)
seu_all <- Seurat::RunUMAP(seu_all, reduction = "harmony", dims = seq_len(SCRNA_N_DIMS), verbose = FALSE)
seu_all <- Seurat::FindNeighbors(seu_all, reduction = "harmony", dims = seq_len(SCRNA_N_DIMS), verbose = FALSE)
seu_all <- Seurat::FindClusters(seu_all, resolution = SCRNA_RESOLUTION, verbose = FALSE)
seu_all <- SeuratObject::JoinLayers(seu_all)   # one data layer, so every marker is visible
message("  Cells: ", ncol(seu_all), " | clusters: ", length(unique(seu_all$seurat_clusters)))

# -----------------------------------------------------------------------------
# Non-malignant cell types: cluster-level marker module scores
# -----------------------------------------------------------------------------
present <- lapply(CELL_TYPE_MARKERS, intersect, rownames(seu_all))
missing <- lapply(names(CELL_TYPE_MARKERS), function(ct) setdiff(CELL_TYPE_MARKERS[[ct]], present[[ct]]))
message("  Markers absent from the data: ", paste(unlist(missing), collapse = ", "))
present <- present[lengths(present) >= 2]
seu_all <- Seurat::AddModuleScore(seu_all, features = present, name = "ms_", seed = GLOBAL_SEED)
ms_cols <- paste0("ms_", seq_along(present))
cl_scores <- aggregate(seu_all@meta.data[, ms_cols], list(cluster = seu_all$seurat_clusters), mean)
cl_type <- setNames(names(present)[apply(cl_scores[, ms_cols], 1, which.max)], cl_scores$cluster)
seu_all$cluster_type <- unname(cl_type[as.character(seu_all$seurat_clusters)])

# -----------------------------------------------------------------------------
# Malignant cells from inferred arm-level copy number
# -----------------------------------------------------------------------------
load_pkg("org.Hs.eg.db",  bioc = TRUE)
load_pkg("AnnotationDbi", bioc = TRUE)
CNV_REFERENCE <- c("Microglia_TAM", "T_NK", "B_plasma", "Endothelial", "Fibroblast")
X <- SeuratObject::GetAssayData(seu_all, layer = "data")
X <- X[Matrix::rowMeans(X) >= 0.1, ]
loc <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = rownames(X),
                                              keytype = "SYMBOL", columns = "MAP"))
loc <- loc[!is.na(loc$MAP) & !duplicated(loc$SYMBOL), ]
loc$arm <- sub("^([0-9]+)([pq]).*$", "\\1\\2", loc$MAP)
loc <- loc[grepl("^[0-9]+[pq]$", loc$arm), ]
arms <- names(which(table(loc$arm) >= 20))
arms <- arms[order(as.integer(sub("[pq]$", "", arms)), arms)]
loc <- loc[loc$arm %in% arms, ]
A <- Matrix::sparseMatrix(i = match(loc$arm, arms), j = match(loc$SYMBOL, rownames(X)), x = 1,
                          dims = c(length(arms), nrow(X)), dimnames = list(arms, rownames(X)))
A <- A / Matrix::rowSums(A)
ref <- seu_all$cluster_type %in% CNV_REFERENCE
message("  CNV: ", nrow(loc), " genes on ", length(arms), " arms; ", sum(ref), " reference cells")
stopifnot(sum(ref) >= 200)
arm_sc <- as.matrix(A %*% X) - as.vector(A %*% Matrix::rowMeans(X[, ref]))
arm_sc <- sweep(arm_sc, 2, colMeans(arm_sc))           # remove per-cell global level
mag <- sqrt(colMeans(arm_sc^2))
# Arm profiles also carry each cell type's expression programme. Subtract the mean
# profile of the cell's cluster (pooled over all tumours, which Harmony mixes), so what
# remains is tumour-specific (copy-number) deviation. A cell is malignant when its
# residual profile correlates with its tumour's consensus residual more strongly than
# 95% of immune/stromal reference cells do (an empirical null from diploid cells).
cl <- as.character(seu_all$seurat_clusters)
cl_mean <- sapply(split(seq_len(ncol(arm_sc)), cl), function(i) rowMeans(arm_sc[, i, drop = FALSE]))
res <- arm_sc - cl_mean[, cl, drop = FALSE]
res_r <- rep(NA_real_, ncol(arm_sc))
profiles <- list()
for (tid in unique(seu_all$tumour_id)) {
  cells <- which(seu_all$tumour_id == tid)
  cand  <- cells[!ref[cells]]
  profiles[[tid]] <- rowMeans(arm_sc[, cand, drop = FALSE])   # for the heatmap
  cons  <- rowMeans(res[, cand, drop = FALSE])
  res_r[cells] <- as.numeric(stats::cor(res[, cells, drop = FALSE], cons))
}
thr_r <- stats::quantile(res_r[ref], 0.95, na.rm = TRUE)
message("  Malignant threshold (95th percentile of reference cells): r > ", round(thr_r, 3))
seu_all$cnv_magnitude <- mag
seu_all$cnv_residual_r <- res_r
seu_all$malignant <- !ref & !is.na(res_r) & res_r > thr_r
print(round(tapply(seu_all$malignant, seu_all$cluster_type, mean), 3))
seu_all$cell_class <- ifelse(seu_all$malignant, "Malignant", seu_all$cluster_type)
# Arm profiles relative to the immune/stromal reference also carry lineage expression:
# some arms (e.g. 10p low, 12p high) deviate in every tumour. For comparison, the same
# profile is computed for non-malignant retinal cells (rod and Muller-glia clusters, not
# called malignant): arms that deviate there too reflect neural expression, not copy number.
retina_nm <- seu_all$cluster_type %in% c("Rod", "Muller_glia") & !seu_all$malignant
message("  Non-malignant retinal cells for the comparison profile: ", sum(retina_nm))
prof_df <- data.frame(tumour_id = names(profiles), do.call(rbind, profiles), check.names = FALSE)
save_tsv(prof_df, file.path(DIR_SC, "cnv_arm_profiles_by_tumour.tsv"))
ref_prof <- data.frame(profile = c("Non-malignant retinal cells (rod, Muller glia)", "Immune/stromal reference"),
                       rbind(rowMeans(arm_sc[, retina_nm, drop = FALSE]), rowMeans(arm_sc[, ref, drop = FALSE])),
                       check.names = FALSE)
save_tsv(ref_prof, file.path(DIR_SC, "cnv_arm_reference_profiles.tsv"))

comp <- as.data.frame.matrix(table(seu_all$tumour_id, seu_all$cell_class))
comp <- data.frame(tumour_id = rownames(comp), n_cells = rowSums(comp), comp, check.names = FALSE)
comp$dataset <- seu_all$dataset[match(comp$tumour_id, seu_all$tumour_id)]
comp$invasion_status <- seu_all$invasion_status[match(comp$tumour_id, seu_all$tumour_id)]
save_tsv(comp, file.path(DIR_SC, "cell_composition_by_tumour.tsv"))
print(comp)

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------
p1 <- Seurat::DimPlot(seu_all, group.by = "tumour_id", pt.size = 0.2) +
  ggplot2::labs(title = "Tumours (RB samples only)") + PLOT_THEME
p2 <- Seurat::DimPlot(seu_all, group.by = "cell_class", label = TRUE, label.size = 3, pt.size = 0.2) +
  ggplot2::labs(title = "Malignant (inferred CNV) and non-malignant cell types") + PLOT_THEME
save_plot(p1, file.path(DIR_FIG, "06_umap_tumour"), width = 9, height = 7)
save_plot(p2, file.path(DIR_FIG, "06_umap_cellclass"), width = 10, height = 8)

prof_long <- tidyr::pivot_longer(rbind(prof_df, data.frame(tumour_id = "Non-malignant retinal cells", ref_prof[1, -1],
                                                           check.names = FALSE)),
                                 -tumour_id, names_to = "arm", values_to = "score")
prof_long$arm <- factor(prof_long$arm, levels = arms)
p3 <- ggplot2::ggplot(prof_long, ggplot2::aes(arm, tumour_id, fill = score)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(low = COL_DOWNREGULATED, high = COL_UPREGULATED, name = "Relative\nexpression") +
  ggplot2::labs(title = "Inferred arm-level copy number (non-reference cells per tumour; bottom: non-malignant retinal cells)",
                subtitle = "Expected RB changes: 1q and 6p gain, 16q loss", x = NULL, y = NULL) +
  PLOT_THEME + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = 7))
save_plot(p3, file.path(DIR_FIG, "06_cnv_arm_profiles"), width = 12, height = 5)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
save_rds(seu_all, file.path(DIR_SC, "integrated_seurat.rds"), "integrated Seurat object (RB only)")
save_tsv(tibble::rownames_to_column(seu_all@meta.data, "cell_id"), file.path(DIR_SC, "cell_metadata.tsv"))
message("\n[06_scrna_integration.R] Complete. Malignant cells: ", sum(seu_all$malignant),
        " of ", ncol(seu_all))
