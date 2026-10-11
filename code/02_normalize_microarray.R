# =============================================================================
# Step 02 — Normalise the four invasion-labelled microarray cohorts
# Method: GEO series-matrix values (already processed by the submitters); primary
#         tumours only (config.R PRIMARY_PATTERN: normal retina, retinocytoma, cell
#         lines and xenografts removed before normalisation); probes mapped to gene
#         symbols, log2 where needed, duplicate probes collapsed to the one with the
#         highest IQR, then quantile normalisation (preprocessCore) per dataset.
# Input : data/raw/microarray/<GSE>.rds (step 01)
# Output: results/01_bulk/normalized_matrices.rds, pheno_list.rds
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))

load_pkg("GEOquery",  bioc = TRUE)
load_pkg("Biobase",   bioc = TRUE)
load_pkg("dplyr")

# -----------------------------------------------------------------------------
# 1. Load cached GEO objects (step 01)
# -----------------------------------------------------------------------------
message("Loading cached GEO objects ...")
geo <- lapply(setNames(GEO_INVASION_MICRO, GEO_INVASION_MICRO),
              function(g) readRDS(file.path(DIR_RAW_MICRO, paste0(g, ".rds"))))

# -----------------------------------------------------------------------------
# 2. Extract + clean expression matrix from a GEO object
#    - Maps probe IDs to gene symbols (platform table, or annotation package for
#      Affymetrix Gene ST arrays)
#    - Removes control / unannotated probes
#    - Log2-transforms if not already on a log scale
#    - Collapses duplicate genes to the probe with the highest IQR
# -----------------------------------------------------------------------------
extract_expr <- function(geo_result, keep) {
  gse_id  <- geo_result$gse_id
  geo_obj <- geo_result$geo_obj
  mat     <- Biobase::exprs(geo_obj)[, keep, drop = FALSE]
  fdata   <- Biobase::fData(geo_obj)

  message("\nProcessing ", gse_id, ":")
  message("  Probes: ", nrow(mat), " | Samples: ", ncol(mat))

  # Identify gene symbol column
  sym_col <- grep("gene.symbol|gene_symbol|Symbol|GENE_SYMBOL",
                  names(fdata), ignore.case = TRUE, value = TRUE)

  if (length(sym_col) > 0) {
    sym_col <- sym_col[1]
    message("  Using symbol column: ", sym_col)
    gene_symbols <- fdata[[sym_col]]
  } else if (!is.null(fdata[["GB_ACC"]]) && any(nchar(fdata[["GB_ACC"]]) > 0, na.rm = TRUE)) {
    # Fallback: map GenBank/RefSeq accessions to HGNC symbols via org.Hs.eg.db
    message("  No symbol column found. Mapping GB_ACC -> symbol via org.Hs.eg.db ...")
    load_pkg("org.Hs.eg.db", bioc = TRUE)
    load_pkg("AnnotationDbi", bioc = TRUE)
    acc <- fdata[["GB_ACC"]]
    # Strip version suffix (e.g. NM_001234.3 -> NM_001234)
    acc_clean <- sub("\\.\\d+$", "", acc)
    gene_symbols <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys      = acc_clean,
      column    = "SYMBOL",
      keytype   = "ACCNUM",
      multiVals = "first"
    )
    gene_symbols <- as.character(gene_symbols)
  } else {
    # Detect Affymetrix Gene ST arrays by platform_id and use annotation DB
    plat <- unique(Biobase::pData(geo_obj)$platform_id)[1]
    annot_pkg <- switch(plat,
      "GPL6244"  = "hugene10sttranscriptcluster.db",
      "GPL11532" = "hugene11sttranscriptcluster.db",
      "GPL16686" = "hugene20sttranscriptcluster.db",
      "GPL17930" = "hugene21sttranscriptcluster.db",
      "GPL5175"  = "huex10sttranscriptcluster.db",
      NULL
    )
    if (!is.null(annot_pkg)) {
      message("  Affymetrix Gene ST platform (", plat, "). Using ", annot_pkg, " ...")
      load_pkg(annot_pkg, bioc = TRUE)
      load_pkg("AnnotationDbi", bioc = TRUE)
      gene_symbols <- AnnotationDbi::mapIds(
        get(annot_pkg),
        keys      = rownames(mat),
        column    = "SYMBOL",
        keytype   = "PROBEID",
        multiVals = "first"
      )
      gene_symbols <- as.character(gene_symbols)
    } else {
      # Last resort: check if row names look like gene symbols (contain letters)
      rn <- rownames(mat)
      if (mean(grepl("[A-Za-z]", rn)) > 0.5) {
        message("  No symbol column found. Using row names as gene symbols.")
        gene_symbols <- rn
      } else {
        stop("Cannot find gene symbol column in fData for ", gse_id,
             ". Available: ", paste(names(fdata), collapse = ", "),
             ". Platform: ", plat)
      }
    }
  }

  # Handle multi-gene probe entries (e.g. "GENE1 /// GENE2") — take first
  gene_symbols <- sub(" ///.*", "", gene_symbols)
  gene_symbols <- trimws(gene_symbols)

  # Remove empty / control probes
  valid <- !is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "---"
  mat   <- mat[valid, , drop = FALSE]
  gene_symbols <- gene_symbols[valid]
  rownames(mat) <- gene_symbols

  message("  After probe filter: ", nrow(mat), " probes")

  # Check if already log2-transformed (median > 30 suggests raw intensities)
  if (median(mat, na.rm = TRUE) > 30) {
    message("  Log2-transforming expression values ...")
    mat <- log2(mat + 1)
  }

  # Deduplicate genes by highest IQR
  mat <- dedup_by_iqr(mat)
  message("  After dedup: ", nrow(mat), " unique genes")

  mat
}

# -----------------------------------------------------------------------------
# 3. Quantile-normalise each dataset separately
# -----------------------------------------------------------------------------
quantile_norm <- function(mat) {
  load_pkg("preprocessCore", bioc = TRUE)
  normed <- preprocessCore::normalize.quantiles(mat)
  dimnames(normed) <- dimnames(mat)
  normed
}


# -----------------------------------------------------------------------------
# 4. Sample metadata (title + GEO characteristics columns)
# -----------------------------------------------------------------------------
extract_pheno <- function(geo_result) {
  pheno <- geo_result$pheno
  # Commonly named columns in GEO pheno data
  title_col <- "title"
  char_cols  <- grep("characteristics_ch", names(pheno), value = TRUE)

  df <- data.frame(
    sample_id = rownames(pheno),
    title     = pheno[[title_col]],
    stringsAsFactors = FALSE
  )
  # Append characteristics columns
  for (col in char_cols) {
    df[[col]] <- pheno[[col]]
  }
  df
}

pheno_all <- lapply(geo, extract_pheno)

# Primary tumours only: the pattern is searched in the title and all characteristics
is_primary <- function(p, pattern) {
  txt <- do.call(paste, c(p[, setdiff(names(p), "sample_id"), drop = FALSE], sep = " | "))
  grepl(pattern, txt, ignore.case = TRUE)
}
pheno_list <- lapply(setNames(names(pheno_all), names(pheno_all)), function(g) {
  p <- pheno_all[[g]]; keep <- is_primary(p, PRIMARY_PATTERN[[g]])
  message("  ", g, ": ", sum(keep), " primary tumours of ", nrow(p), " samples")
  stopifnot(sum(keep) >= 6)
  p[keep, , drop = FALSE]
})
excluded <- dplyr::bind_rows(lapply(names(pheno_all), function(g) {
  p <- pheno_all[[g]]; p <- p[!p$sample_id %in% pheno_list[[g]]$sample_id, , drop = FALSE]
  if (nrow(p)) data.frame(dataset = g, sample_id = p$sample_id, title = p$title) else NULL
}))
save_tsv(excluded, file.path(DIR_BULK, "excluded_non_primary_samples.tsv"))

normalized <- lapply(setNames(names(geo), names(geo)), function(g)
  quantile_norm(extract_expr(geo[[g]], pheno_list[[g]]$sample_id)))

# -----------------------------------------------------------------------------
# 5. Save
# -----------------------------------------------------------------------------
save_rds(normalized, file.path(DIR_BULK, "normalized_matrices.rds"), "normalized microarray")
save_rds(pheno_list, file.path(DIR_BULK, "pheno_list.rds"),          "phenotype tables")

for (nm in names(normalized)) {
  message("    ", nm, ": ", nrow(normalized[[nm]]), " genes, ",
          ncol(normalized[[nm]]), " samples")
}
message("\n[02_normalize_microarray.R] Complete.")
