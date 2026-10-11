# =============================================================================
# Step 01 — Load (or download) the four invasion-labelled microarray cohorts
# Input : data/raw/microarray/<GSE>.rds  (GEOquery caches shipped with the data)
#         or NCBI GEO if a cache is missing
# Output: results/01_bulk/dataset_metadata.tsv
# The RNA-seq (step 03) and scRNA-seq (step 06) inputs are read directly there.
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))

load_pkg("GEOquery", bioc = TRUE)
options(timeout = 600)

# Platform annotation is needed for probe -> gene symbol mapping in step 02.
# GSE110811 (GPL16686, HuGene 2.0 ST) is mapped from GB_ACC accessions with org.Hs.eg.db
# instead, so its platform table is not downloaded (this matches the shipped caches).
GET_GPL <- c(GSE97508 = TRUE, GSE29683 = TRUE, GSE59983 = TRUE, GSE110811 = FALSE)

download_geo_series <- function(gse_id, destdir) {
  cache_path <- file.path(destdir, paste0(gse_id, ".rds"))
  if (file.exists(cache_path)) {
    message("  [CACHE] Loading ", gse_id, " from ", cache_path)
    return(readRDS(cache_path))
  }
  message("  [DOWNLOAD] Fetching ", gse_id, " ...")
  geo_list <- GEOquery::getGEO(gse_id, destdir = destdir, GSEMatrix = TRUE,
                               getGPL = GET_GPL[[gse_id]], AnnotGPL = FALSE)
  geo_obj <- if (is.list(geo_list) && !is(geo_list, "ExpressionSet")) geo_list[[1]] else geo_list
  result <- list(gse_id  = gse_id,
                 geo_obj = geo_obj,
                 pheno   = Biobase::pData(geo_obj),
                 expr    = Biobase::exprs(geo_obj))
  saveRDS(result, cache_path)
  message("  Cached to: ", cache_path)
  result
}

message("\n=== Invasion-labelled microarray cohorts ===")
geo <- lapply(setNames(GEO_INVASION_MICRO, GEO_INVASION_MICRO),
              download_geo_series, destdir = DIR_RAW_MICRO)

metadata_summary <- dplyr::bind_rows(lapply(geo, function(g) data.frame(
  gse_id    = g$gse_id,
  n_samples = nrow(g$pheno),
  platform  = unique(g$pheno$platform_id)[1]
)))
print(metadata_summary)
save_tsv(metadata_summary, file.path(DIR_BULK, "dataset_metadata.tsv"))

message("\n[01_download_geo.R] Complete.")
