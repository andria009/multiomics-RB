# =============================================================================
# Step 00 (run once, not part of run_all.sh) — What do the other public RB
# methylation series contain?
# Checks GSE57362 and GSE58783: platform, number of samples, and every sample
# characteristic, to decide whether either carries an aggressive-phenotype label
# or a paired transcriptome that could replicate the GSE267015 analysis.
# Output: results/08_methylation/scout_<GSE>_samples.tsv and a printed summary
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
load_pkg("GEOquery", bioc = TRUE)
options(timeout = 1200)

for (gse in c("GSE57362", "GSE58783")) {
  message("\n=== ", gse, " ===")
  g <- tryCatch(GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE, destdir = tempdir()),
                error = function(e) { message("  could not download: ", conditionMessage(e)); NULL })
  if (is.null(g)) next
  for (es in g) {
    p <- Biobase::pData(es)
    keep <- grep("^title$|^source_name|^characteristics|^platform_id$|^geo_accession$", names(p), value = TRUE)
    p <- p[, keep, drop = FALSE]
    message("  platform ", unique(p$platform_id), ": ", nrow(p), " samples; ",
            nrow(Biobase::exprs(es)), " features in the series matrix")
    for (col in grep("^source_name|^characteristics", names(p), value = TRUE))
      print(utils::head(sort(table(p[[col]]), decreasing = TRUE), 15))
    save_tsv(p, file.path(DIR_METH, paste0("scout_", gse, "_", unique(p$platform_id), "_samples.tsv")))
  }
}
message("\n[00_scout_methylation_cohorts.R] Complete.")
