# =============================================================================
# Step 19 — Supporting Information tables and session information
# Output: results/supporting_info/S1-S9 *.xlsx, sessionInfo.txt
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
if (length(INVASION_EXCLUDE) > 0) stop("Unset INVASION_EXCLUDE.")
load_pkg("openxlsx", bioc = FALSE)

rd <- function(...) read.delim(file.path(...), stringsAsFactors = FALSE, check.names = FALSE)
excl_dir <- function(x) paste0(DIR_INV, "_excl_", gsub(",", "_", x))
pb <- readRDS(file.path(DIR_SC, "pseudobulk_malignant.rds"))
cand <- rd(DIR_INV, "invasion_candidates.tsv")
pb_c <- function(m) data.frame(gene = intersect(cand$gene, rownames(m)),
                               m[intersect(cand$gene, rownames(m)), , drop = FALSE], check.names = FALSE)

tables <- list(
  S1_Table_sample_labels = list(
    labels = rd(DIR_INV, "sample_labels.tsv"),
    excluded_non_primary = rd(DIR_BULK, "excluded_non_primary_samples.tsv"),
    GSE97508_QC = rd(DIR_INV, "gse97508_qc.tsv")),
  S2_Table_per_contrast_effects = list(per_contrast = rd(DIR_INV, "per_contrast_effects.tsv")),
  S3_Table_meta_analysis = list(all_genes = rd(DIR_INV, "meta_invasion_all.tsv"),
                                hksj_sensitivity = rd(DIR_INV, "hksj_sensitivity.tsv")),
  S4_Table_sensitivity = setNames(lapply(SENS_EXCLUDES, function(x) rd(excl_dir(x), "meta_invasion_all.tsv")),
                                  paste0("without_", gsub(",", "_", SENS_EXCLUDES))),
  S5_Table_priority = list(
    candidates = rd(DIR_PRI, "priority_all_candidates.tsv"),
    scenarios = rd(DIR_PRI, "shortlist_scenarios.tsv"),
    bootstrap = rd(DIR_PRI, "bootstrap_shortlist_probability.tsv")),
  S6_Table_single_cell = list(
    composition = rd(DIR_SC, "cell_composition_by_tumour.tsv"),
    cnv_arm_profiles = rd(DIR_SC, "cnv_arm_profiles_by_tumour.tsv"),
    cnv_arm_reference_profiles = rd(DIR_SC, "cnv_arm_reference_profiles.tsv"),
    copykat_agreement = rd(DIR_SC, "cnv_caller_agreement.tsv"),
    copykat_by_celltype = rd(DIR_SC, "cnv_caller_by_celltype.tsv"),
    GSE249995_candidates = rd(DIR_INV, "scrna_gse249995_candidates.tsv"),
    pseudobulk_mean = pb_c(pb$mean), pseudobulk_pct = pb_c(pb$pct),
    lineage_mean = pb_c(pb$lineage_mean), lineage_pct = pb_c(pb$lineage_pct)),
  S7_Table_immune = list(
    mcp_counter_tests = rd(DIR_IMM, "deconv_wilcox_invasion.tsv"),
    purity = rd(DIR_IMM, "estimate_purity.tsv"),
    score_vs_celltypes = rd(DIR_INT, "invasion_immune_corr_celltype.tsv"),
    gene_correlations = rd(DIR_INT, "gene_invasion_immune_corr.tsv"),
    covarying_genes = rd(DIR_INT, "covarying_genes.tsv"),
    # optional: needs the CIBERSORTx web-tool outputs in data/cibersortx/
    cibersortx_fit = if (file.exists(file.path(DIR_IMM, "cibersortx_fit_summary.tsv")))
      rd(DIR_IMM, "cibersortx_fit_summary.tsv") else data.frame(note = "data/cibersortx/ not present")),
  S8_Table_subtype_benchmark = list(
    label_auc = rd(DIR_SUB, "subtype_label_auc.tsv"),
    adjusted_effects = rd(DIR_SUB, "subtype_adjusted_effects.tsv"),
    label_collinearity = rd(DIR_SUB, "subtype_label_collinearity.tsv"),
    signature_genes = rd(DIR_SUB, "signatures.tsv")),
  S10_Table_methylation = list(
    samples = rd(DIR_METH, "meth_samples.tsv"),
    dmp_significant = { d <- data.table::fread(file.path(DIR_METH, "dmp_phenotype.tsv.gz"), data.table = FALSE); d[d$dmp, ] },
    dmp_context = rd(DIR_METH, "dmp_context.tsv"),
    dmr = if (file.exists(file.path(DIR_METH, "dmr_phenotype.tsv"))) rd(DIR_METH, "dmr_phenotype.tsv") else data.frame(note = "DMRcate not run"),
    cpg_regions = rd(DIR_METH, "cpg_regions.tsv"),
    cpg_panel = rd(DIR_METH, "cpg_panel.tsv"),
    cpg_panel_cv = rd(DIR_METH, "cpg_panel_cv.tsv"),
    cpg_panel_scores = rd(DIR_METH, "cpg_panel_scores.tsv")),
  S11_Table_multiomic = list(
    candidate_methylation = rd(DIR_MO, "candidate_methylation_support.tsv"),
    enrichment = rd(DIR_MO, "methylation_support_enrichment.tsv"),
    multiomic_priority = rd(DIR_MO, "multiomic_priority.tsv"),
    multiomic_shortlist = rd(DIR_MO, "multiomic_shortlist.tsv")),
  S12_Table_liu_validation = list(
    samples = rd(DIR_VAL, "liu_samples.tsv"),
    subtype_vs_histology = rd(DIR_VAL, "liu_subtype_vs_histology.tsv"),
    expression_tests = rd(DIR_VAL, "liu_expression_tests.tsv"),
    methylation_tests = rd(DIR_VAL, "liu_methylation_tests.tsv"),
    liu9_in_gse267015 = rd(DIR_VAL, "liu9_in_gse267015.tsv")),
  S9_Table_classifier = list(
    lasso_stability = rd(DIR_INV, "lasso_stability.tsv"),
    loco_auc = rd(DIR_INV, "loco_auc.tsv"),
    score_validation = rd(DIR_INT, "score_validation.tsv"))
)
for (nm in names(tables)) {
  openxlsx::write.xlsx(tables[[nm]], file.path(DIR_SI, paste0(nm, ".xlsx")), overwrite = TRUE)
  message("  Saved: ", nm, ".xlsx (", length(tables[[nm]]), " sheets)")
}
writeLines(capture.output(sessionInfo()), file.path(DIR_SI, "sessionInfo.txt"))
# Versions of every package the analysis uses (quoted in Methods by the manuscript build)
pkgs <- c("GEOquery", "Biobase", "limma", "preprocessCore", "AnnotationDbi", "org.Hs.eg.db",
          "glmnet", "pROC", "Seurat", "SeuratObject", "harmony", "copykat",
          "immunedeconv", "DMRcate", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "minfi", "data.table",
          "readxl", "openxlsx", "dplyr", "ggplot2", "Matrix")
ver <- data.frame(package = pkgs, version = vapply(pkgs, function(p)
  tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_), character(1)))
ver <- rbind(data.frame(package = "R", version = paste(R.version$major, R.version$minor, sep = ".")), ver)
save_tsv(ver, file.path(DIR_SI, "package_versions.tsv"))
message("\n[19_supporting_info.R] Complete.")
