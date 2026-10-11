# =============================================================================
# Step 08 — MCP-counter immune/stromal estimates and ESTIMATE tumour purity (immunedeconv)
# Input : results/01_bulk/ normalised microarrays + GSE268136 log2(TPM+1)
#         results/02_invasion/pheno_{29683,59983}_annotated.rds (step 04)
# Output: results/04_immune/deconv_long_all_datasets.{rds,tsv},
#         deconv_summary_by_group.tsv, deconv_wilcox_invasion.tsv; figures/08_*
#         CIBERSORTx check (not used further): cibersortx_input/mixture_<GSE>.txt (upload
#         files) and cibersortx_fit_summary.tsv, summarising the web-tool results kept in
#         data/cibersortx/ (LM22, 100 permutations; see README.md, Data files)
# Datasets: GSE268136 (Ryl clusters), GSE59983, GSE29683, GSE110811 (invasion)
# Microarray log2 values and RNA-seq log2(TPM+1) are back-transformed to linear scale.
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
if (length(INVASION_EXCLUDE) > 0)
  stop("INVASION_EXCLUDE is set: step 08 uses the main invasion labels. Unset it.")

load_pkg("immunedeconv", bioc = FALSE)
load_pkg("dplyr")
load_pkg("tidyr")
load_pkg("purrr")
load_pkg("ggplot2")
load_pkg("ggpubr",         bioc = FALSE)
load_pkg("tibble")

# =============================================================================
# 1. Load expression matrices and annotate
# =============================================================================
message("Loading expression matrices ...")

normalized <- readRDS(file.path(DIR_BULK, "normalized_matrices.rds"))
pheno_list <- readRDS(file.path(DIR_BULK, "pheno_list.rds"))

# GSE268136 RNA-seq
mat_rnaseq   <- readRDS(file.path(DIR_BULK, "rnaseq_log2tpm_GSE268136.rds"))
pheno_rnaseq <- readRDS(file.path(DIR_BULK, "rnaseq_pheno_GSE268136.rds"))

# Invasion-annotated phenos from step 04 (k-means for GSE29683, Kapatai for GSE59983)
p29_file <- file.path(DIR_INV, "pheno_29683_annotated.rds")
p59_file <- file.path(DIR_INV, "pheno_59983_annotated.rds")
if (!file.exists(p29_file) || !file.exists(p59_file))
  stop("Missing annotated phenos — run code/04_invasion_meta_analysis.R first.")
pheno_29683 <- readRDS(p29_file)
pheno_59983 <- readRDS(p59_file)

# =============================================================================
# 2. Build dataset registry with expression + group annotations
# =============================================================================
# Microarray back-transform: 2^x → linear scale (comparable to TPM-like input)
bt <- function(m) 2^m

build_pheno_110811 <- function() {
  p <- pheno_list$GSE110811
  p$group <- label_invasion_gse110811(p)
  p
}

datasets <- list(
  GSE268136 = list(
    mat       = 2^mat_rnaseq - 1,    # values are log2(TPM+1) despite the 'log1p' file name
    pheno     = pheno_rnaseq |>
                  dplyr::mutate(
                    sample_id = rownames(pheno_rnaseq),
                    group     = dplyr::case_when(
                      ryl_cluster == "A" ~ "Cluster_A",
                      ryl_cluster == "B" ~ "Cluster_B",
                      ryl_cluster == "C" ~ "Cluster_C",
                      is_wt              ~ "WT",
                      TRUE               ~ "unknown"
                    )
                  ),
    group_col = "group",
    type      = "rnaseq",
    label     = "Ryl Cluster (A/B/C)",
    methods   = c("mcp_counter", "estimate")
  ),
  GSE59983  = list(
    mat       = bt(normalized$GSE59983),
    pheno     = pheno_59983 |>
                  dplyr::mutate(group = invasion_status),
    group_col = "group",
    type      = "microarray",
    label     = "Invasion Status (Kapatai groups)",
    methods   = c("mcp_counter", "estimate")
  ),
  GSE29683  = list(
    mat       = bt(normalized$GSE29683),
    pheno     = pheno_29683 |>
                  dplyr::filter(!is.na(invasion_status)) |>
                  dplyr::mutate(group = invasion_status),
    group_col = "group",
    type      = "microarray",
    label     = "Invasion Status (k-means k=3)",
    methods   = c("mcp_counter", "estimate")
  ),
  GSE110811 = list(
    mat       = bt(normalized$GSE110811),
    pheno     = build_pheno_110811(),
    group_col = "group",
    type      = "microarray",
    label     = "Anaplasia Severity",
    methods   = c("mcp_counter", "estimate")
  )
)

# =============================================================================
# 3. Run deconvolution for each dataset
# =============================================================================
run_deconv_dataset <- function(ds_name, ds) {
  message("\n=== ", ds_name, " (", ds$type, ") ===")
  mat   <- ds$mat
  mat[is.na(mat) | mat < 0] <- 0

  results <- list()
  for (method in ds$methods) {
    message("  Running: ", method)
    res <- tryCatch(
      immunedeconv::deconvolute(mat, method),
      error = function(e) { message("    FAILED: ", conditionMessage(e)); NULL }
    )
    if (!is.null(res)) {
      message("    OK — ", nrow(res), " cell types")
      results[[method]] <- res
    }
  }
  results
}

all_deconv <- purrr::imap(datasets, function(ds, ds_name) run_deconv_dataset(ds_name, ds))

# =============================================================================
# 4. Tidy into long format with group annotations
# =============================================================================
tidy_dataset <- function(ds_name, ds, deconv_res) {
  if (length(deconv_res) == 0) return(NULL)

  pheno <- ds$pheno
  # align sample_id column
  if (!"sample_id" %in% names(pheno)) pheno$sample_id <- rownames(pheno)

  purrr::imap_dfr(deconv_res, function(res_df, method_name) {
    res_df |>
      tidyr::pivot_longer(-cell_type, names_to = "sample_id", values_to = "score") |>
      dplyr::mutate(method = method_name, dataset = ds_name) |>
      dplyr::left_join(
        pheno |> dplyr::select(sample_id, group = dplyr::all_of(ds$group_col)),
        by = "sample_id"
      )
  })
}

long_all <- purrr::pmap_dfr(
  list(names(datasets), datasets, all_deconv),
  tidy_dataset
) |> dplyr::filter(!is.na(group))

# MCP-counter's "Monocyte" row repeats the monocytic-lineage score already reported as
# "Macrophage/Monocyte"; keep one so it is not counted twice in BH families.
mono <- long_all |> dplyr::filter(method == "mcp_counter", cell_type %in% c("Monocyte", "Macrophage/Monocyte")) |>
  tidyr::pivot_wider(id_cols = c(dataset, sample_id), names_from = cell_type, values_from = score)
if (all(c("Monocyte", "Macrophage/Monocyte") %in% names(mono)) &&
    isTRUE(all.equal(mono$Monocyte, mono$`Macrophage/Monocyte`))) {
  long_all <- long_all |> dplyr::filter(!(method == "mcp_counter" & cell_type == "Monocyte"))
  message("  Dropped duplicate MCP-counter 'Monocyte' (identical to 'Macrophage/Monocyte')")
}
message("\nTotal long-format rows: ", nrow(long_all))

# =============================================================================
# 5. Per-dataset plots: immune fractions by group
# =============================================================================
plot_deconv_by_group <- function(ds_name, ds, deconv_res) {
  if (length(deconv_res) == 0) return(invisible(NULL))

  pheno <- ds$pheno
  if (!"sample_id" %in% names(pheno)) pheno$sample_id <- rownames(pheno)

  df <- purrr::imap_dfr(deconv_res, function(res_df, method_name) {
    res_df |>
      tidyr::pivot_longer(-cell_type, names_to = "sample_id", values_to = "score") |>
      dplyr::mutate(method = method_name) |>
      dplyr::left_join(
        pheno |> dplyr::select(sample_id, group = dplyr::all_of(ds$group_col)),
        by = "sample_id"
      ) |>
      dplyr::filter(!is.na(group))
  })

  groups       <- sort(unique(df$group))
  comparisons  <- if (length(groups) >= 2) utils::combn(groups, 2, simplify = FALSE) else list()
  n_methods    <- dplyr::n_distinct(df$method)
  n_celltypes  <- dplyr::n_distinct(df$cell_type)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = score, fill = group)) +
    ggplot2::geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.12, size = 0.8, alpha = 0.5) +
    (if (length(comparisons) > 0)
      ggpubr::stat_compare_means(comparisons = comparisons,
                                  method = "wilcox.test", label = "p.signif", size = 2.5)
     else NULL) +
    ggplot2::facet_grid(cell_type ~ method, scales = "free_y") +
    ggplot2::labs(
      title = paste0(ds_name, " — Immune Deconvolution by ", ds$label),
      x = NULL, y = "Score / Fraction"
    ) +
    ggplot2::theme(
      legend.position  = "none",
      strip.text       = ggplot2::element_text(size = 7),
      axis.text.x      = ggplot2::element_text(angle = 30, hjust = 1, size = 7),
      axis.text.y      = ggplot2::element_text(size = 6)
    ) +
    PLOT_THEME

  save_plot(p,
            file.path(DIR_FIG, paste0("08_", ds_name, "_immune_by_group")),
            width  = min(48, max(10, n_methods * 3.5)),
            height = min(48, max(8,  n_celltypes * 1.2)))
  message("  Saved figure: 08_", ds_name, "_immune_by_group")
}

purrr::pwalk(list(names(datasets), datasets, all_deconv), plot_deconv_by_group)

# =============================================================================
# 6. Cross-dataset comparison: MCP-counter T-cell & Macrophage across datasets
#    Focus on invasion-annotated datasets (GSE59983, GSE29683, GSE97508)
# =============================================================================
message("\nBuilding cross-dataset invasion immune comparison ...")

invasion_long <- long_all |>
  dplyr::filter(
    method  == "mcp_counter",
    dataset %in% c("GSE59983", "GSE29683", "GSE110811"),
    group   %in% c("invasive", "non-invasive")
  ) |>
  dplyr::filter(grepl("T cell|Macrophage|NK|B cell|Monocyte|Neutrophil|Dendritic",
                        cell_type, ignore.case = TRUE))

if (nrow(invasion_long) > 0) {
  p_inv_cross <- ggplot2::ggplot(
    invasion_long,
    ggplot2::aes(x = group, y = score, fill = group)
  ) +
    ggplot2::geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.12, size = 1, alpha = 0.6) +
    ggpubr::stat_compare_means(
      comparisons = list(c("invasive","non-invasive")),
      method = "wilcox.test", label = "p.signif", size = 3
    ) +
    ggplot2::facet_grid(cell_type ~ dataset, scales = "free_y") +
    ggplot2::scale_fill_manual(
      values = c("invasive" = COL_INVASIVE, "non-invasive" = COL_NONINVASIVE)
    ) +
    ggplot2::labs(
      title = "MCP-counter: Immune Infiltration — Invasive vs Non-Invasive (cross-dataset)",
      x = NULL, y = "MCP-counter Score"
    ) +
    ggplot2::theme(
      legend.position = "none",
      strip.text      = ggplot2::element_text(size = 8),
      axis.text.x     = ggplot2::element_text(angle = 20, hjust = 1)
    ) +
    PLOT_THEME
  save_plot(p_inv_cross,
            file.path(DIR_FIG, "08_invasion_immune_crossdataset"),
            width = 12, height = 14)
  message("  Saved: 08_invasion_immune_crossdataset")
}

# =============================================================================
# 7. Summary table — median immune fractions per group per dataset (MCP-counter)
# =============================================================================
summary_tbl <- long_all |>
  dplyr::filter(method == "mcp_counter") |>
  dplyr::group_by(dataset, cell_type, group) |>
  dplyr::summarise(
    median_score = round(median(score, na.rm = TRUE), 4),
    n            = dplyr::n(),
    .groups      = "drop"
  )

# Wilcoxon p-values for invasion comparisons
wilcox_tbl <- long_all |>
  dplyr::filter(
    method == "mcp_counter",
    group  %in% c("invasive","non-invasive")
  ) |>
  dplyr::group_by(dataset, cell_type) |>
  dplyr::summarise(
    median_inv    = median(score[group == "invasive"],     na.rm = TRUE),
    median_noninv = median(score[group == "non-invasive"], na.rm = TRUE),
    wilcox_p = tryCatch(
      wilcox.test(score[group=="invasive"], score[group=="non-invasive"])$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) |>
  # direction is reported explicitly: a p-value alone does not say which group is higher
  dplyr::mutate(wilcox_padj = p.adjust(wilcox_p, method = "BH"),
                higher_in   = ifelse(median_inv > median_noninv, "invasive", "non-invasive")) |>
  dplyr::arrange(wilcox_padj)

# Sensitivity: score ~ group + ESTIMATE StromalScore (linear model). The ESTIMATE
# ImmuneScore (and so tumour purity, which is computed from it) shares genes with the
# MCP-counter immune signatures, so adjusting for purity removes part of the signal being
# tested. The StromalScore does not overlap the immune signatures; it does overlap the
# fibroblast/endothelial ones, so those two are not adjusted. Unadjusted tests are primary.
est <- long_all |> dplyr::filter(method == "estimate", cell_type %in% c("tumor purity", "stroma score")) |>
  dplyr::select(dataset, sample_id, cell_type, score) |>
  tidyr::pivot_wider(names_from = cell_type, values_from = score) |>
  dplyr::rename(purity = `tumor purity`, stroma_score = `stroma score`)
STROMAL_TYPES <- c("Endothelial cell", "Cancer associated fibroblast")
adj_tbl <- long_all |>
  dplyr::filter(method == "mcp_counter", group %in% c("invasive", "non-invasive")) |>
  dplyr::inner_join(est, by = c("dataset", "sample_id")) |>
  dplyr::group_by(dataset, cell_type) |>
  dplyr::group_modify(function(d, key) {
    if (key$cell_type %in% STROMAL_TYPES) return(data.frame(beta_invasive_adj = NA_real_, p_adj_stromal = NA_real_))
    f <- tryCatch(summary(stats::lm(score ~ I(group == "invasive") + stroma_score, data = d))$coefficients,
                  error = function(e) NULL)
    if (is.null(f) || nrow(f) < 3) return(data.frame(beta_invasive_adj = NA_real_, p_adj_stromal = NA_real_))
    data.frame(beta_invasive_adj = f[2, 1], p_adj_stromal = f[2, 4])
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(padj_adj_stromal = p.adjust(p_adj_stromal, method = "BH"))
wilcox_tbl <- wilcox_tbl |> dplyr::left_join(adj_tbl, by = c("dataset", "cell_type"))
purity_tbl <- est |>
  dplyr::left_join(long_all |> dplyr::distinct(dataset, sample_id, group), by = c("dataset", "sample_id"))
message("\nESTIMATE purity range per dataset:")
print(purity_tbl |> dplyr::group_by(dataset) |>
        dplyr::summarise(min = min(purity), median = median(purity), max = max(purity), .groups = "drop"))

message("\nTop significant immune differences (MCP-counter, invasive vs non-invasive):")
print(head(wilcox_tbl |> dplyr::filter(!is.na(wilcox_p)), 15))

# =============================================================================
# 8. Save
# =============================================================================
save_rds(all_deconv,   file.path(DIR_IMM, "deconv_results_all_datasets.rds"),
         "MCP-counter results — all datasets")
save_rds(long_all,     file.path(DIR_IMM, "deconv_long_all_datasets.rds"),
         "immune deconvolution long-format — all datasets")
save_tsv(long_all,     file.path(DIR_IMM, "deconv_long_all_datasets.tsv"))
save_tsv(summary_tbl,  file.path(DIR_IMM, "deconv_summary_by_group.tsv"))
save_tsv(wilcox_tbl,   file.path(DIR_IMM, "deconv_wilcox_invasion.tsv"))
save_tsv(purity_tbl,   file.path(DIR_IMM, "estimate_purity.tsv"))

# =============================================================================
# 9. CIBERSORTx check (web tool; reported in Methods, not used further)
#    Mixture files are the same linear-scale matrices given to MCP-counter.
#    The downloaded results are kept in data/cibersortx/.
# =============================================================================
cbx_in <- file.path(DIR_IMM, "cibersortx_input")
dir.create(cbx_in, showWarnings = FALSE)
for (ds in names(datasets)) {
  m <- datasets[[ds]]$mat
  m[is.na(m) | m < 0] <- 0
  readr::write_tsv(tibble::rownames_to_column(as.data.frame(m), "GeneSymbol"),
                   file.path(cbx_in, paste0("mixture_", ds, ".txt")))
}
message("  CIBERSORTx mixtures written to: ", cbx_in)

cbx_fit <- dplyr::bind_rows(lapply(names(datasets), function(ds) {
  f <- here::here("data", "cibersortx", paste0("CIBERSORTx_", ds, ".txt"))
  if (!file.exists(f)) return(NULL)
  x <- read.delim(f, check.names = FALSE)
  stopifnot(all(c("Mixture", "P-value", "Correlation") %in% names(x)))
  data.frame(dataset = ds, n_samples = nrow(x), n_p_below_0.05 = sum(x$`P-value` < 0.05),
             median_correlation = round(median(x$Correlation), 3))
}))
if (nrow(cbx_fit) > 0) {
  print(cbx_fit)
  message("  CIBERSORTx: ", sum(cbx_fit$n_p_below_0.05), " of ", sum(cbx_fit$n_samples),
          " samples with p < 0.05")
  save_tsv(cbx_fit, file.path(DIR_IMM, "cibersortx_fit_summary.tsv"))
}

message("\n[08_immune_deconvolution.R] Complete.")
message("  Datasets processed: ", paste(names(all_deconv), collapse = ", "))
for (ds_name in names(all_deconv)) {
  message("  ", ds_name, ": ",
          paste(names(all_deconv[[ds_name]]), collapse=", "))
}
message("  Outputs in: ", DIR_IMM)
