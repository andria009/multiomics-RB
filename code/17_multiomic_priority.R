# =============================================================================
# Step 17 — Multi-omic prioritisation: expression evidence (step 11) + methylation (step 15)
# Input : results/07_priority/ (step 11), results/09_multiomic/candidate_methylation_support.tsv
#         (step 15), results/08_methylation/cpg_panel.tsv and cpg_regions.tsv (step 16)
# Output: results/09_multiomic/
#   multiomic_priority.tsv   every gene passing F1-F7 with expression and methylation evidence
#   multiomic_shortlist.tsv  genes in the expression shortlist, the bootstrap core, the
#                            multi-omic shortlist or the exploratory tier, with all evidence
# Multi-omic shortlist (rule fixed before the results were seen): among genes passing F1-F7,
# those with methylation support (step 15) are ranked by the lower 95% bound of |pooled g|;
# the top PRI_N_PER_DIR in each direction form the multi-omic shortlist.
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
if (length(INVASION_EXCLUDE) > 0) stop("Unset INVASION_EXCLUDE.")
load_pkg("dplyr")

pri <- read.delim(file.path(DIR_PRI, "priority_all_candidates.tsv"), stringsAsFactors = FALSE)
cs  <- read.delim(file.path(DIR_MO, "candidate_methylation_support.tsv"), stringsAsFactors = FALSE)
reg <- read.delim(file.path(DIR_METH, "cpg_regions.tsv"), stringsAsFactors = FALSE)

meth_cols <- c("gene", "prom_n_cpg", "prom_mean_beta", "prom_rho_expr", "prom_padj_expr", "prom_rho_pheno",
               "prom_padj_pheno", "body_rho_expr", "body_padj_expr", "meth_testable", "meth_support")
mo <- pri |> dplyr::filter(pass) |> dplyr::left_join(cs[, meth_cols], by = "gene")
mo$meth_support <- mo$meth_support %in% TRUE

# CpG regions (step 16 candidates) annotated to each gene; panel membership
rg <- reg[reg$genes != "" & !is.na(reg$genes), ]
gene_regions <- dplyr::bind_rows(lapply(seq_len(nrow(rg)), function(i)
  data.frame(gene = strsplit(rg$genes[i], ";")[[1]], region = rg$region[i], in_panel = rg$in_panel[i])))
gr <- gene_regions |> dplyr::group_by(gene) |>
  dplyr::summarise(n_candidate_cpg_regions = dplyr::n(), in_cpg_panel = any(in_panel), .groups = "drop")
mo <- dplyr::left_join(mo, gr, by = "gene")
mo$n_candidate_cpg_regions[is.na(mo$n_candidate_cpg_regions)] <- 0L
mo$in_cpg_panel <- mo$in_cpg_panel %in% TRUE

ms <- mo[mo$meth_support, ]
ms <- ms[order(-ms$lower_bound), ]
MO_SHORT <- c(head(ms$gene[ms$g_re > 0], PRI_N_PER_DIR), head(ms$gene[ms$g_re < 0], PRI_N_PER_DIR))
mo$multiomic_shortlist <- mo$gene %in% MO_SHORT
mo <- mo[order(-mo$multiomic_shortlist, -mo$lower_bound), ]
save_tsv(mo, file.path(DIR_MO, "multiomic_priority.tsv"))
message("Genes passing F1-F7: ", nrow(mo), "; testable for methylation: ", sum(mo$meth_testable %in% TRUE),
        "; methylation-supported: ", sum(mo$meth_support), " (", sum(mo$meth_support & mo$g_re > 0), " up, ",
        sum(mo$meth_support & mo$g_re < 0), " down)")
message("Multi-omic shortlist: ", paste(MO_SHORT, collapse = ", "))

adj <- read.delim(file.path(DIR_SUB, "subtype_adjusted_effects.tsv"), stringsAsFactors = FALSE)
adj <- adj[adj$signature == "LIU2021_subtype2" & adj$padj_adjusted < 0.05, ]
EXPLORATORY <- intersect(adj$gene[sign(adj$z_adjusted) == sign(pri$g_re[match(adj$gene, pri$gene)])], mo$gene)
message("Exploratory (subtype-independent) tier: ", paste(EXPLORATORY, collapse = ", "))
keep <- mo$shortlist %in% TRUE | mo$core %in% TRUE | mo$multiomic_shortlist | mo$gene %in% EXPLORATORY
tab <- mo[keep, c("gene", "g_re", "ci_lo", "ci_hi", "lower_bound", "p_top5", "core", "shortlist", "multiomic_shortlist",
                  "prom_n_cpg", "prom_rho_expr", "prom_padj_expr", "prom_rho_pheno", "prom_padj_pheno", "meth_support",
                  "n_candidate_cpg_regions", "in_cpg_panel")]
save_tsv(tab, file.path(DIR_MO, "multiomic_shortlist.tsv"))
print(tab)
message("\n[17_multiomic_priority.R] Complete.")
