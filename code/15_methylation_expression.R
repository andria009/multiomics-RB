# =============================================================================
# Step 15 — Does each gene's methylation track its expression and the phenotype?
# Input : results/08_methylation/ (steps 13-14), results/01_bulk/rnaseq_log2tpm_GSE268136.rds
#         and limma_clusterC_vs_AB_GSE268136.tsv (step 03),
#         results/07_priority/priority_all_candidates.tsv (step 11)
# Output: results/09_multiomic/
#   gene_methylation_expression.tsv    every expressed gene with annotated CpGs: promoter and
#                                      gene-body methylation vs expression and vs phenotype score
#   candidate_methylation_support.tsv  candidate genes (steps 04/11) with the methylation evidence
#   methylation_support_enrichment.tsv candidates vs other expressed genes
#   figures/15_shortlist_methylation_expression
#
# Paired tumours only (methylation + RNA-seq + phenotype score; step 14). Per gene and
# tumour: mean beta of its promoter CpGs (TSS1500, TSS200, 5'UTR, 1st exon) and of its
# gene-body CpGs (body, 3'UTR). Spearman rho with log2(TPM+1) and with the phenotype
# score; BH over genes.
# Methylation support (rule fixed before the results were seen): promoter methylation
#   (i) anti-correlates with expression (rho <= -METH_RHO_MIN, BH < 0.05) and
#   (ii) changes with the phenotype score in the direction that explains the expression
#        change (higher in aggressive tumours for down-regulated genes and vice versa; BH < 0.05).
# =============================================================================

source(here::here("code", "R", "config.R"))
source(here::here("code", "R", "utils.R"))
if (length(INVASION_EXCLUDE) > 0) stop("Unset INVASION_EXCLUDE.")

load_pkg("dplyr")
load_pkg("tidyr")
load_pkg("ggplot2")
load_pkg("Matrix", bioc = FALSE)

beta <- readRDS(file.path(DIR_METH, "meth_beta.rds"))
pg   <- readRDS(file.path(DIR_METH, "meth_probe_gene.rds"))
ps   <- read.delim(file.path(DIR_METH, "meth_phenotype_samples.tsv"), stringsAsFactors = FALSE)
rna  <- readRDS(file.path(DIR_BULK, "rnaseq_log2tpm_GSE268136.rds"))
lim  <- read.delim(file.path(DIR_BULK, "limma_clusterC_vs_AB_GSE268136.tsv"), stringsAsFactors = FALSE)
pri  <- read.delim(file.path(DIR_PRI, "priority_all_candidates.tsv"), stringsAsFactors = FALSE)

B <- beta[, ps$title, drop = FALSE]
E <- rna[, ps$rnaseq_sample, drop = FALSE]
expressed <- intersect(lim$gene, rownames(E))        # genes kept by the step 03 expression filter
message("Paired tumours: ", ncol(B), "; expressed genes: ", length(expressed))

# Gene x tumour mean beta over a CpG set (sparse averaging matrix, no loops)
gene_means <- function(pairs) {
  pairs <- pairs[pairs$gene %in% expressed & pairs$cpg %in% rownames(B), ]
  genes <- sort(unique(pairs$gene))
  A <- Matrix::sparseMatrix(i = match(pairs$gene, genes), j = match(pairs$cpg, rownames(B)), x = 1,
                            dims = c(length(genes), nrow(B)), dimnames = list(genes, rownames(B)))
  n <- Matrix::rowSums(A)
  list(mean = as.matrix(Matrix::Diagonal(x = 1 / n) %*% (A %*% B)), n = n, genes = genes)
}
rank_rows <- function(X) t(apply(X, 1, rank))
row_spearman <- function(X, Y) {           # rho between row i of X and row i of Y (or vector Y)
  rx <- rank_rows(X); ry <- if (is.matrix(Y)) rank_rows(Y) else matrix(rank(Y), nrow(X), length(Y), byrow = TRUE)
  rx <- rx - rowMeans(rx); ry <- ry - rowMeans(ry)
  r <- rowSums(rx * ry) / sqrt(rowSums(rx^2) * rowSums(ry^2))
  n <- ncol(X); tt <- r * sqrt((n - 2) / pmax(1e-12, 1 - r^2))
  data.frame(rho = r, p = 2 * stats::pt(-abs(tt), n - 2))
}
stopifnot(abs(row_spearman(matrix(c(1, 2, 3, 4, 5), 1), c(5, 4, 3, 2, 1))$rho + 1) < 1e-12)

assoc <- function(gm, label) {
  g <- gm$genes; M <- gm$mean
  keep <- apply(M, 1, stats::sd) > 0 & apply(E[g, , drop = FALSE], 1, stats::sd) > 0
  g <- g[keep]; M <- M[keep, , drop = FALSE]
  ex <- row_spearman(M, E[g, , drop = FALSE]); ph <- row_spearman(M, ps$score)
  out <- data.frame(gene = g, n_cpg = gm$n[keep], mean_beta = rowMeans(M),
                    rho_expr = ex$rho, p_expr = ex$p, padj_expr = p.adjust(ex$p, "BH"),
                    rho_pheno = ph$rho, p_pheno = ph$p, padj_pheno = p.adjust(ph$p, "BH"))
  names(out)[-1] <- paste0(label, "_", names(out)[-1])
  out
}
prom <- assoc(gene_means(pg[pg$promoter, ]), "prom")
body <- assoc(gene_means(pg[pg$region %in% c("Body", "3'UTR"), ]), "body")
gx <- dplyr::full_join(prom, body, by = "gene")
gx$expr_rho_pheno <- row_spearman(E[gx$gene, , drop = FALSE], ps$score)$rho
gx$prom_coupled <- (gx$prom_rho_expr <= -METH_RHO_MIN & gx$prom_padj_expr < 0.05) %in% TRUE
save_tsv(gx, file.path(DIR_MO, "gene_methylation_expression.tsv"))
message("Genes with promoter CpGs: ", sum(!is.na(gx$prom_rho_expr)),
        "; promoter methylation anti-correlated with expression (rho <= -", METH_RHO_MIN, ", BH < 0.05): ",
        sum(gx$prom_coupled))

# -----------------------------------------------------------------------------
# Candidates: methylation support
# -----------------------------------------------------------------------------
cand <- pri[, intersect(names(pri), c("gene", "g_re", "pass", "shortlist", "core", "p_top5", "lower_bound"))]
cs <- dplyr::left_join(cand, gx, by = "gene")
dir_up <- cs$g_re > 0
cs$prom_pheno_concordant <- (cs$prom_padj_pheno < 0.05 &
                               ifelse(dir_up, cs$prom_rho_pheno < 0, cs$prom_rho_pheno > 0)) %in% TRUE
cs$meth_support <- cs$prom_coupled & cs$prom_pheno_concordant
cs$meth_testable <- !is.na(cs$prom_rho_expr)
save_tsv(cs, file.path(DIR_MO, "candidate_methylation_support.tsv"))
message("Candidates with promoter CpGs and expression in GSE268136: ", sum(cs$meth_testable), " of ", nrow(cs),
        "; methylation-supported: ", sum(cs$meth_support), " (", sum(cs$meth_support & dir_up), " up, ",
        sum(cs$meth_support & !dir_up), " down)")
message("Among genes passing F1-F7: ", sum(cs$meth_support & cs$pass %in% TRUE), " of ",
        sum(cs$meth_testable & cs$pass %in% TRUE), " testable")
print(cs[cs$shortlist %in% TRUE | cs$core %in% TRUE,
         c("gene", "g_re", "prom_n_cpg", "prom_rho_expr", "prom_padj_expr", "prom_rho_pheno", "prom_padj_pheno", "meth_support")])

# Is promoter methylation-expression coupling more frequent among candidates?
bg <- gx[!is.na(gx$prom_rho_expr), ]
bg$candidate <- bg$gene %in% cand$gene
enr <- dplyr::bind_rows(lapply(list(all = bg, up = bg[!bg$candidate | bg$gene %in% cand$gene[cand$g_re > 0], ],
                                    down = bg[!bg$candidate | bg$gene %in% cand$gene[cand$g_re < 0], ]), function(d) {
  tab <- table(factor(d$candidate, c(TRUE, FALSE)), factor(d$prom_coupled, c(TRUE, FALSE)))
  ft <- stats::fisher.test(tab)
  data.frame(n_candidates = sum(d$candidate), pct_coupled_candidates = 100 * mean(d$prom_coupled[d$candidate]),
             n_other = sum(!d$candidate), pct_coupled_other = 100 * mean(d$prom_coupled[!d$candidate]),
             odds_ratio = unname(ft$estimate), p = ft$p.value)
}), .id = "candidate_set")
save_tsv(enr, file.path(DIR_MO, "methylation_support_enrichment.tsv")); print(enr)

# -----------------------------------------------------------------------------
# Figure: promoter methylation vs expression for the shortlisted and core genes
# -----------------------------------------------------------------------------
pm <- gene_means(pg[pg$promoter, ])$mean
show <- intersect(pri$gene[pri$shortlist %in% TRUE | pri$core %in% TRUE], rownames(pm))
if (length(show)) {
  d <- dplyr::bind_rows(lapply(show, function(g) data.frame(gene = g, beta = pm[g, ], expr = E[g, ], score = ps$score)))
  lab <- setNames(paste0(show, " (rho = ", round(gx$prom_rho_expr[match(show, gx$gene)], 2), ")"), show)
  p <- ggplot2::ggplot(d, ggplot2::aes(beta, expr, colour = score)) + ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_colour_gradient2(low = COL_NONINVASIVE, mid = "grey80", high = COL_INVASIVE, name = "Phenotype\nscore") +
    ggplot2::facet_wrap(~ gene, scales = "free", ncol = 5, labeller = ggplot2::as_labeller(lab)) +
    ggplot2::labs(x = "Mean promoter beta", y = "log2(TPM + 1)") + PLOT_THEME
  save_plot(p, file.path(DIR_FIG, "15_shortlist_methylation_expression"), width = 14, height = 6)
}
message("\n[15_methylation_expression.R] Complete.")
