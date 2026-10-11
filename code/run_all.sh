#!/usr/bin/env bash
# Run the full multi-omic (MO) analysis in order. Stops at the first failing step.
# Usage (from the project root):  bash code/run_all.sh
# Uses the first Rscript on PATH (CRAN R 4.6.1 was used for the paper).
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs
LOG=logs/run_all.log
: > "$LOG"

run() {  # run <script> [ENV=value]
  echo "=== $(date '+%F %T')  $1 ${2:-}" | tee -a "$LOG"
  env ${2:-} Rscript "code/$1" >> "$LOG" 2>&1 || { echo "FAILED: $1 (see $LOG)"; exit 1; }
}

run 01_download_geo.R
run 02_normalize_microarray.R
run 03_rnaseq_gse268136.R
run 04_invasion_meta_analysis.R
# Sensitivity: without the expression-labelled cohorts (needed by step 11)
run 04_invasion_meta_analysis.R INVASION_EXCLUDE=GSE29683
run 04_invasion_meta_analysis.R INVASION_EXCLUDE=GSE29683,GSE59983
run 05_invasion_signature.R
run 06_scrna_integration.R
run 06b_copykat_validation.R
run 07_scrna_validation.R
run 08_immune_deconvolution.R
run 09_invasion_immune.R
run 10_subtype_benchmark.R
run 11_priority_ranking.R
run 13_methylation_preprocess.R
run 14_methylation_phenotype.R
run 15_methylation_expression.R
run 16_cpg_marker_panel.R
run 17_multiomic_priority.R
run 18_liu_validation.R
run 19_supporting_info.R
echo "=== $(date '+%F %T')  done" | tee -a "$LOG"
