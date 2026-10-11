#!/usr/bin/env bash
# Download the raw data from NCBI GEO into data/raw/ and verify checksums.
# Usage (from the data/ folder):  bash download_raw_data.sh
# Only GEO original files are downloaded. The three *.rds files in MANIFEST.md5
# (GSE*.rds, GSE268136_rnaseq.rds) are local GEOquery caches; the analysis code
# recreates them from GEO on first run, so they are skipped here.
set -euo pipefail
cd "$(dirname "$0")"

GEO=https://ftp.ncbi.nlm.nih.gov/geo/series
dir_of() { echo "${1:0:${#1}-3}nnn/$1"; }   # GSE97508 -> GSE97nnn/GSE97508

get() {  # get <subdir> <GSE> <matrix|suppl> <file>; re-downloads a .gz that fails its integrity test
  mkdir -p "raw/$1"
  if [ -f "raw/$1/$4" ]; then
    if [[ "$4" != *.gz ]] || gzip -t "raw/$1/$4" 2>/dev/null; then echo "exists: raw/$1/$4"; return; fi
    echo "incomplete or corrupt, downloading again: raw/$1/$4"; rm -f "raw/$1/$4"
  fi
  echo "downloading $4"
  curl -fL --retry 3 -o "raw/$1/$4.part" "$GEO/$(dir_of "$2")/$3/$4" && mv "raw/$1/$4.part" "raw/$1/$4"
}

# Bulk microarray invasion cohorts (series matrices)
for g in GSE97508 GSE29683 GSE59983 GSE110811; do
  get microarray "$g" matrix "${g}_series_matrix.txt.gz"
done

# Bulk RNA-seq (Ryl et al. 2024): tumour log2(TPM+1) table and sample annotation
get rnaseq GSE268136 suppl  GSE268136_ryl_et_al_rb_tumors_log1ptpm.txt.gz
get rnaseq GSE268136 matrix GSE268136-GPL24676_series_matrix.txt.gz

# DNA methylation (Ryl et al. 2024): processed tumour beta matrix (450K/EPIC common probes)
# and the two platform sub-series matrices with sample annotation
get methylation GSE267015 suppl  GSE267015_ryl_et_al_beta_matrix_tumors.txt.gz
get methylation GSE267015 matrix GSE267015-GPL13534_series_matrix.txt.gz
get methylation GSE267015 matrix GSE267015-GPL21145_series_matrix.txt.gz

# Validation cohort (Liu et al. 2021, Institut Curie; SuperSeries GSE58785): expression
# (HG-U133 Plus 2.0, Brainarray Entrez-gene CDF) and 450K methylation (processed betas);
# histology and subtypes come from Liu et al. Supplementary Data 1 (data/supplementary)
get liu GSE58780 matrix GSE58780_series_matrix.txt.gz
get liu GSE58783 matrix GSE58783_series_matrix.txt.gz
get liu GSE58783 suppl  GSE58783_processed_data_i_methylation_retino_all_samples.csv.gz

# Single-cell RNA-seq (10x MEX files inside RAW.tar)
for g in GSE168434 GSE249995 GSE166173; do
  get scrna "$g" suppl "${g}_RAW.tar"
done

# Verify everything that was downloaded (skips the local .rds caches).
[ -f MANIFEST.md5 ] || { echo "MANIFEST.md5 not found: checksum verification skipped"; exit 0; }
# md5sum on Linux, md5 on macOS.
fail=0
while read -r sum file; do
  [[ "$file" == *.rds ]] && continue
  got=$( (command -v md5sum >/dev/null && md5sum "$file" | cut -d' ' -f1) || md5 -q "$file")
  if [ "$got" = "$sum" ]; then echo "OK    $file"; else echo "FAIL  $file"; fail=1; fi
done < MANIFEST.md5
exit $fail
