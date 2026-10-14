# =============================================================================
# gene_lists.R — curated gene sets used by the analysis
# =============================================================================

# Kapatai-style subgroup markers, used to label GSE29683 tumours (04_invasion_meta_analysis.R).
# These labels are derived from expression, hence the GSE29683 sensitivity analysis.
SUBTYPE_A_MARKERS <- c(  # Photoreceptor-rich (Cluster A)
  "RHO", "SAG", "CNGA1", "CNGB1", "ROM1", "RLBP1", "RDH12",
  "ARR3", "RCVRN", "GUCA1A", "GUCA1B", "PDE6B", "PRPH2"
)

SUBTYPE_C_MARKERS <- c(  # MYCN-driven/invasive (Cluster C)
  "MYCN", "CDK1", "BUB1", "TOP2A", "CCNB1", "KIF11",
  "TTK", "AURKB", "PLK1", "BIRC5"
)

# Residual normal retina in bulk tumours (step 04 GSE97508 QC)
RETINA_PHOTORECEPTOR_MARKERS <- c("RHO", "SAG", "GNAT1", "PDE6A", "NRL", "OPN1LW", "OPN1MW",
                                  "OPN1SW", "ARR3", "GUCA1C", "CNGB3")
RETINA_MULLER_MARKERS <- c("RLBP1", "GLUL", "AQP4", "SLC1A3", "CRYM")

# Cell-type marker sets for scRNA-seq annotation of non-malignant cells (step 06).
# Malignant cells are called separately from inferred copy-number changes.
CELL_TYPE_MARKERS <- list(
  Cone_like     = c("ARR3", "PDE6H", "GNGT2", "GUCA1C", "OPN1SW", "OPN1MW", "CRX", "RXRG", "THRB"),
  Rod           = c("RHO", "NRL", "GNAT1", "PDE6A", "SAG", "CNGA1"),
  Bipolar       = c("VSX2", "OTX2", "PRKCA", "GRM6", "CABP5", "TRPM1"),
  Muller_glia   = c("RLBP1", "GLUL", "AQP4", "SLC1A3", "CRYM", "APOE"),
  Astrocyte     = c("GFAP", "S100B", "PAX2", "ALDH1L1"),
  Microglia_TAM = c("AIF1", "C1QA", "C1QB", "CX3CR1", "P2RY12", "TMEM119", "CD68", "CD163",
                    "MRC1", "CSF1R", "TYROBP"),
  T_NK          = c("CD3D", "CD3E", "CD2", "NKG7", "GNLY", "IL7R"),
  B_plasma      = c("MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"),
  Endothelial   = c("PECAM1", "CDH5", "VWF", "CLDN5", "ESAM", "KDR"),
  Fibroblast    = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRB"),
  RPE           = c("RPE65", "BEST1", "TTR", "PMEL", "TYRP1"),
  Amacrine      = c("GAD1", "GAD2", "TFAP2A", "TFAP2B", "SLC6A9"),
  Horizontal    = c("ONECUT1", "ONECUT2", "LHX1", "PROX1")
)

# Bulk composition and sex proxies (round-2 revision; steps 21 and 24): mean z of each set
COMPOSITION_PROXIES <- list(
  stroma       = CELL_TYPE_MARKERS$Fibroblast,
  endothelial  = CELL_TYPE_MARKERS$Endothelial,
  immune       = c("PTPRC", CELL_TYPE_MARKERS$Microglia_TAM),
  rpe_melanocyte = c(CELL_TYPE_MARKERS$RPE, "TYR", "MLANA", "DCT"),
  rod          = CELL_TYPE_MARKERS$Rod,
  muller       = CELL_TYPE_MARKERS$Muller_glia,
  proliferation = setdiff(SUBTYPE_C_MARKERS, "MYCN")
)
SEX_FEMALE_GENES <- c("XIST", "TSIX")
SEX_MALE_GENES   <- c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y")

# scRNA cluster types that make up the tumour lineage (neuronal / photoreceptor-like);
# glia, RPE, immune and stromal types are excluded (step 07)
TUMOUR_LINEAGE <- c("Cone_like", "Rod", "Bipolar", "Horizontal", "Amacrine")
