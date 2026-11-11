# ==============================================================================
# 00_download_GSE193336.R
# Download and process GSE193336 — external validation dataset
# Bulk RNA-seq: unstimulated and LPS-stimulated human macrophages
#
# GEO: GSE193336
# Design: monocyte-derived macrophages (M-CSF, 5 days), 4 donors
#   U = Unstimulated (donor 1-4)
#   L = LPS-stimulated 50 ng/ml, 24h (donor 1-4)
#
# NOTE: Processed data are normalized counts (non-integer, likely CPM/FPKM)
#   — DESeq2 is NOT appropriate. Use limma-voom in 04_validation_GSE193336.R.
#
# Data: single combined matrix (not per-sample STAR files like GSE130567)
#   GSE193336_GEO_Processed_Data_Raw_Count_Macrophage.csv.gz
#   Columns: geneid, donor.1_unstim, donor.1_LPS, ..., donor.4_LPS
#   Gene IDs: bare Ensembl (no version suffix — no stripping needed)
#
# Outputs:
#   data/raw/GSE193336_counts.csv.gz          — downloaded count matrix
#   data/processed/GSE193336_metadata.csv     — sample metadata
#   data/processed/GSE193336_counts_clean.csv — filtered, column-renamed matrix
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

library(GEOquery)
library(tidyverse)
library(edgeR)

setwd("~/Desktop/Pathway-enrichment-comparison")
dir.create("data/raw",       showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 2. DOWNLOAD COUNT MATRIX
# ==============================================================================

cat("=== Step 1: Downloading GSE193336 count matrix ===\n")

COUNT_URL  <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/",
                     "GSE193336/suppl/",
                     "GSE193336_GEO_Processed_Data_Raw_Count_Macrophage.csv.gz")
COUNT_PATH <- "data/raw/GSE193336_counts.csv.gz"

if (!file.exists(COUNT_PATH)) {
  cat("Downloading count matrix...\n")
  download.file(COUNT_URL, COUNT_PATH, mode = "wb")
  cat("  Done\n")
} else {
  cat("  Count matrix already exists, skipping\n")
}


# ==============================================================================
# 3. LOAD COUNT MATRIX
# ==============================================================================

cat("\n=== Step 2: Loading count matrix ===\n")

counts_raw <- read.csv(COUNT_PATH, row.names = 1)

cat("Dimensions:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")
cat("Column names:\n")
print(colnames(counts_raw))
cat("\nFirst 3 rows:\n")
print(head(counts_raw, 3))

# Check for non-integer values — confirms normalized counts, not raw
cat("\nAre values integers?", all(counts_raw == floor(counts_raw), na.rm = TRUE), "\n")
cat("(If FALSE: normalized counts — use limma-voom, not DESeq2)\n")


# ==============================================================================
# 4. BUILD METADATA
# ==============================================================================

cat("\n=== Step 3: Building metadata ===\n")

# 8 samples: 4 donors x 2 conditions
# Column names in matrix: donor.1_unstim, donor.1_LPS, ..., donor.4_LPS
metadata <- tribble(
  ~gsm_id,      ~col_name,        ~condition, ~donor, ~sample_id,
  "GSM5792562", "donor.1_unstim", "U",        1,      "U_1",
  "GSM5792563", "donor.1_LPS",    "L",        1,      "L_1",
  "GSM5792564", "donor.2_unstim", "U",        2,      "U_2",
  "GSM5792565", "donor.2_LPS",    "L",        2,      "L_2",
  "GSM5792566", "donor.3_unstim", "U",        3,      "U_3",
  "GSM5792567", "donor.3_LPS",    "L",        3,      "L_3",
  "GSM5792568", "donor.4_unstim", "U",        4,      "U_4",
  "GSM5792569", "donor.4_LPS",    "L",        4,      "L_4",
) %>%
  mutate(
    condition = factor(condition, levels = c("U", "L")),
    donor     = factor(donor)
  )

cat("Metadata (8 samples):\n")
print(metadata)
cat("\nCondition breakdown:\n")
print(table(metadata$condition))


# ==============================================================================
# 5. RENAME COLUMNS AND ALIGN TO METADATA
# ==============================================================================

cat("\n=== Step 4: Renaming columns to sample IDs ===\n")

# Verify all expected columns are present
stopifnot(all(metadata$col_name %in% colnames(counts_raw)))

counts_clean <- counts_raw[, metadata$col_name]
colnames(counts_clean) <- metadata$sample_id

cat("Column order after renaming:\n")
print(colnames(counts_clean))


# ==============================================================================
# 6. FILTER LOW-ABUNDANCE GENES
# ==============================================================================
#
# NOTE: filterByExpr expects count-like values. Since these are normalized
# counts (not raw integers), we use a simple rowSums filter instead:
# keep genes with total expression > 0 across all samples.
# A more stringent filter (e.g. rowMeans > 1 in at least one group) is
# applied in the DE analysis script (04_validation_GSE193336.R).

keep <- rowSums(counts_clean) > 0
counts_filtered <- counts_clean[keep, ]

cat("\nAfter removing zero-count genes:", nrow(counts_filtered), "genes retained\n")


# ==============================================================================
# 7. SAVE
# ==============================================================================

cat("\n=== Step 5: Saving processed files ===\n")

write_csv(metadata,                                   "data/processed/GSE193336_metadata.csv")
write.csv(counts_filtered, "data/processed/GSE193336_counts_clean.csv")

cat("Saved:\n")
cat("  data/processed/GSE193336_metadata.csv     —", nrow(metadata), "samples\n")
cat("  data/processed/GSE193336_counts_clean.csv —",
    nrow(counts_filtered), "genes x", ncol(counts_filtered), "samples\n")