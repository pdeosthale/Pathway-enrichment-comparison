# ==============================================================================
# 01_differential_expression.R
# Differential expression — GSE193336
# LPS vs Unstimulated, primary human macrophages
#
# Part of: benchmark-pathway-enrichment
# Purpose: Establish DE results and GO BP enrichment used as input
#          for method benchmarking in 02_ora_gsea.R
#
# Design:
#   Dataset  : GSE193336 — monocyte-derived macrophages, 4 donors (paired)
#   Contrast : LPS 50 ng/ml 24h vs Unstimulated
#   Method   : limma-voom with paired donor design (~ donor + condition)
#   Note     : Values are normalized counts (non-integer) — DESeq2 not appropriate
#
# Inputs:
#   data/processed/GSE193336_counts_clean.csv   — from 00_download_GSE193336.R
#   data/processed/GSE193336_metadata.csv       — from 00_download_GSE193336.R
#
# Outputs:
#   results/GSE193336_limma_results.csv         — full DE results, all genes
#   results/GSE193336_GO_LPS_up.rds             — GO BP enrichment, LPS upregulated
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

library(tidyverse)
library(limma)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

setwd("~/Desktop/Pathway-enrichment-comparison")
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)


# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

cat("=== Loading data ===\n")

counts <- read.csv("data/processed/GSE193336_counts_clean.csv",
                   row.names = 1) %>% as.matrix()

meta <- read_csv("data/processed/GSE193336_metadata.csv",
                 show_col_types = FALSE) %>%
  mutate(condition = factor(condition, levels = c("U", "L")),
         donor     = factor(donor))

counts <- counts[, meta$sample_id]
stopifnot(all(colnames(counts) == meta$sample_id))

cat("Counts:", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("Conditions:\n")
print(table(meta$condition))


# ==============================================================================
# 3. EXPRESSION FILTER + LIMMA-VOOM
# ==============================================================================

cat("\n=== Step 1: limma-voom DE (LPS vs Unstimulated) ===\n")

# Keep genes with mean normalized count > 1 in at least one condition
u_mean <- rowMeans(counts[, meta$sample_id[meta$condition == "U"]])
l_mean <- rowMeans(counts[, meta$sample_id[meta$condition == "L"]])
counts_filtered <- counts[u_mean > 1 | l_mean > 1, ]
cat("After expression filter (mean > 1 in either condition):",
    nrow(counts_filtered), "genes\n")

# Paired donor design: ~ donor + condition
# donor term absorbs inter-donor variability; condition captures LPS effect
dge    <- DGEList(counts = counts_filtered) %>% calcNormFactors()
design <- model.matrix(~ donor + condition, data = meta)
cat("Design matrix columns:", colnames(design), "\n")

v   <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design) %>% eBayes()

# Extract LPS vs Unstimulated (conditionL coefficient)
res <- topTable(fit, coef = "conditionL", number = Inf, sort.by = "none") %>%
  rownames_to_column("ensembl_id") %>%
  as_tibble() %>%
  mutate(
    symbol = mapIds(org.Hs.eg.db,
                    keys      = ensembl_id,
                    column    = "SYMBOL",
                    keytype   = "ENSEMBL",
                    multiVals = "first"),
    # Ranking metric for GSEA (used in 02_ora_gsea.R)
    rank_metric = sign(logFC) * -log10(adj.P.Val)
  )

write_csv(res, "results/GSE193336_limma_results.csv")

cat("\nDE summary:\n")
cat("  Total genes tested:", nrow(res), "\n")
cat("  Upregulated (adj.P<0.05, logFC>1):",
    sum(res$adj.P.Val < 0.05 & res$logFC > 1, na.rm = TRUE), "\n")
cat("  Downregulated (adj.P<0.05, logFC< -1):",
    sum(res$adj.P.Val < 0.05 & res$logFC < -1, na.rm = TRUE), "\n")

# Sanity check: canonical LPS-response genes
cat("\nCanonical LPS genes (expected upregulated):\n")
res %>%
  filter(symbol %in% c("TNF", "IL6", "IL1B", "CXCL8", "CCL2", "PTGS2")) %>%
  dplyr::select(symbol, logFC, adj.P.Val) %>%
  arrange(desc(logFC)) %>%
  as.data.frame() %>%
  print()


# ==============================================================================
# 4. GO BP ENRICHMENT — LPS UPREGULATED DEGS
# ==============================================================================
#
# ORA on upregulated DEGs only — used as input for method comparison
# in 02_ora_gsea.R alongside GSEA (which uses the full ranked list).

cat("\n=== Step 2: GO BP enrichment — LPS upregulated DEGs ===\n")

lps_up <- res %>%
  filter(adj.P.Val < 0.05 & logFC > 1) %>%
  pull(ensembl_id)

cat("LPS upregulated DEGs:", length(lps_up), "\n")

bg_entrez <- mapIds(org.Hs.eg.db,
                    keys      = rownames(counts_filtered),
                    column    = "ENTREZID",
                    keytype   = "ENSEMBL",
                    multiVals = "first") %>%
  na.omit() %>% as.character()

lps_up_entrez <- mapIds(org.Hs.eg.db,
                        keys      = lps_up,
                        column    = "ENTREZID",
                        keytype   = "ENSEMBL",
                        multiVals = "first") %>%
  na.omit() %>% as.character()

cat("LPS upregulated DEGs with Entrez IDs:", length(lps_up_entrez), "\n")

go_lps_up <- enrichGO(
  gene          = lps_up_entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

saveRDS(go_lps_up, "results/GSE193336_GO_LPS_up.rds")

cat("Significant GO BP terms (adj.P < 0.05):",
    sum(go_lps_up@result$p.adjust < 0.05), "\n")

cat("\nTop 10 GO BP terms — LPS upregulated:\n")
go_lps_up@result %>%
  filter(p.adjust < 0.05) %>%
  slice_min(p.adjust, n = 10) %>%
  dplyr::select(Description, p.adjust) %>%
  mutate(p.adjust = signif(p.adjust, 3)) %>%
  as.data.frame() %>%
  print()
# Save voom object and design matrix for use in 03_camera.R
saveRDS(v,      "results/voom_object.rds")
saveRDS(design, "results/design_matrix.rds")
cat("\n=== Differential expression complete ===\n")
cat("Outputs:\n")
cat("  results/GSE193336_limma_results.csv  —", nrow(res), "genes\n")
cat("  results/GSE193336_GO_LPS_up.rds      — GO BP enrichment\n")
cat("\nNext step: run 02_ora_gsea.R\n")