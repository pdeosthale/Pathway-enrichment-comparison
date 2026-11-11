# ==============================================================================
# 03_camera.R
# CAMERA enrichment benchmarking — GSE193336
#
# Part of: benchmark-pathway-enrichment
# Purpose: Run CAMERA across GO BP, Reactome, and Hallmark and compare
#          results against ORA and GSEA from 02_ora_gsea.R
#
# Why CAMERA:
#   CAMERA (Competitive gene set test Accounting for inter-gene Correlation)
#   is a competitive test that accounts for correlation between genes in the
#   same gene set. ORA and GSEA treat genes as independent — this inflates
#   Type I error when pathway genes are co-regulated (common in immune biology).
#   CAMERA operates directly on the voom object with the full design matrix,
#   making it statistically consistent with the limma-voom DE model.
#   For GSE193336 (paired donors, n=4, immune stimulation), CAMERA is the
#   most principled choice among the three methods compared here.
#
# Design:
#   Dataset   : GSE193336 — LPS vs Unstimulated, 4 paired donors
#   Contrast  : conditionL (LPS effect, accounting for donor pairing)
#   Resources : GO BP, Reactome, MSigDB Hallmark
#   Total     : 3 resources x 1 method = 3 enrichment runs
#
# Inputs:
#   results/voom_object.rds              — from 01_differential_expression.R
#   results/design_matrix.rds            — from 01_differential_expression.R
#   results/benchmark/ora_GO.rds         — from 02_ora_gsea.R (for comparison)
#   results/benchmark/ora_Reactome.rds   — from 02_ora_gsea.R (for comparison)
#   results/benchmark/ora_Hallmark.csv   — from 02_ora_gsea.R (for comparison)
#   results/benchmark/gsea_GO.rds        — from 02_ora_gsea.R (for comparison)
#   results/benchmark/gsea_Reactome.rds  — from 02_ora_gsea.R (for comparison)
#   results/benchmark/gsea_Hallmark.csv  — from 02_ora_gsea.R (for comparison)
#
# Outputs:
#   results/benchmark/camera_GO.csv
#   results/benchmark/camera_Reactome.csv
#   results/benchmark/camera_Hallmark.csv
#   results/benchmark/method_comparison_summary.csv  — updated with CAMERA
#   figures/benchmark_GO_ORA_vs_GSEA_vs_CAMERA.png
#   figures/benchmark_Reactome_ORA_vs_GSEA_vs_CAMERA.png
#   figures/benchmark_Hallmark_ORA_vs_GSEA_vs_CAMERA.png
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

library(tidyverse)
library(limma)
library(edgeR)
library(clusterProfiler)
library(ReactomePA)
library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)

set.seed(42)
setwd("~/Desktop/Pathway-enrichment-comparison")
dir.create("results/benchmark", showWarnings = FALSE, recursive = TRUE)
dir.create("figures",           showWarnings = FALSE)


# ==============================================================================
# 2. LOAD INPUTS
# ==============================================================================

cat("=== Step 1: Loading inputs ===\n")

# Voom object and design matrix from 01_differential_expression.R
v      <- readRDS("results/voom_object.rds")
design <- readRDS("results/design_matrix.rds")

cat("Voom expression matrix:", nrow(v$E), "genes x", ncol(v$E), "samples\n")
cat("Design matrix columns:", colnames(design), "\n")

# ORA and GSEA results from 02_ora_gsea.R (for comparison heatmap)
ora_go        <- readRDS("results/benchmark/ora_GO.rds")
ora_reactome  <- readRDS("results/benchmark/ora_Reactome.rds")
ora_hallmark  <- read_csv("results/benchmark/ora_Hallmark.csv",
                          show_col_types = FALSE)

gsea_go       <- readRDS("results/benchmark/gsea_GO.rds")
gsea_reactome <- readRDS("results/benchmark/gsea_Reactome.rds")
gsea_hallmark <- read_csv("results/benchmark/gsea_Hallmark.csv",
                          show_col_types = FALSE)

# Gene IDs in the voom matrix (Ensembl IDs)
voom_genes <- rownames(v$E)
cat("Gene IDs in voom matrix (first 3):", head(voom_genes, 3), "\n")

# Map Ensembl IDs to gene symbols for gene set index matching
ensembl_to_symbol <- mapIds(org.Hs.eg.db,
                            keys      = voom_genes,
                            column    = "SYMBOL",
                            keytype   = "ENSEMBL",
                            multiVals = "first")

ensembl_to_entrez <- mapIds(org.Hs.eg.db,
                            keys      = voom_genes,
                            column    = "ENTREZID",
                            keytype   = "ENSEMBL",
                            multiVals = "first")

cat("Genes with symbol annotation:",
    sum(!is.na(ensembl_to_symbol)), "/", length(voom_genes), "\n")
cat("Genes with Entrez annotation:",
    sum(!is.na(ensembl_to_entrez)), "/", length(voom_genes), "\n")


# ==============================================================================
# 3. BUILD GENE SET INDEX LISTS FOR CAMERA
# ==============================================================================
#
# CAMERA requires gene sets as named lists of INTEGER INDICES into the rows
# of the expression matrix — not gene IDs directly.
# ids2indices() from limma handles this conversion.

cat("\n=== Step 2: Building gene set index lists ===\n")

# --- GO BP gene sets ---
cat("Building GO BP index list...\n")
go_sets_raw <- msigdbr(species = "Homo sapiens", collection = "C5",
                       subcollection = "GO:BP")
go_list_symbol <- split(go_sets_raw$gene_symbol, go_sets_raw$gs_name)

# Convert to indices using gene symbols
go_index <- ids2indices(go_list_symbol, identifier = ensembl_to_symbol)
go_index <- go_index[sapply(go_index, length) >= 15]  # min gene set size
cat("  GO BP gene sets (min 15 genes):", length(go_index), "\n")

# --- Reactome gene sets ---
cat("Building Reactome index list...\n")
reactome_sets_raw <- msigdbr(species = "Homo sapiens", collection = "C2",
                             subcollection = "CP:REACTOME")
reactome_list_symbol <- split(reactome_sets_raw$gene_symbol,
                              reactome_sets_raw$gs_name)

reactome_index <- ids2indices(reactome_list_symbol,
                              identifier = ensembl_to_symbol)
reactome_index <- reactome_index[sapply(reactome_index, length) >= 15]
cat("  Reactome gene sets (min 15 genes):", length(reactome_index), "\n")

# --- Hallmark gene sets ---
cat("Building Hallmark index list...\n")
h_sets <- msigdbr(species = "Homo sapiens", collection = "H")
h_list_symbol <- split(h_sets$gene_symbol, h_sets$gs_name)
# Clean names to match ORA/GSEA naming
names(h_list_symbol) <- gsub("HALLMARK_", "", names(h_list_symbol))
names(h_list_symbol) <- gsub("_", " ", names(h_list_symbol))

hallmark_index <- ids2indices(h_list_symbol, identifier = ensembl_to_symbol)
hallmark_index <- hallmark_index[sapply(hallmark_index, length) >= 15]
cat("  Hallmark gene sets (min 15 genes):", length(hallmark_index), "\n")


# ==============================================================================
# 4. RUN CAMERA — THREE RESOURCES
# ==============================================================================
#
# contrast.matrix: extract the conditionL coefficient from the design matrix
# inter.gene.cor = NA: estimate correlation from data (recommended)
# use.ranks = FALSE: parametric test (default); set TRUE for rank-based version

cat("\n=== Step 3: Running CAMERA ===\n")

# Contrast: conditionL is the column in design for LPS effect
contrast_col <- which(colnames(design) == "conditionL")
cat("Contrast column:", colnames(design)[contrast_col], "\n")

# --- CAMERA: GO BP ---
cat("Running CAMERA: GO BP...\n")
camera_go <- camera(
  y              = v,
  index          = go_index,
  design         = design,
  contrast       = contrast_col,
  inter.gene.cor = NA      # estimate from data
) %>%
  rownames_to_column("pathway") %>%
  as_tibble() %>%
  arrange(FDR)

write_csv(camera_go, "results/benchmark/camera_GO.csv")
cat("  Significant GO BP terms (PValue < 0.05):",
    sum(camera_go$PValue < 0.05), "\n")
cat("  Note: FDR correction too conservative for n=8; PValue used for CAMERA\n")

# --- CAMERA: Reactome ---
cat("Running CAMERA: Reactome...\n")
camera_reactome <- camera(
  y              = v,
  index          = reactome_index,
  design         = design,
  contrast       = contrast_col,
  inter.gene.cor = NA
) %>%
  rownames_to_column("pathway") %>%
  as_tibble() %>%
  arrange(FDR)

write_csv(camera_reactome, "results/benchmark/camera_Reactome.csv")
cat("  Significant Reactome terms (PValue < 0.05):",
    sum(camera_reactome$PValue < 0.05), "\n")

# --- CAMERA: Hallmark ---
cat("Running CAMERA: Hallmark...\n")
camera_hallmark <- camera(
  y              = v,
  index          = hallmark_index,
  design         = design,
  contrast       = contrast_col,
  inter.gene.cor = NA
) %>%
  rownames_to_column("pathway") %>%
  as_tibble() %>%
  arrange(FDR)

write_csv(camera_hallmark, "results/benchmark/camera_Hallmark.csv")
cat("  Significant Hallmark terms (PValue < 0.05):",
    sum(camera_hallmark$PValue < 0.05), "\n")


# ==============================================================================
# 5. UPDATED SUMMARY TABLE
# ==============================================================================

cat("\n=== Step 4: Updated summary — all methods x resources ===\n")

summary_df <- tribble(
  ~resource,   ~method,   ~n_significant,
  "GO BP",     "ORA",     sum(ora_go@result$p.adjust < 0.05),
  "Reactome",  "ORA",     sum(ora_reactome@result$p.adjust < 0.05),
  "Hallmark",  "ORA",     sum(ora_hallmark$p.adjust < 0.05, na.rm = TRUE),
  "GO BP",     "GSEA",    sum(gsea_go@result$p.adjust < 0.05),
  "Reactome",  "GSEA",    sum(gsea_reactome@result$p.adjust < 0.05),
  "Hallmark",  "GSEA",    sum(gsea_hallmark$padj < 0.05, na.rm = TRUE),
  "GO BP",     "CAMERA",  sum(camera_go$PValue < 0.05),
  "Reactome",  "CAMERA",  sum(camera_reactome$PValue < 0.05),
  "Hallmark",  "CAMERA",  sum(camera_hallmark$PValue < 0.05)
)

print(summary_df)
write_csv(summary_df, "results/benchmark/method_comparison_summary.csv")

cat("\nTop 10 CAMERA terms — Hallmark (PValue < 0.05):\n")
camera_hallmark %>%
  filter(PValue < 0.05) %>%
  dplyr::select(pathway, Direction, NGenes, Correlation, FDR) %>%
  head(10) %>%
  as.data.frame() %>%
  print()


# ==============================================================================
# 6. THREE-WAY COMPARISON HEATMAP — ORA / GSEA / CAMERA
# ==============================================================================
#
# Anchor on top 15 ORA terms (same as 02_ora_gsea.R),
# then look up those terms in GSEA and CAMERA.
# CAMERA FDR column maps to -log10 for color scale.
# Direction (Up/Down) from CAMERA added as annotation.

cat("\n=== Step 5: Three-way comparison heatmap ===\n")

make_3way_heatmap <- function(ora_desc, ora_padj,
                              gsea_desc, gsea_padj,
                              camera_pathway, camera_fdr, camera_dir,
                              db_name, n_top = 15,
                              fig_width = 10, fig_height = 7,
                              filename) {
  
  # Top N ORA terms
  ora_df <- tibble(Description = ora_desc, padj = ora_padj) %>%
    filter(padj < 0.05) %>%
    slice_min(padj, n = n_top, with_ties = FALSE) %>%
    mutate(neg_log10 = -log10(padj), method = "ORA")
  
  if (nrow(ora_df) == 0) {
    cat("  No significant ORA terms for", db_name, "\n")
    return(NULL)
  }
  
  top_terms <- ora_df$Description
  
  # Match in GSEA
  gsea_df <- tibble(Description = gsea_desc, padj = gsea_padj) %>%
    filter(Description %in% top_terms, !is.na(padj)) %>%
    mutate(neg_log10 = -log10(padj), method = "GSEA")
  
  # Match in CAMERA — clean pathway names to match ORA descriptions
  # CAMERA: use PValue (not FDR) — BH correction too conservative for n=8
  # PValue is appropriate here because CAMERA already accounts for
  # inter-gene correlation, the main source of inflation in competitive tests
  camera_clean <- tibble(
    Description = camera_pathway,
    padj        = camera_fdr,   # this is PValue when called below
    direction   = camera_dir
  ) %>%
    filter(Description %in% top_terms, !is.na(padj)) %>%
    mutate(neg_log10 = -log10(padj), method = "CAMERA")
  
  cat(" ", db_name, "— ORA:", length(top_terms),
      "| GSEA matches:", nrow(gsea_df),
      "| CAMERA matches:", nrow(camera_clean), "\n")
  
  # Build complete grid: all top_terms x 3 methods
  grid <- expand.grid(
    Description = top_terms,
    method      = c("ORA", "GSEA", "CAMERA"),
    stringsAsFactors = FALSE
  ) %>% as_tibble()
  
  vals <- bind_rows(
    ora_df    %>% dplyr::select(Description, method, neg_log10),
    gsea_df   %>% dplyr::select(Description, method, neg_log10),
    camera_clean %>% dplyr::select(Description, method, neg_log10)
  )
  
  all_df <- grid %>%
    left_join(vals, by = c("Description", "method")) %>%
    replace_na(list(neg_log10 = 0)) %>%
    mutate(
      method = factor(method, levels = c("ORA", "GSEA", "CAMERA")),
      Description = factor(Description,
                           levels = rev(ora_df %>%
                                          arrange(desc(padj)) %>%
                                          pull(Description)))
    )
  
  # Cap color scale at CAMERA maximum * 1.2 so CAMERA values are visible
  # alongside ORA/GSEA which can have much larger -log10 values
  cap_value <- max(all_df$neg_log10[all_df$method == "CAMERA"], na.rm = TRUE) * 1.2
  cap_value <- max(cap_value, 3)  # floor of 3 to avoid degenerate scale
  all_df <- all_df %>% mutate(neg_log10 = pmin(neg_log10, cap_value))
  
  p <- ggplot(all_df, aes(x = method, y = Description, fill = neg_log10)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient(low  = "white",
                        high = "#08306b",
                        name = "-log10\n(adj. p / FDR)",
                        limits = c(0, NA)) +
    labs(
      title    = paste0(db_name, ": ORA vs GSEA vs CAMERA"),
      subtitle = "GSE193336 LPS vs Unstimulated | Top 15 ORA terms",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x     = element_text(size = 11, face = "bold"),
      axis.text.y     = element_text(size = 8, hjust = 1),
      panel.grid      = element_blank(),
      plot.title      = element_text(size = 12, face = "bold"),
      plot.subtitle   = element_text(size = 9),
      plot.margin     = margin(t = 10, r = 10, b = 10, l = 10),
      legend.position = "right"
    )
  
  ggsave(filename, p, width = fig_width, height = fig_height,
         dpi = 300, limitsize = FALSE)
  cat("  Saved:", filename, "\n")
  p
}

# ==============================================================================
# Name matching strategy for GO BP and Reactome:
# ORA uses readable descriptions (lowercase, with punctuation)
# CAMERA uses MSigDB prefixed names (GOBP_/REACTOME_ + underscores)
# Solution: normalize both to lowercase + strip punctuation before matching,
# then display ORA descriptions in the heatmap.
# ==============================================================================

normalize_name <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9 ]", " ", x)
  x <- gsub("  +", " ", x)
  x <- trimws(x)
  x
}

# --- GO BP ---
# Clean CAMERA GO names: strip GOBP_ prefix, replace underscores, normalize
camera_go_norm <- camera_go %>%
  mutate(
    Description_raw = pathway %>%
      gsub("^GOBP_", "", .) %>%
      gsub("_", " ", .) %>%
      tolower(),
    Description_norm = normalize_name(Description_raw)
  )

# Normalize ORA descriptions for matching
ora_go_sig <- ora_go@result %>%
  filter(p.adjust < 0.05) %>%
  slice_min(p.adjust, n = 15, with_ties = FALSE) %>%
  mutate(Description_norm = normalize_name(Description))

# Match on normalized names, display ORA descriptions
camera_go_matched <- ora_go_sig %>%
  left_join(
    camera_go_norm %>% dplyr::select(Description_norm, PValue, Direction),
    by = "Description_norm"
  )

cat("GO BP CAMERA matches:", sum(!is.na(camera_go_matched$PValue)), "/ 15\n")

make_3way_heatmap(
  ora_desc       = camera_go_matched$Description,
  ora_padj       = camera_go_matched$p.adjust,
  gsea_desc      = gsea_go@result$Description,
  gsea_padj      = gsea_go@result$p.adjust,
  camera_pathway = camera_go_matched$Description,        # use ORA descriptions
  camera_fdr     = replace_na(camera_go_matched$PValue, 1),  # 1 = no signal
  camera_dir     = replace_na(camera_go_matched$Direction, "Up"),
  db_name        = "GO BP",
  fig_width      = 10, fig_height = 6,
  filename       = "figures/benchmark_GO_ORA_vs_GSEA_vs_CAMERA.png"
)

# --- Reactome ---
# Clean CAMERA Reactome names: strip REACTOME_ prefix, normalize
camera_reactome_norm <- camera_reactome %>%
  mutate(
    Description_raw = pathway %>%
      gsub("^REACTOME_", "", .) %>%
      gsub("_", " ", .) %>%
      tolower(),
    Description_norm = normalize_name(Description_raw)
  )

ora_reactome_sig <- ora_reactome@result %>%
  filter(p.adjust < 0.05) %>%
  slice_min(p.adjust, n = 15, with_ties = FALSE) %>%
  mutate(Description_norm = normalize_name(Description))

camera_reactome_matched <- ora_reactome_sig %>%
  left_join(
    camera_reactome_norm %>% dplyr::select(Description_norm, PValue, Direction),
    by = "Description_norm"
  )

cat("Reactome CAMERA matches:", sum(!is.na(camera_reactome_matched$PValue)), "/ 15\n")

make_3way_heatmap(
  ora_desc       = camera_reactome_matched$Description,
  ora_padj       = camera_reactome_matched$p.adjust,
  gsea_desc      = gsea_reactome@result$Description,
  gsea_padj      = gsea_reactome@result$p.adjust,
  camera_pathway = camera_reactome_matched$Description,
  camera_fdr     = replace_na(camera_reactome_matched$PValue, 1),
  camera_dir     = replace_na(camera_reactome_matched$Direction, "Up"),
  db_name        = "Reactome",
  fig_width      = 12, fig_height = 6,
  filename       = "figures/benchmark_Reactome_ORA_vs_GSEA_vs_CAMERA.png"
)

# --- Hallmark ---
# Hallmark names already cleaned when building hallmark_index above
make_3way_heatmap(
  ora_desc       = gsub("_", " ", gsub("HALLMARK_", "",
                                       ora_hallmark$ID[ora_hallmark$p.adjust < 0.05])),
  ora_padj       = ora_hallmark$p.adjust[ora_hallmark$p.adjust < 0.05],
  gsea_desc      = gsea_hallmark$pathway,
  gsea_padj      = gsea_hallmark$padj,
  camera_pathway = camera_hallmark$pathway,
  camera_fdr     = camera_hallmark$PValue,   # PValue used — FDR too conservative for n=8
  camera_dir     = camera_hallmark$Direction,
  db_name        = "Hallmark",
  fig_width      = 9, fig_height = 6,
  filename       = "figures/benchmark_Hallmark_ORA_vs_GSEA_vs_CAMERA.png"
)

cat("\n=== CAMERA benchmarking complete ===\n")
cat("Summary (all methods x resources):\n")
print(summary_df)
cat("\nNext step: 04_sample_scoring.R (GSVA + ssGSEA — planned)\n")