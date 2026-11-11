# ==============================================================================
# 05_additional_resources.R
# ORA, GSEA, and CAMERA with C7 Immunologic Signatures and WikiPathways
#
# Part of: benchmark-pathway-enrichment
# Purpose: Extend the method comparison to two additional gene-set resources
#          with different design philosophy from GO/Reactome/Hallmark:
#
#   C7 Immunologic Signatures (MSigDB C7):
#     ~5,219 gene sets derived from immune cell perturbation experiments.
#     Highest specificity for immune biology — designed for datasets like
#     GSE193336 where the perturbation is immunological (LPS stimulation).
#     Expected to outperform general resources for recovering LPS macrophage
#     programs.
#
#   WikiPathways:
#     ~925 community-curated pathways, fully open access, frequently updated.
#     Complementary to Reactome — similar mechanistic focus but broader
#     contributor base and faster annotation of emerging biology.
#
# Design:
#   Dataset   : GSE193336 — LPS vs Unstimulated, 4 paired donors
#   Methods   : ORA, GSEA, CAMERA
#   Resources : C7 Immunologic Signatures, WikiPathways
#   Total     : 3 methods x 2 resources = 6 enrichment runs
#
# Inputs:
#   results/GSE193336_limma_results.csv   — from 01_differential_expression.R
#   results/voom_object.rds               — from 01_differential_expression.R
#   results/design_matrix.rds             — from 01_differential_expression.R
#
# Outputs:
#   results/benchmark/ora_C7.csv
#   results/benchmark/ora_WikiPathways.csv
#   results/benchmark/gsea_C7.csv
#   results/benchmark/gsea_WikiPathways.csv
#   results/benchmark/camera_C7.csv
#   results/benchmark/camera_WikiPathways.csv
#   results/benchmark/method_comparison_summary_full.csv  — all 5 resources
#   figures/benchmark_C7_ORA_vs_GSEA_vs_CAMERA.png
#   figures/benchmark_WikiPathways_ORA_vs_GSEA_vs_CAMERA.png
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

library(tidyverse)
library(limma)
library(clusterProfiler)
library(fgsea)
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

res <- read_csv("results/GSE193336_limma_results.csv",
                show_col_types = FALSE) %>%
  filter(!is.na(symbol), !is.na(adj.P.Val))

v      <- readRDS("results/voom_object.rds")
design <- readRDS("results/design_matrix.rds")

# Background Entrez IDs
bg_entrez <- mapIds(org.Hs.eg.db,
                    keys      = res$ensembl_id,
                    column    = "ENTREZID",
                    keytype   = "ENSEMBL",
                    multiVals = "first") %>%
  na.omit() %>% as.character()

# ORA input: LPS upregulated DEGs
deg_up <- res %>% filter(adj.P.Val < 0.05, logFC > 1)
cat("ORA input — upregulated DEGs:", nrow(deg_up), "\n")

# GSEA ranked list (symbol-keyed)
res_ranked <- res %>%
  filter(!is.na(rank_metric)) %>%
  group_by(symbol) %>%
  slice_max(abs(rank_metric), n = 1) %>%
  ungroup()

ranked_genes <- setNames(res_ranked$rank_metric, res_ranked$symbol)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Ensembl-to-symbol map for CAMERA
voom_genes         <- rownames(v$E)
ensembl_to_symbol  <- mapIds(org.Hs.eg.db,
                             keys      = voom_genes,
                             column    = "SYMBOL",
                             keytype   = "ENSEMBL",
                             multiVals = "first")

contrast_col <- which(colnames(design) == "conditionL")

cat("GSEA input — ranked genes:", length(ranked_genes), "\n")
cat("CAMERA contrast column:", colnames(design)[contrast_col], "\n")


# ==============================================================================
# 3. LOAD GENE SET RESOURCES
# ==============================================================================

cat("\n=== Step 2: Loading gene set resources ===\n")

# --- C7 Immunologic Signatures ---
cat("Loading C7...\n")
c7_sets <- msigdbr(species = "Homo sapiens", collection = "C7")
cat("  Total C7 gene sets:", n_distinct(c7_sets$gs_name), "\n")

c7_t2g   <- c7_sets %>% dplyr::select(gs_name, gene_symbol)
c7_list  <- split(c7_sets$gene_symbol, c7_sets$gs_name)

c7_index <- ids2indices(c7_list, identifier = ensembl_to_symbol)
c7_index <- c7_index[sapply(c7_index, length) >= 15 &
                       sapply(c7_index, length) <= 500]
cat("  C7 gene sets for CAMERA (15-500 genes):", length(c7_index), "\n")

# --- WikiPathways ---
cat("Loading WikiPathways...\n")
wp_sets <- msigdbr(species = "Homo sapiens",
                   collection = "C2", subcollection = "CP:WIKIPATHWAYS")
cat("  Total WikiPathways gene sets:", n_distinct(wp_sets$gs_name), "\n")

wp_t2g  <- wp_sets %>% dplyr::select(gs_name, gene_symbol)
wp_list <- split(wp_sets$gene_symbol, wp_sets$gs_name)
# Clean names: strip WP prefix and underscores
names(wp_list) <- names(wp_list) %>%
  gsub("^WIKIPATHWAYS_", "", .) %>%
  gsub("_", " ", .)

wp_index <- ids2indices(wp_list, identifier = ensembl_to_symbol)
wp_index <- wp_index[sapply(wp_index, length) >= 15 &
                       sapply(wp_index, length) <= 500]
cat("  WikiPathways for CAMERA (15-500 genes):", length(wp_index), "\n")


# ==============================================================================
# 4. ORA — C7 AND WIKIPATHWAYS
# ==============================================================================

cat("\n=== Step 3: ORA ===\n")

# --- ORA: C7 ---
cat("Running ORA: C7...\n")
ora_c7 <- enricher(
  gene          = deg_up$symbol,
  TERM2GENE     = c7_t2g,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
write_csv(ora_c7@result, "results/benchmark/ora_C7.csv")
cat("  Significant C7 terms (FDR < 0.05):",
    sum(ora_c7@result$p.adjust < 0.05), "\n")

# --- ORA: WikiPathways ---
cat("Running ORA: WikiPathways...\n")
wp_t2g_clean <- wp_sets %>%
  mutate(gs_name = gs_name %>%
           gsub("^WIKIPATHWAYS_", "", .) %>%
           gsub("_", " ", .)) %>%
  dplyr::select(gs_name, gene_symbol)

ora_wp <- enricher(
  gene          = deg_up$symbol,
  TERM2GENE     = wp_t2g_clean,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
write_csv(ora_wp@result, "results/benchmark/ora_WikiPathways.csv")
cat("  Significant WikiPathways terms (FDR < 0.05):",
    sum(ora_wp@result$p.adjust < 0.05), "\n")


# ==============================================================================
# 5. GSEA — C7 AND WIKIPATHWAYS
# ==============================================================================

cat("\n=== Step 4: GSEA ===\n")

# --- GSEA: C7 ---
cat("Running GSEA: C7 (5219 sets — may take several minutes)...\n")
gsea_c7 <- fgsea(
  pathways    = c7_list,
  stats       = ranked_genes,
  minSize     = 15,
  maxSize     = 500,
  nPermSimple = 10000
)
write_csv(gsea_c7 %>% dplyr::select(-leadingEdge),
          "results/benchmark/gsea_C7.csv")
cat("  Significant C7 terms (FDR < 0.05):",
    sum(gsea_c7$padj < 0.05, na.rm = TRUE), "\n")

# --- GSEA: WikiPathways ---
cat("Running GSEA: WikiPathways...\n")
gsea_wp <- fgsea(
  pathways    = wp_list,
  stats       = ranked_genes,
  minSize     = 15,
  maxSize     = 500,
  nPermSimple = 10000
)
write_csv(gsea_wp %>% dplyr::select(-leadingEdge),
          "results/benchmark/gsea_WikiPathways.csv")
cat("  Significant WikiPathways terms (FDR < 0.05):",
    sum(gsea_wp$padj < 0.05, na.rm = TRUE), "\n")


# ==============================================================================
# 6. CAMERA — C7 AND WIKIPATHWAYS
# ==============================================================================

cat("\n=== Step 5: CAMERA ===\n")

# --- CAMERA: C7 ---
cat("Running CAMERA: C7 (may take several minutes)...\n")
camera_c7 <- camera(
  y              = v,
  index          = c7_index,
  design         = design,
  contrast       = contrast_col,
  inter.gene.cor = NA
) %>%
  rownames_to_column("pathway") %>%
  as_tibble() %>%
  arrange(PValue)

write_csv(camera_c7, "results/benchmark/camera_C7.csv")
cat("  Significant C7 terms (PValue < 0.05):",
    sum(camera_c7$PValue < 0.05), "\n")

# --- CAMERA: WikiPathways ---
cat("Running CAMERA: WikiPathways...\n")
camera_wp <- camera(
  y              = v,
  index          = wp_index,
  design         = design,
  contrast       = contrast_col,
  inter.gene.cor = NA
) %>%
  rownames_to_column("pathway") %>%
  as_tibble() %>%
  arrange(PValue)

write_csv(camera_wp, "results/benchmark/camera_WikiPathways.csv")
cat("  Significant WikiPathways terms (PValue < 0.05):",
    sum(camera_wp$PValue < 0.05), "\n")


# ==============================================================================
# 7. FULL SUMMARY TABLE — ALL 5 RESOURCES
# ==============================================================================

cat("\n=== Step 6: Full summary — all methods x resources ===\n")

# Load existing results from 02 and 03
ora_go       <- readRDS("results/benchmark/ora_GO.rds")
ora_reactome <- readRDS("results/benchmark/ora_Reactome.rds")
ora_hallmark <- read_csv("results/benchmark/ora_Hallmark.csv",
                         show_col_types = FALSE)
gsea_go       <- readRDS("results/benchmark/gsea_GO.rds")
gsea_reactome <- readRDS("results/benchmark/gsea_Reactome.rds")
gsea_hallmark <- read_csv("results/benchmark/gsea_Hallmark.csv",
                          show_col_types = FALSE)
camera_go       <- read_csv("results/benchmark/camera_GO.csv",
                            show_col_types = FALSE)
camera_reactome <- read_csv("results/benchmark/camera_Reactome.csv",
                            show_col_types = FALSE)
camera_hallmark <- read_csv("results/benchmark/camera_Hallmark.csv",
                            show_col_types = FALSE)

summary_full <- tribble(
  ~resource,        ~method,   ~n_significant,
  "GO BP",          "ORA",     sum(ora_go@result$p.adjust < 0.05),
  "Reactome",       "ORA",     sum(ora_reactome@result$p.adjust < 0.05),
  "Hallmark",       "ORA",     sum(ora_hallmark$p.adjust < 0.05, na.rm=TRUE),
  "C7",             "ORA",     sum(ora_c7@result$p.adjust < 0.05),
  "WikiPathways",   "ORA",     sum(ora_wp@result$p.adjust < 0.05),
  "GO BP",          "GSEA",    sum(gsea_go@result$p.adjust < 0.05),
  "Reactome",       "GSEA",    sum(gsea_reactome@result$p.adjust < 0.05),
  "Hallmark",       "GSEA",    sum(gsea_hallmark$padj < 0.05, na.rm=TRUE),
  "C7",             "GSEA",    sum(gsea_c7$padj < 0.05, na.rm=TRUE),
  "WikiPathways",   "GSEA",    sum(gsea_wp$padj < 0.05, na.rm=TRUE),
  "GO BP",          "CAMERA",  sum(camera_go$PValue < 0.05),
  "Reactome",       "CAMERA",  sum(camera_reactome$PValue < 0.05),
  "Hallmark",       "CAMERA",  sum(camera_hallmark$PValue < 0.05),
  "C7",             "CAMERA",  sum(camera_c7$PValue < 0.05),
  "WikiPathways",   "CAMERA",  sum(camera_wp$PValue < 0.05)
) %>%
  mutate(resource = factor(resource,
                           levels = c("GO BP","Reactome","Hallmark",
                                      "C7","WikiPathways")),
         method   = factor(method, levels = c("ORA","GSEA","CAMERA")))

print(summary_full)
write_csv(summary_full, "results/benchmark/method_comparison_summary_full.csv")


# ==============================================================================
# 8. THREE-WAY HEATMAPS — C7 AND WIKIPATHWAYS
# ==============================================================================

cat("\n=== Step 7: Comparison heatmaps ===\n")

normalize_name <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9 ]", " ", x)
  x <- gsub("  +", " ", x)
  x <- trimws(x)
  x
}

make_3way_heatmap <- function(ora_desc, ora_padj,
                              gsea_desc, gsea_padj,
                              camera_pathway, camera_pvalue,
                              db_name, n_top = 15,
                              fig_width = 10, fig_height = 7,
                              filename) {
  
  ora_df <- tibble(Description = ora_desc, padj = ora_padj) %>%
    filter(padj < 0.05) %>%
    slice_min(padj, n = n_top, with_ties = FALSE) %>%
    mutate(neg_log10 = -log10(padj), method = "ORA")
  
  if (nrow(ora_df) == 0) {
    cat("  No significant ORA terms for", db_name, "\n")
    return(NULL)
  }
  
  top_terms <- ora_df$Description
  
  gsea_df <- tibble(Description = gsea_desc, padj = gsea_padj) %>%
    filter(Description %in% top_terms, !is.na(padj)) %>%
    mutate(neg_log10 = -log10(padj), method = "GSEA")
  
  camera_df <- tibble(Description = camera_pathway,
                      padj        = camera_pvalue) %>%
    filter(Description %in% top_terms, !is.na(padj)) %>%
    mutate(neg_log10 = -log10(padj), method = "CAMERA")
  
  cat(" ", db_name, "— ORA:", length(top_terms),
      "| GSEA:", nrow(gsea_df),
      "| CAMERA:", nrow(camera_df), "\n")
  
  grid <- expand.grid(Description = top_terms,
                      method = c("ORA","GSEA","CAMERA"),
                      stringsAsFactors = FALSE) %>% as_tibble()
  
  vals <- bind_rows(
    ora_df    %>% dplyr::select(Description, method, neg_log10),
    gsea_df   %>% dplyr::select(Description, method, neg_log10),
    camera_df %>% dplyr::select(Description, method, neg_log10)
  )
  
  all_df <- grid %>%
    left_join(vals, by = c("Description","method")) %>%
    replace_na(list(neg_log10 = 0)) %>%
    mutate(
      method = factor(method, levels = c("ORA","GSEA","CAMERA")),
      Description = factor(Description,
                           levels = rev(ora_df %>%
                                          arrange(desc(padj)) %>%
                                          pull(Description)))
    )
  
  cap_value <- max(all_df$neg_log10[all_df$method == "CAMERA"],
                   na.rm = TRUE) * 1.2
  cap_value <- max(cap_value, 3)
  all_df <- all_df %>% mutate(neg_log10 = pmin(neg_log10, cap_value))
  
  p <- ggplot(all_df, aes(x = method, y = Description, fill = neg_log10)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient(low  = "white", high = "#08306b",
                        name = "-log10\n(adj. p / p)",
                        limits = c(0, NA)) +
    labs(title    = paste0(db_name, ": ORA vs GSEA vs CAMERA"),
         subtitle = "GSE193336 LPS vs Unstimulated | Top 15 ORA terms",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x     = element_text(size = 11, face = "bold"),
      axis.text.y     = element_text(size = 8, hjust = 1),
      panel.grid      = element_blank(),
      plot.title      = element_text(size = 12, face = "bold"),
      plot.subtitle   = element_text(size = 9),
      plot.margin     = margin(t=10, r=10, b=10, l=10),
      legend.position = "right"
    )
  
  ggsave(filename, p, width = fig_width, height = fig_height,
         dpi = 300, limitsize = FALSE)
  cat("  Saved:", filename, "\n")
  p
}

# --- C7 ---
# C7 names kept as raw MSigDB strings (GSE##### or AUTHORNAME_ format)
# These are long but informative — they encode source experiment and direction
c7_ora_sig <- ora_c7@result %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust) %>%
  mutate(Description_norm = normalize_name(ID))

camera_c7_norm <- camera_c7 %>%
  mutate(Description_norm = normalize_name(pathway))

c7_matched <- c7_ora_sig %>%
  left_join(camera_c7_norm %>%
              dplyr::select(Description_norm, PValue, Direction),
            by = "Description_norm")

cat("C7 CAMERA matches:", sum(!is.na(c7_matched$PValue)), "/",
    nrow(c7_ora_sig), "\n")

make_3way_heatmap(
  ora_desc       = c7_matched$ID,
  ora_padj       = c7_matched$p.adjust,
  gsea_desc      = gsea_c7$pathway,
  gsea_padj      = gsea_c7$padj,
  camera_pathway = c7_matched$ID,
  camera_pvalue  = replace_na(c7_matched$PValue, 1),
  db_name        = "C7 Immunologic Signatures",
  fig_width      = 13, fig_height = 7,
  filename       = "figures/benchmark_C7_ORA_vs_GSEA_vs_CAMERA.png"
)

# --- WikiPathways ---
wp_ora_sig <- ora_wp@result %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust)

camera_wp_norm <- camera_wp %>%
  mutate(Description_norm = normalize_name(pathway))

wp_ora_sig <- wp_ora_sig %>%
  mutate(Description_norm = normalize_name(ID))

wp_matched <- wp_ora_sig %>%
  left_join(camera_wp_norm %>%
              dplyr::select(Description_norm, PValue, Direction),
            by = "Description_norm")

cat("WikiPathways CAMERA matches:", sum(!is.na(wp_matched$PValue)), "/",
    nrow(wp_ora_sig), "\n")

make_3way_heatmap(
  ora_desc       = wp_matched$ID,
  ora_padj       = wp_matched$p.adjust,
  gsea_desc      = gsea_wp$pathway,
  gsea_padj      = gsea_wp$padj,
  camera_pathway = wp_matched$ID,
  camera_pvalue  = replace_na(wp_matched$PValue, 1),
  db_name        = "WikiPathways",
  fig_width      = 10, fig_height = 7,
  filename       = "figures/benchmark_WikiPathways_ORA_vs_GSEA_vs_CAMERA.png"
)

cat("\n=== Additional resources complete ===\n")
cat("Full summary (all 5 resources x 3 methods):\n")
print(summary_full %>% arrange(resource, method), n = 15)