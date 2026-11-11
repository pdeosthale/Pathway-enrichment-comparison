# ==============================================================================
# 02_ora_gsea.R
# ORA and GSEA enrichment benchmarking — GSE193336
#
# Part of: benchmark-pathway-enrichment
# Purpose: Compare ORA and GSEA across three gene-set resources
#          (GO BP, Reactome, Hallmark) on the same DE results
#
# Design:
#   Dataset  : GSE193336 — LPS vs Unstimulated, 4 paired donors
#   Methods  : ORA (clusterProfiler/ReactomePA) vs GSEA (fgsea)
#   Resources: GO BP, Reactome, MSigDB Hallmark
#   Total    : 2 methods x 3 resources = 6 enrichment runs
#
# Input for ORA  : upregulated DEGs (adj.P < 0.05, logFC > 1)
# Input for GSEA : all genes ranked by sign(logFC) x -log10(adj.P)
#
# Inputs:
#   results/GSE193336_limma_results.csv   — from 01_differential_expression.R
#
# Outputs:
#   results/benchmark/ora_GO.rds
#   results/benchmark/ora_Reactome.rds
#   results/benchmark/ora_Hallmark.csv
#   results/benchmark/gsea_GO.rds
#   results/benchmark/gsea_Reactome.rds
#   results/benchmark/gsea_Hallmark.csv
#   results/benchmark/method_comparison_summary.csv
#   figures/benchmark_GO_ORA_vs_GSEA.png
#   figures/benchmark_Reactome_ORA_vs_GSEA.png
#   figures/benchmark_Hallmark_ORA_vs_GSEA.png
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

library(tidyverse)
library(clusterProfiler)
library(ReactomePA)
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)

set.seed(42)
setwd("~/Desktop/Pathway-enrichment-comparison")
dir.create("results/benchmark", showWarnings = FALSE, recursive = TRUE)
dir.create("figures",           showWarnings = FALSE)


# ==============================================================================
# 2. LOAD DE RESULTS + BUILD INPUTS
# ==============================================================================

cat("=== Step 1: Loading DE results and building inputs ===\n")

res <- read_csv("results/GSE193336_limma_results.csv",
                show_col_types = FALSE) %>%
  filter(!is.na(symbol), !is.na(adj.P.Val))

# Background: all tested genes
bg_entrez <- mapIds(org.Hs.eg.db,
                    keys      = res$ensembl_id,
                    column    = "ENTREZID",
                    keytype   = "ENSEMBL",
                    multiVals = "first") %>%
  na.omit() %>% as.character()

# ORA input: LPS upregulated DEGs
deg_up <- res %>% filter(adj.P.Val < 0.05, logFC > 1)
cat("ORA input — LPS upregulated DEGs:", nrow(deg_up), "\n")

deg_up_entrez <- mapIds(org.Hs.eg.db,
                        keys      = deg_up$ensembl_id,
                        column    = "ENTREZID",
                        keytype   = "ENSEMBL",
                        multiVals = "first") %>%
  na.omit() %>% as.character()
cat("ORA input — with Entrez IDs:", length(deg_up_entrez), "\n")

# GSEA input: all genes ranked by sign(logFC) x -log10(adj.P)
# rank_metric precomputed in 01_differential_expression.R
res_ranked <- res %>%
  filter(!is.na(rank_metric)) %>%
  group_by(symbol) %>%
  slice_max(abs(rank_metric), n = 1) %>%
  ungroup()

ranked_genes <- setNames(res_ranked$rank_metric, res_ranked$symbol)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Entrez-keyed ranked list for Reactome GSEA
entrez_map <- mapIds(org.Hs.eg.db,
                     keys      = res_ranked$ensembl_id,
                     column    = "ENTREZID",
                     keytype   = "ENSEMBL",
                     multiVals = "first")

res_ranked <- res_ranked %>%
  mutate(entrez = entrez_map[ensembl_id]) %>%
  filter(!is.na(entrez))

ranked_entrez <- setNames(res_ranked$rank_metric, res_ranked$entrez)
ranked_entrez <- sort(ranked_entrez, decreasing = TRUE)

cat("GSEA input — ranked genes (symbol):", length(ranked_genes), "\n")
cat("GSEA input — ranked genes (entrez):", length(ranked_entrez), "\n")


# ==============================================================================
# 3. ORA — THREE RESOURCES
# ==============================================================================

cat("\n=== Step 2: ORA — GO BP, Reactome, Hallmark ===\n")

# --- ORA: GO BP ---
cat("Running ORA: GO BP...\n")
ora_go <- enrichGO(
  gene          = deg_up_entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
saveRDS(ora_go, "results/benchmark/ora_GO.rds")
cat("  Significant GO BP terms:", sum(ora_go@result$p.adjust < 0.05), "\n")

# --- ORA: Reactome ---
cat("Running ORA: Reactome...\n")
ora_reactome <- enrichPathway(
  gene          = deg_up_entrez,
  universe      = bg_entrez,
  organism      = "human",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
saveRDS(ora_reactome, "results/benchmark/ora_Reactome.rds")
cat("  Significant Reactome terms:", sum(ora_reactome@result$p.adjust < 0.05), "\n")

# --- ORA: Hallmark ---
cat("Running ORA: Hallmark...\n")
h_sets <- msigdbr(species = "Homo sapiens", collection = "H")
h_t2g  <- h_sets %>% dplyr::select(gs_name, gene_symbol)

ora_hallmark <- enricher(
  gene          = deg_up$symbol,
  TERM2GENE     = h_t2g,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
write_csv(ora_hallmark@result, "results/benchmark/ora_Hallmark.csv")
cat("  Significant Hallmark terms:", sum(ora_hallmark@result$p.adjust < 0.05), "\n")


# ==============================================================================
# 4. GSEA — THREE RESOURCES
# ==============================================================================

cat("\n=== Step 3: GSEA — GO BP, Reactome, Hallmark ===\n")

# --- GSEA: GO BP ---
cat("Running GSEA: GO BP...\n")
gsea_go <- gseGO(
  geneList      = ranked_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  minGSSize     = 15,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)
saveRDS(gsea_go, "results/benchmark/gsea_GO.rds")
cat("  Significant GO BP terms:", sum(gsea_go@result$p.adjust < 0.05), "\n")

# --- GSEA: Reactome ---
cat("Running GSEA: Reactome...\n")
gsea_reactome <- gsePathway(
  geneList      = ranked_entrez,
  organism      = "human",
  minGSSize     = 15,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)
saveRDS(gsea_reactome, "results/benchmark/gsea_Reactome.rds")
cat("  Significant Reactome terms:", sum(gsea_reactome@result$p.adjust < 0.05), "\n")

# --- GSEA: Hallmark ---
cat("Running GSEA: Hallmark...\n")
h_list <- split(h_sets$gene_symbol, h_sets$gs_name)
names(h_list) <- gsub("HALLMARK_", "", names(h_list))
names(h_list) <- gsub("_", " ", names(h_list))

gsea_hallmark <- fgsea(
  pathways    = h_list,
  stats       = ranked_genes,
  minSize     = 15,
  maxSize     = 500,
  nPermSimple = 10000
)
write_csv(gsea_hallmark %>% dplyr::select(-leadingEdge),
          "results/benchmark/gsea_Hallmark.csv")
cat("  Significant Hallmark terms:", sum(gsea_hallmark$padj < 0.05), "\n")


# ==============================================================================
# 5. SUMMARY TABLE
# ==============================================================================

cat("\n=== Step 4: Summary — significant terms per method/resource ===\n")

summary_df <- tribble(
  ~resource,  ~method, ~n_significant,
  "GO BP",    "ORA",   sum(ora_go@result$p.adjust < 0.05),
  "Reactome", "ORA",   sum(ora_reactome@result$p.adjust < 0.05),
  "Hallmark", "ORA",   sum(ora_hallmark@result$p.adjust < 0.05),
  "GO BP",    "GSEA",  sum(gsea_go@result$p.adjust < 0.05),
  "Reactome", "GSEA",  sum(gsea_reactome@result$p.adjust < 0.05),
  "Hallmark", "GSEA",  sum(gsea_hallmark$padj < 0.05)
)

print(summary_df)
write_csv(summary_df, "results/benchmark/method_comparison_summary.csv")


# ==============================================================================
# 6. TOP TERMS PER METHOD — COMPARISON TABLE
# ==============================================================================

cat("\n=== Step 5: Top 10 terms per method/resource ===\n")

get_top <- function(result, method, resource, n = 10) {
  if (inherits(result, "enrichResult") || inherits(result, "gseaResult")) {
    result@result %>%
      filter(p.adjust < 0.05) %>%
      slice_min(p.adjust, n = n) %>%
      mutate(method = method, resource = resource) %>%
      dplyr::select(method, resource, Description, p.adjust)
  } else {
    result %>%
      filter(padj < 0.05) %>%
      slice_min(padj, n = n) %>%
      mutate(method = method, resource = resource,
             Description = pathway, p.adjust = padj) %>%
      dplyr::select(method, resource, Description, p.adjust)
  }
}

top_terms <- bind_rows(
  get_top(ora_go,        "ORA",  "GO BP"),
  get_top(ora_reactome,  "ORA",  "Reactome"),
  get_top(ora_hallmark,  "ORA",  "Hallmark"),
  get_top(gsea_go,       "GSEA", "GO BP"),
  get_top(gsea_reactome, "GSEA", "Reactome"),
  get_top(gsea_hallmark, "GSEA", "Hallmark")
)

cat("\nTop 5 per method/resource:\n")
top_terms %>%
  group_by(method, resource) %>%
  slice_min(p.adjust, n = 5) %>%
  mutate(p.adjust = signif(p.adjust, 3)) %>%
  as.data.frame() %>%
  print()


# ==============================================================================
# 7. COMPARISON HEATMAP — TOP 15 ORA TERMS, SHOW GSEA FOR SAME TERMS
# ==============================================================================

cat("\n=== Step 6: Comparison heatmap (ORA vs GSEA per resource) ===\n")

make_heatmap <- function(ora_desc, ora_padj,
                         gsea_desc, gsea_padj,
                         db_name, n_top = 15,
                         fig_width = 8, fig_height = 7,
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
  
  cat(" ", db_name, "— ORA:", length(top_terms),
      "terms | GSEA matches:", nrow(gsea_df), "\n")
  
  grid <- expand.grid(
    Description = top_terms,
    method      = c("ORA", "GSEA"),
    stringsAsFactors = FALSE
  ) %>% as_tibble()
  
  vals <- bind_rows(
    ora_df  %>% dplyr::select(Description, method, neg_log10),
    gsea_df %>% dplyr::select(Description, method, neg_log10)
  )
  
  all_df <- grid %>%
    left_join(vals, by = c("Description", "method")) %>%
    replace_na(list(neg_log10 = 0)) %>%
    mutate(
      method = factor(method, levels = c("ORA", "GSEA")),
      Description = factor(Description,
                           levels = rev(ora_df %>%
                                          arrange(desc(padj)) %>%
                                          pull(Description)))
    )
  
  p <- ggplot(all_df, aes(x = method, y = Description, fill = neg_log10)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient(low  = "white",
                        high = "#08306b",
                        name = "-log10\n(adj. p)",
                        limits = c(0, NA)) +
    labs(
      title    = paste0(db_name, ": ORA vs GSEA"),
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

make_heatmap(
  ora_desc  = ora_go@result$Description,
  ora_padj  = ora_go@result$p.adjust,
  gsea_desc = gsea_go@result$Description,
  gsea_padj = gsea_go@result$p.adjust,
  db_name   = "GO BP",
  fig_width = 11, fig_height = 6,
  filename  = "figures/benchmark_GO_ORA_vs_GSEA.png"
)

make_heatmap(
  ora_desc  = ora_reactome@result$Description,
  ora_padj  = ora_reactome@result$p.adjust,
  gsea_desc = gsea_reactome@result$Description,
  gsea_padj = gsea_reactome@result$p.adjust,
  db_name   = "Reactome",
  fig_width = 11, fig_height = 6,
  filename  = "figures/benchmark_Reactome_ORA_vs_GSEA.png"
)

h_ora_desc <- gsub("_", " ", gsub("HALLMARK_", "", ora_hallmark@result$ID))
h_ora_padj <- ora_hallmark@result$p.adjust

make_heatmap(
  ora_desc  = h_ora_desc,
  ora_padj  = h_ora_padj,
  gsea_desc = gsea_hallmark$pathway,
  gsea_padj = gsea_hallmark$padj,
  db_name   = "Hallmark",
  fig_width = 8, fig_height = 6,
  filename  = "figures/benchmark_Hallmark_ORA_vs_GSEA.png"
)

cat("\n=== Benchmarking complete ===\n")
cat("Summary:\n")
print(summary_df)
cat("\nNext step: 03_camera.R (planned)\n")