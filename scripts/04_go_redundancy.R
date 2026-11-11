# ==============================================================================
# 04_go_redundancy.R
# GO BP redundancy quantification via semantic similarity clustering
#
# Part of: benchmark-pathway-enrichment
# Purpose: Quantify how many of the 1,007 significant GO BP ORA terms
#          represent genuinely distinct biological themes vs redundant
#          annotations of the same biology.
#
# Method:
#   1. Extract significant GO BP terms from ORA results
#   2. Calculate pairwise semantic similarity (Wang method, GOSemSim)
#   3. Cluster terms using rrvgo::reduceSimMatrix()
#   4. Report: X terms → Y unique biological themes
#   5. Visualize: treemap of clustered terms + scatter plot of semantic space
#
# Why this matters:
#   GO BP has ~7,500 terms with hierarchical structure — child terms inherit
#   parent annotations, creating redundancy. "regulation of inflammatory
#   response", "positive regulation of inflammatory response", and
#   "positive regulation of defense response" all capture overlapping biology.
#   Semantic similarity clustering collapses these into representative themes,
#   giving a cleaner picture of the actual biological content.
#
# Inputs:
#   results/benchmark/ora_GO.rds   — from 02_ora_gsea.R
#
# Outputs:
#   results/benchmark/go_similarity_matrix.rds   — pairwise similarity matrix
#   results/benchmark/go_reduced_terms.csv       — representative terms per cluster
#   results/benchmark/go_redundancy_summary.csv  — before/after counts
#   figures/go_treemap.png                       — treemap of clustered terms
#   figures/go_scatter.png                       — semantic space scatter plot
# ==============================================================================


# ==============================================================================
# 1. SETUP
# ==============================================================================

library(tidyverse)
library(GOSemSim)
library(rrvgo)
library(org.Hs.eg.db)

set.seed(42)
setwd("~/Desktop/Pathway-enrichment-comparison")
dir.create("results/benchmark", showWarnings = FALSE, recursive = TRUE)
dir.create("figures",           showWarnings = FALSE)


# ==============================================================================
# 2. LOAD GO BP ORA RESULTS
# ==============================================================================

cat("=== Step 1: Loading GO BP ORA results ===\n")

ora_go <- readRDS("results/benchmark/ora_GO.rds")

# Extract significant terms
sig_terms <- ora_go@result %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust)

cat("Significant GO BP terms (FDR < 0.05):", nrow(sig_terms), "\n")
cat("GO IDs (first 5):", head(sig_terms$ID, 5), "\n")


# ==============================================================================
# 3. BUILD GOSEM SIM DATA OBJECT
# ==============================================================================

cat("\n=== Step 2: Building GOSemSim data object ===\n")
cat("Computing information content for GO BP terms (may take ~60 seconds)...\n")

godata_bp <- godata("org.Hs.eg.db", ont = "BP", computeIC = TRUE)
cat("Done\n")


# ==============================================================================
# 4. CALCULATE PAIRWISE SEMANTIC SIMILARITY
# ==============================================================================
#
# Wang method: accounts for the contribution of all ancestor terms,
# weighted by their semantic contribution to each term.
# Most appropriate for comparing terms at different levels of GO hierarchy.
# Alternative: "Rel" (relevance) — similar results, slightly faster.

cat("\n=== Step 3: Calculating pairwise semantic similarity ===\n")
cat("Terms to compare:", nrow(sig_terms), "\n")
cat("Computing similarity matrix (may take several minutes for >500 terms)...\n")

go_ids <- sig_terms$ID

sim_matrix <- mgoSim(
  GO1     = go_ids,
  GO2     = go_ids,
  semData = godata_bp,
  measure = "Wang",
  combine = NULL       # return full matrix, not single score
)

cat("Similarity matrix dimensions:", dim(sim_matrix), "\n")
cat("Value range:", round(min(sim_matrix, na.rm=TRUE), 3),
    "to", round(max(sim_matrix, na.rm=TRUE), 3), "\n")

saveRDS(sim_matrix, "results/benchmark/go_similarity_matrix.rds")
cat("Saved: results/benchmark/go_similarity_matrix.rds\n")


# ==============================================================================
# 5. REDUCE REDUNDANT TERMS VIA RRVGO
# ==============================================================================
#
# reduceSimMatrix() clusters GO terms by semantic similarity.
# threshold: terms with similarity > threshold are collapsed into one cluster.
# 0.7 is standard (used in published benchmarks); higher = fewer clusters.
# Representative term per cluster = term with highest score (here: -log10 padj)

cat("\n=== Step 4: Reducing redundant terms ===\n")

scores <- setNames(-log10(sig_terms$p.adjust), sig_terms$ID)

reduced <- reduceSimMatrix(
  simMatrix  = sim_matrix,
  scores     = scores,
  threshold  = 0.7,
  orgdb      = "org.Hs.eg.db"
)

cat("Before reduction:", nrow(sig_terms), "terms\n")
cat("After reduction: ", length(unique(reduced$parentTerm)), "unique themes\n")
cat("Redundancy:      ",
    round(100 * (1 - length(unique(reduced$parentTerm)) / nrow(sig_terms)), 1),
    "%\n")

# Add p.adjust from original results
reduced <- reduced %>%
  left_join(sig_terms %>% dplyr::select(ID, p.adjust, Description),
            by = c("go" = "ID"))

write_csv(reduced, "results/benchmark/go_reduced_terms.csv")
cat("Saved: results/benchmark/go_reduced_terms.csv\n")

# Summary table
redundancy_summary <- tibble(
  resource        = "GO BP",
  method          = "ORA",
  n_sig_terms     = nrow(sig_terms),
  n_unique_themes = length(unique(reduced$parentTerm)),
  redundancy_pct  = round(100 * (1 - length(unique(reduced$parentTerm)) /
                                   nrow(sig_terms)), 1)
)

print(redundancy_summary)
write_csv(redundancy_summary, "results/benchmark/go_redundancy_summary.csv")

cat("\nTop 20 representative themes (by significance):\n")
reduced %>%
  group_by(parentTerm) %>%
  slice_min(p.adjust, n = 1) %>%
  ungroup() %>%
  arrange(p.adjust) %>%
  head(20) %>%
  dplyr::select(parentTerm, p.adjust) %>%
  mutate(p.adjust = signif(p.adjust, 3)) %>%
  as.data.frame() %>%
  print()


# ==============================================================================
# 6. BAR CHART — TOP REPRESENTATIVE THEMES + REDUNDANCY SUMMARY
# ==============================================================================

cat("\n=== Step 5: Bar chart of representative themes ===\n")

n_terms     <- nrow(sig_terms)
n_themes    <- length(unique(reduced$parentTerm))
redundancy  <- round(100 * (1 - n_themes / n_terms), 1)

# Top 25 representative themes by significance
top_themes <- reduced %>%
  group_by(parentTerm) %>%
  summarise(
    p.adjust    = min(p.adjust, na.rm = TRUE),
    n_terms     = n(),           # how many GO terms this theme absorbed
    .groups     = "drop"
  ) %>%
  slice_min(p.adjust, n = 25, with_ties = FALSE) %>%
  mutate(
    neg_log10_p = -log10(p.adjust),
    parentTerm  = str_wrap(parentTerm, width = 45),
    parentTerm  = fct_reorder(parentTerm, neg_log10_p)
  )

# Subtitle with redundancy summary
subtitle_text <- paste0(
  n_terms, " significant GO BP terms (FDR < 0.05)  →  ",
  n_themes, " unique biological themes after semantic similarity clustering",
  "  (", redundancy, "% redundancy)"
)

ggplot(top_themes,
       aes(x = neg_log10_p, y = parentTerm, fill = neg_log10_p)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0("(", n_terms, " terms)")),
            hjust = -0.1, size = 2.8, color = "grey40") +
  scale_fill_gradient(low = "#c6dbef", high = "#08306b",
                      guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title    = "GO BP: Representative biological themes after redundancy reduction",
    subtitle = subtitle_text,
    x        = "-log10(adj. p-value)",
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y   = element_text(size = 9),
    axis.text.x   = element_text(size = 9),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 8, color = "grey30"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin   = margin(t = 10, r = 20, b = 10, l = 10)
  )

ggsave("figures/go_redundancy_barchart.png",
       width = 10, height = 9, dpi = 300)
cat("Saved: figures/go_redundancy_barchart.png\n")

cat("\n=== GO redundancy quantification complete ===\n")
cat("Key result:", n_terms, "significant GO BP terms →",
    n_themes, "unique biological themes (",
    redundancy, "% redundancy)\n")
cat("\nNext step: 05_additional_resources.R (C7, WikiPathways — planned)\n")
