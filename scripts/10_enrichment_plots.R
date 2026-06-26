# ============================================
# Step 10. Enrichment visualisation
# ============================================
# Dot plot of the GSEA GO:BP results using enrichplot::dotplot.
# split = ".sign" separates terms enriched in HCM (NES > 0) from terms
# enriched in AS (NES < 0); the panels are relabelled for the reader.
#   - x position = gene ratio (leading-edge genes / set size)
#   - dot size   = number of leading-edge genes
#   - dot colour = NES (enrichment strength and direction)

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# 1. Load the saved GSEA result object
gsea_bp <- readRDS("data/processed/gsea_bp.rds")

# 2. Output directory for figures
dir.create("results/figures/10_enrichment", recursive = TRUE, showWarnings = FALSE)

# 3. Build the dot plot: top N terms per side, coloured by NES,
#    with the activated/suppressed panels relabelled to HCM/AS
n_show <- 10
p <- dotplot(gsea_bp, showCategory = n_show, split = ".sign", color = "NES") +
  facet_grid(
    . ~ .sign,
    labeller = as_labeller(c(activated  = "Higher in HCM",
                             suppressed = "Higher in AS"))
  ) +
  labs(title = "GSEA (GO Biological Process): HCM vs AS")

# 4. Save the figure
ggsave(
  "results/figures/10_enrichment/gsea_dotplot.png",
  plot = p, width = 11, height = 8, dpi = 150
)
cat("Saved: results/figures/10_enrichment/gsea_dotplot.png\n")
