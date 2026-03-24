library(pheatmap)

# 0. Load DE results
vsd <- readRDS("data/processed/vsd.rds")
res_df <- readRDS("data/processed/deseq_results.rds")

# 1. Keep only significant genes
sig_res <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]

# 2. Sort genes by adjusted p-value
sig_res <- sig_res[order(sig_res$padj), ]

# 3. Select top 30 DE genes (or fewer if less than 30 are available)
n_top <- min(30, nrow(sig_res))
top_genes <- rownames(sig_res)[1:n_top]

# 4. Extract VST-normalized expression values for selected genes
mat <- assay(vsd)[top_genes, ]

# 5. Prepare sample annotation
annotation_col <- as.data.frame(colData(vsd)[, "condition_short", drop = FALSE])

# 6. Define annotation colors
annotation_colors <- list(
  condition_short = c(
    "AS" = "#F4B6C2",
    "HCM" = "#76D7EA"
  )
)

# 7. Plot heatmap
pheatmap(
  mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  border_color = NA,
  main = paste0("Top ", n_top, " DE genes: HCM vs AS")
)

# 8. Save heatmap to file
png("results/figures/06_heatmap_top_de_genes.png", width = 1000, height = 900, res = 150)

pheatmap(
  mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  border_color = NA,
  main = paste0("Top ", n_top, " DE genes: HCM vs AS")
)

dev.off()