# ============================================
# Step 06. Heatmap of top DE genes
# ============================================
# Selects the top strict DE genes (padj < 0.05 & |log2FC| > 1), extracts
# their VST expression, and draws row-scaled heatmaps across all samples,
# labelled first by Ensembl ID and then by gene symbol.

library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)

# 1. Load DE results
vsd <- readRDS("data/processed/vsd.rds")
res_df <- readRDS("data/processed/deseq_results.rds")

# Convert the DESeqResults object to a plain data.frame so that the
# row subsetting below uses standard data.frame behaviour
res_df <- as.data.frame(res_df)

# 2. Create output directory
dir.create("results/figures/06_heatmap_top_genes", recursive = TRUE, showWarnings = FALSE)

# 3. Keep only significant genes

# non-strict filter
# sig_res <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]

sig_res <- res_df[
  !is.na(res_df$padj) &
  res_df$padj < 0.05 &
  abs(res_df$log2FoldChange) > 1,
]

# 4. Sort genes by adjusted p-value
sig_res <- sig_res[order(sig_res$padj), ]

# 5. Select top 30 DE genes (or fewer if less than 30 are available)
n_top <- min(30, nrow(sig_res))
top_genes <- rownames(sig_res)[1:n_top]

# 6. Extract VST-normalized expression values for selected genes
mat <- assay(vsd)[top_genes, ]

# 7. Prepare sample annotation
annotation_col <- as.data.frame(colData(vsd)[, "condition_short", drop = FALSE])

# 8. Define annotation colors
annotation_colors <- list(
  condition_short = c(
    "AS" = "#F4B6C2",
    "HCM" = "#76D7EA"
  )
)

# 9. Heatmap with Ensembl IDs
png("results/figures/06_heatmap_top_genes/heatmap_top_de_genes_strict.png", width = 1000, height = 900, res = 150)

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


# 10. Map Ensembl IDs → gene symbols

# Map Ensembl IDs with version numbers removed
top_genes_clean <- sub("\\..*$", "", top_genes)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = top_genes_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# 11. Replace rownames with gene symbols (fallback to Ensembl if NA)
rownames(mat) <- ifelse(is.na(gene_symbols), top_genes, gene_symbols)


# 12. Heatmap with gene symbols
png("results/figures/06_heatmap_top_genes/heatmap_top_de_genes_symbols_strict.png", width = 1000, height = 900, res = 150)

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
  main = paste0("Top ", n_top, " DE genes (symbols): HCM vs AS")
)

dev.off()
