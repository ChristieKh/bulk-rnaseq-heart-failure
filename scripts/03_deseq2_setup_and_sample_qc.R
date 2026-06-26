# ============================================
# Step 03. DESeq2 setup and sample-level QC
# ============================================
# Builds the DESeq2 dataset (design = ~ condition_short), estimates size
# factors, applies the variance-stabilizing transformation, and produces
# sample-level QC plots: library sizes, VST boxplot, PCA, and a
# sample-to-sample distance heatmap.

library(DESeq2)

# 2. Load filtered input objects
counts <- readRDS("data/processed/counts_filtered.rds")
meta   <- readRDS("data/processed/meta_filtered.rds")

# 3. Quick check
dim(counts)
dim(meta)
head(counts[, 1:min(5, ncol(counts))])
head(meta)
table(meta$condition_short)

# 3.1 Check and delete all additional rownames such as "N_ambiguous" "N_multimapping" "N_noFeature" "N_unmapped"
tail(rownames(counts), 20)
counts <- counts[!grepl("^N_", rownames(counts)), ]
dim(counts)
tail(rownames(counts), 10)

# 3.2 Define levels 
meta$condition_short <- factor(meta$condition_short, levels = c("AS", "HCM"))
levels(meta$condition_short)

# 4. Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ condition_short
)

dds

# 5. Check group sizes inside DESeq2 object
table(colData(dds)$condition_short)

# 6. Estimate size factors (normalization factors)
dds <- estimateSizeFactors(dds)

sizeFactors(dds)

# 7. Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

dim(normalized_counts)
head(normalized_counts[, 1:min(5, ncol(normalized_counts))])


# 8. Apply variance stabilizing transformation
vsd <- vst(dds, blind = TRUE)

vsd

# 9. Save DESeq2 objects for downstream steps
saveRDS(dds, "data/processed/dds.rds")
saveRDS(vsd, "data/processed/vsd.rds")
saveRDS(normalized_counts, "data/processed/normalized_counts.rds")

cat("\nSaved processed objects:\n")
cat("- data/processed/dds.rds\n")
cat("- data/processed/vsd.rds\n")
cat("- data/processed/normalized_counts.rds\n")

#11. Library sizes (raw counts)
library_sizes <- colSums(counts)
library_sizes

png("results/figures/03_sample_qc/library_sizes.png", width = 1200, height = 800)
barplot(
  library_sizes / 1e6,
  las = 2,
  main = "Library sizes across samples",
  ylab = "Total raw counts (millions)"
)
dev.off()

# 12. Boxplot of VST-transformed values
png("results/figures/03_sample_qc/vst_boxplot.png", width = 1200, height = 800)
boxplot(
  assay(vsd),
  las = 2,
  main = "Distribution of VST-transformed expression values",
  ylab = "VST expression"
)
dev.off()

# 13. PCA plot
png("results/figures/03_sample_qc/pca_plot.png", width = 1200, height = 800)
plotPCA(vsd, intgroup = "condition_short")
dev.off()

# 14. Sample-to-sample distance heatmap
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)

png("results/figures/03_sample_qc/sample_distance_heatmap.png", width = 1200, height = 1000)
heatmap(
  sample_dist_matrix,
  symm = TRUE,
  main = "Sample-to-sample distance heatmap"
)
dev.off()