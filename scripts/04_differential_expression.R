# ============================================
# Step 04. Differential expression analysis
# ============================================

library(DESeq2)


# 1. Load DESeq2 object
dds <- readRDS("data/processed/dds.rds")


# 2. Create output directories
dir.create("results/tables/04_differential_expression", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/04_differential_expression", recursive = TRUE, showWarnings = FALSE)


# 3. Run DESeq2 model
dds <- DESeq(dds)

dds


# 4. Extract results
res <- results(dds)

res
summary(res)


# 5. Order results and define significant genes
res_ordered <- res[order(res$padj), ]

res_clean <- res_ordered[!is.na(res_ordered$padj), ]

sig_res <- subset(res_clean, padj < 0.05)
sig_res_fc <- subset(res_clean, padj < 0.05 & abs(log2FoldChange) > 1)

head(res_ordered)
nrow(sig_res)
nrow(sig_res_fc)


# 6. Define up/downregulated gene sets
up_genes <- subset(res_clean, padj < 0.05 & log2FoldChange > 0)
down_genes <- subset(res_clean, padj < 0.05 & log2FoldChange < 0)

nrow(up_genes)
nrow(down_genes)


# 7. Save result tables
write.csv(as.data.frame(res_ordered),
          "results/tables/04_differential_expression/deseq2_results_all.csv")
write.csv(as.data.frame(sig_res),
          "results/tables/04_differential_expression/deseq2_results_significant.csv")
write.csv(as.data.frame(sig_res_fc),
          "results/tables/04_differential_expression/deseq2_results_significant_fc.csv")
write.csv(as.data.frame(up_genes),
          "results/tables/04_differential_expression/upregulated_genes.csv")
write.csv(as.data.frame(down_genes),
          "results/tables/04_differential_expression/downregulated_genes.csv")


# 8. Differential expression summary
total_tested <- nrow(res)
significant_padj_005 <- nrow(sig_res)
significant_padj_005_fc1 <- nrow(sig_res_fc)
upregulated_005 <- nrow(up_genes)
downregulated_005 <- nrow(down_genes)

percent_significant <- round(100 * significant_padj_005 / total_tested, 2)
percent_significant_fc1 <- round(100 * significant_padj_005_fc1 / total_tested, 2)
percent_up <- round(100 * upregulated_005 / total_tested, 2)
percent_down <- round(100 * downregulated_005 / total_tested, 2)

cat("\n========== DIFFERENTIAL EXPRESSION SUMMARY ==========\n")
cat("Genes tested:", total_tested, "\n")
cat("Significant genes (padj < 0.05):", significant_padj_005, "\n")
cat("Significant genes (padj < 0.05 & |log2FC| > 1):", significant_padj_005_fc1, "\n")
cat("Upregulated in HCM (padj < 0.05):", upregulated_005, "\n")
cat("Downregulated in HCM / higher in AS (padj < 0.05):", downregulated_005, "\n")
cat("Percent significant:", percent_significant, "%\n")
cat("Percent significant with |log2FC| > 1:", percent_significant_fc1, "%\n")
cat("Percent upregulated in HCM:", percent_up, "%\n")
cat("Percent downregulated in HCM:", percent_down, "%\n")
cat("Comparison: HCM vs AS\n")
cat("Reference group: AS\n")
cat("====================================================\n")


# 9. Save DESeq2 object with model results
saveRDS(dds, "data/processed/dds_deseq.rds")
saveRDS(res, "data/processed/deseq_results.rds")
saveRDS(res_ordered, "data/processed/deseq_results_ordered.rds")


# 10. MA plot
png("results/figures/04_differential_expression/ma_plot.png", width = 1200, height = 800)
plotMA(res, main = "MA plot: HCM vs AS", ylim = c(-4, 4),
        colSig = "red")
dev.off()
