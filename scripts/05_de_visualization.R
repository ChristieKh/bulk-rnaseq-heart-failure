# ============================================
# Step 05. Differential expression visualization
# ============================================

library(ggplot2)


# 1. Load DE results
res <- readRDS("data/processed/deseq_results.rds")
res_ordered <- readRDS("data/processed/deseq_results_ordered.rds")


# 2. Create output directory
dir.create("results/figures/05_de_visualization", recursive = TRUE, showWarnings = FALSE)


# 3. Prepare dataframe for plotting
res_df <- as.data.frame(res_ordered)
res_df <- res_df[!is.na(res_df$padj), ]


# 4. Define significance groups
res_df$significance <- "NS"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up in HCM"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Up in AS"

table(res_df$significance)


# 5. Add -log10 adjusted p-value
res_df$neg_log10_padj <- -log10(res_df$padj)


# 6. Volcano plot
png("results/figures/05_de_visualization/volcano_plot.png", width = 1200, height = 800)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(
    data = subset(res_df, significance == "NS"),
    color = "grey80", alpha = 0.5, size = 1.2
  ) +
  geom_point(
    data = subset(res_df, significance != "NS"),
    aes(fill = significance),
    shape = 21, color = "black", stroke = 0.35, size = 2.4, alpha = 0.95
  ) +
  scale_fill_manual(values = c(
    "Up in AS" = "#3775C1",
    "Up in HCM" = "#E64B35"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey35") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey35") +
  coord_cartesian(xlim = c(-3, 3)) +
  labs(
    title = "Volcano plot: HCM vs AS",
    x = "log2 fold change",
    y = "-log10 adjusted p-value"
  ) +
  theme_minimal(base_size = 12)+
  theme(panel.grid.minor = element_blank())

dev.off()


# 7. Save plotting dataframe
write.csv(res_df,
          "results/tables/04_differential_expression/deseq_results_for_volcano.csv",
          row.names = TRUE)