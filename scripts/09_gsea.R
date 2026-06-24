# ============================================
# Step 09. Gene Set Enrichment Analysis (GSEA) on GO Biological Process
# ============================================
# Unlike ORA (script 08), GSEA uses ALL tested genes ranked by
# log2FoldChange, with no significance cutoff. It detects GO terms whose
# genes are collectively shifted toward one end of the ranking.

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


# 1. Load full DE results and convert to a plain data.frame
res <- readRDS("data/processed/deseq_results.rds")
res_df <- as.data.frame(res)

# 2. Drop genes that were not tested (padj = NA) and strip Ensembl versions
res_df <- res_df[!is.na(res_df$padj), ]
res_df$ensembl <- sub("\\..*$", "", rownames(res_df))

# 3. Map Ensembl -> Entrez (GSEA on GO needs Entrez IDs)
res_df$entrez <- mapIds(
  org.Hs.eg.db, keys = res_df$ensembl,
  column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first"
)

# 4. Keep only genes with a usable, unique Entrez ID
#    (GSEA cannot have missing or duplicated names in the ranked list)
res_df <- res_df[!is.na(res_df$entrez) & !duplicated(res_df$entrez), ]

# 5. Build the ranked gene list: a named numeric vector of log2FoldChange,
#    named by Entrez ID, sorted from most up-in-HCM to most up-in-AS
gene_list <- res_df$log2FoldChange
names(gene_list) <- res_df$entrez
gene_list <- sort(gene_list, decreasing = TRUE)

# 6. Inspect the top (HCM) and bottom (AS) of the ranked list
cat("Length of ranked list:", length(gene_list), "\n")
cat("\nTop of list (most up in HCM):\n");  print(head(gene_list))
cat("\nBottom of list (most up in AS):\n"); print(tail(gene_list))

# 7. Run GSEA on GO Biological Process
#    (permutation-based, so fix the seed for reproducible results)
set.seed(42)
gsea_bp <- gseGO(
  geneList      = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

# 8. Total number of enriched GO:BP terms
gsea_df <- as.data.frame(gsea_bp)
cat("\nGSEA GO:BP terms found:", nrow(gsea_df), "\n")

# 9. Split results by direction of enrichment (sign of NES)
#    NES > 0 -> shifted toward the top of the ranking    -> enriched in HCM
#    NES < 0 -> shifted toward the bottom of the ranking  -> enriched in AS
gsea_hcm <- gsea_df[gsea_df$NES > 0, ]
gsea_as  <- gsea_df[gsea_df$NES < 0, ]

# 10. Sort each side by adjusted p-value and show the top terms
gsea_hcm <- gsea_hcm[order(gsea_hcm$p.adjust), ]
gsea_as  <- gsea_as[order(gsea_as$p.adjust), ]

cat("\nHCM side (NES > 0):", nrow(gsea_hcm), "terms\n")
print(head(gsea_hcm[, c("Description", "NES", "p.adjust", "setSize")], 10))

cat("\nAS side (NES < 0):", nrow(gsea_as), "terms\n")
print(head(gsea_as[, c("Description", "NES", "p.adjust", "setSize")], 10))

# 11. Save the full GSEA table and the result object for plotting
dir.create("results/tables/09_gsea", recursive = TRUE, showWarnings = FALSE)
write.csv(gsea_df, "results/tables/09_gsea/GSEA_GO_BP.csv", row.names = FALSE)
saveRDS(gsea_bp, "data/processed/gsea_bp.rds")
