# ============================================
# Step 08. Over-representation analysis (ORA) of DE genes
# ============================================
# GO Biological Process enrichment for genes higher in HCM and higher in AS,
# tested against the background of all expressed (tested) genes.

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# 1. Load full DE results (all tested genes, not only significant)
res <- readRDS("data/processed/deseq_results.rds")
res_df <- as.data.frame(res)

# 2. Drop genes that were not tested (padj = NA)
res_df <- res_df[!is.na(res_df$padj), ]

# 3. Clean Ensembl IDs (remove version suffix like ".15")
res_df$ensembl <- sub("\\..*$", "", rownames(res_df))

# 4. Define the background universe = ALL tested genes (as Entrez IDs)
universe_entrez <- mapIds(
  org.Hs.eg.db,
  keys     = res_df$ensembl,
  column   = "ENTREZID",
  keytype  = "ENSEMBL",
  multiVals = "first"
)
universe_entrez <- unique(na.omit(universe_entrez))

# 5. Define the two gene sets of interest (padj < 0.05, split by direction)
up_in_hcm_ensembl <- res_df$ensembl[res_df$padj < 0.05 & res_df$log2FoldChange > 0]
up_in_as_ensembl  <- res_df$ensembl[res_df$padj < 0.05 & res_df$log2FoldChange < 0]

# 6. Convert each gene set to Entrez IDs
up_in_hcm_entrez <- unique(na.omit(mapIds(
  org.Hs.eg.db, keys = up_in_hcm_ensembl,
  column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first"
)))
up_in_as_entrez <- unique(na.omit(mapIds(
  org.Hs.eg.db, keys = up_in_as_ensembl,
  column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first"
)))

# 7. Sanity check: how many genes survived ID conversion?
cat("Universe (tested) genes:        ", length(universe_entrez), "\n")
cat("Up in HCM: ensembl =", length(up_in_hcm_ensembl), " -> entrez =", length(up_in_hcm_entrez), "\n")
cat("Up in AS:  ensembl =", length(up_in_as_ensembl),  " -> entrez =", length(up_in_as_entrez),  "\n")

# 8. Run GO Biological Process over-representation analysis for each gene set
ego_hcm <- enrichGO(
  gene          = up_in_hcm_entrez,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",            # Biological Process branch of GO
  pAdjustMethod = "BH",            # multiple-testing correction (same idea as padj)
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE             # show gene symbols (not Entrez) in the results
)

ego_as <- enrichGO(
  gene          = up_in_as_entrez,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# 9. How many enriched GO terms were found?
cat("\nEnriched GO:BP terms (HCM-up):", nrow(as.data.frame(ego_hcm)), "\n")
cat("Enriched GO:BP terms (AS-up): ", nrow(as.data.frame(ego_as)),  "\n")

# 10. Peek at the top terms
cat("\n--- Top terms higher in HCM ---\n")
print(head(as.data.frame(ego_hcm)[, c("Description", "GeneRatio", "p.adjust", "Count")], 10))
cat("\n--- Top terms higher in AS ---\n")
print(head(as.data.frame(ego_as)[, c("Description", "GeneRatio", "p.adjust", "Count")], 10))

# 11. Create output directory and save full result tables
dir.create("results/tables/08_enrichment", recursive = TRUE, showWarnings = FALSE)
write.csv(as.data.frame(ego_hcm), "results/tables/08_enrichment/GO_BP_up_in_HCM.csv", row.names = FALSE)
write.csv(as.data.frame(ego_as),  "results/tables/08_enrichment/GO_BP_up_in_AS.csv",  row.names = FALSE)

# 12. Save the enrichResult objects for plotting in the next step
saveRDS(ego_hcm, "data/processed/ego_hcm.rds")
saveRDS(ego_as,  "data/processed/ego_as.rds")
