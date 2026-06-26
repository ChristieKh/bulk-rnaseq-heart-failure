# ============================================
# Step 07. Annotate and split DE genes
# ============================================
# Annotates the strict DE genes with gene symbols and full names, splits
# them into "upregulated in HCM" and "higher in AS", and saves the
# annotated tables.

library(org.Hs.eg.db)
library(AnnotationDbi)

# 1. Load DE results
res_df <- readRDS("data/processed/deseq_results.rds")

# Convert the DESeqResults object to a plain data.frame so that
# row/column subsetting below uses standard data.frame behaviour
res_df <- as.data.frame(res_df)

# 2. Create output directory
dir.create("results/tables/07_annotated_de_genes", recursive = TRUE, showWarnings = FALSE)

# 3. Keep significant genes using the strict threshold
res_sig <- res_df[
  !is.na(res_df$padj) &
  res_df$padj < 0.05 &
  abs(res_df$log2FoldChange) > 1,
]

# 4. Sort genes by adjusted p-value
res_sig <- res_sig[order(res_sig$padj), ]

# 5. Extract Ensembl IDs (cleaned)
ensembl_ids <- sub("\\..*$", "", rownames(res_sig))

# 6. Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# 7. Map Ensembl IDs to full gene names
gene_names <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "GENENAME",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# 8. Build annotated DE table
annotated_table <- data.frame(
    gene_symbol = ifelse(is.na(gene_symbols), "", gene_symbols),
  gene_name = ifelse(is.na(gene_names), "", gene_names),
  ensembl_id = ensembl_ids,
  log2FoldChange = res_sig$log2FoldChange,
  padj = res_sig$padj,
  baseMean = res_sig$baseMean,
  stringsAsFactors = FALSE
)

# 9. Split into genes upregulated in HCM and genes higher in AS
up_in_hcm <- annotated_table[annotated_table$log2FoldChange > 1, ]
higher_in_as <- annotated_table[annotated_table$log2FoldChange < -1, ]

# 10. Sort each table by adjusted p-value
up_in_hcm <- up_in_hcm[order(up_in_hcm$padj), ]
higher_in_as <- higher_in_as[order(higher_in_as$padj), ]

# 11. Save full annotated strict DE table
write.csv(
  annotated_table,
  "results/tables/07_annotated_de_genes/all_de_genes_strict_annotated.csv",
  row.names = FALSE
)

# 12. Save genes upregulated in HCM
write.csv(
  up_in_hcm,
  "results/tables/07_annotated_de_genes/upregulated_in_HCM_annotated.csv",
  row.names = FALSE
)

# 13. Save genes higher in AS
write.csv(
  higher_in_as,
  "results/tables/07_annotated_de_genes/higher_in_AS_annotated.csv",
  row.names = FALSE
)

# 14. Save top 30 genes for quick review
top_n <- min(30, nrow(annotated_table))
top_30_annotated <- annotated_table[1:top_n, , drop = FALSE]

write.csv(
  top_30_annotated,
  "results/tables/07_annotated_de_genes/top_30_de_genes_annotated.csv",
  row.names = FALSE
)

# 15. Print summary
cat("========== ANNOTATED DE GENE TABLE SUMMARY ==========\n")
cat("Total strict DE genes:", nrow(annotated_table), "\n")
cat("Upregulated in HCM:", nrow(up_in_hcm), "\n")
cat("Higher in AS:", nrow(higher_in_as), "\n")
cat("Top annotated genes saved:", nrow(top_30_annotated), "\n")