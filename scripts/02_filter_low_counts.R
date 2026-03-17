# --------------------------------------------
# 1. Load checked input objects
# --------------------------------------------

counts <- readRDS("data/processed/counts_checked.rds")
meta <- readRDS("data/processed/meta_checked.rds")


# --------------------------------------------
# 2. Quick check
# --------------------------------------------

dim(counts)
dim(meta)

head(counts[, 1:min(5, ncol(counts))])
head(meta)

table(meta$condition)


# --------------------------------------------
# 3. Explore gene-level count support
# --------------------------------------------

# For each gene, count in how many samples it has count >= 10
samples_with_count_ge_10 <- rowSums(counts >= 10)

# Quick summary
summary(samples_with_count_ge_10)

# Show how many genes are supported in 0, 1, 2, ... samples
table(samples_with_count_ge_10)


# --------------------------------------------
# 4. Define filtering rule
# --------------------------------------------

# Keep genes with count >= 10 in at least 3 samples
keep_genes <- samples_with_count_ge_10 >= 3

# How many genes pass the filter?
table(keep_genes)

genes_before <- nrow(counts)
genes_after  <- sum(keep_genes)

genes_before
genes_after


# --------------------------------------------
# 5. Apply the filter
# --------------------------------------------

counts_filtered <- counts[keep_genes, ]

dim(counts_filtered)

# Check that metadata stays the same
dim(meta)


# --------------------------------------------
# 6. Summarize filtering result
# --------------------------------------------

genes_removed <- genes_before - genes_after
percent_kept  <- round(100 * genes_after / genes_before, 2)
percent_removed <- round(100 * genes_removed / genes_before, 2)

cat("\n========== LOW-COUNT FILTER SUMMARY ==========\n")
cat("Genes before filtering:", genes_before, "\n")
cat("Genes after filtering:", genes_after, "\n")
cat("Genes removed:", genes_removed, "\n")
cat("Percent kept:", percent_kept, "%\n")
cat("Percent removed:", percent_removed, "%\n")
cat("Filtering rule: keep genes with count >= 10 in at least 3 samples\n")
cat("=============================================\n")


# --------------------------------------------
# 7. Save filtered objects
# --------------------------------------------

saveRDS(counts_filtered, "data/processed/counts_filtered.rds")
saveRDS(meta, "data/processed/meta_filtered.rds")

cat("\nFiltered objects saved:\n")
cat("- data/processed/counts_filtered.rds\n")
cat("- data/processed/meta_filtered.rds\n")