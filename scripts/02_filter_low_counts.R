# ============================================
# Step 02. Filter low-count genes
# ============================================
# Removes genes with too few reads to be informative, using group-aware
# filtering: a gene is kept if it has count >= 10 in at least 3 samples
# of EITHER condition (so genes specific to one group are not lost).

# 1. Load checked input
counts <- readRDS("data/processed/counts_checked.rds")
meta <- readRDS("data/processed/meta_checked.rds")


# 2. Quick check
dim(counts)
dim(meta)

head(counts[, 1:min(5, ncol(counts))])
head(meta)

table(meta$condition_short)


# 3. Explore gene-level count support by condition

# Logical vectors for sample groups
is_hcm <- meta$condition_short == "HCM"
is_as  <- meta$condition_short == "AS"

# For each gene, count in how many HCM samples it has count >= 10
hcm_support <- rowSums(counts[, is_hcm] >= 10)

# For each gene, count in how many AS samples it has count >= 10
as_support <- rowSums(counts[, is_as] >= 10)

# Quick summaries
summary(hcm_support)
summary(as_support)

# Show how many genes are supported in 0, 1, 2, ... samples within each group
table(hcm_support)
table(as_support)


# 4. Define group-aware filtering rule

# Keep genes with count >= 10 in at least 3 HCM samples
# OR in at least 3 AS samples
keep_genes <- (hcm_support >= 3) |
              (as_support >= 3)

# How many genes pass the filter?
table(keep_genes)

# Optional: create filtered count matrix
counts_filtered <- counts[keep_genes, ]

# Check dimensions after filtering
dim(counts_filtered)


genes_before <- nrow(counts)
genes_after  <- sum(keep_genes)

genes_before
genes_after

# 5. Apply the filter

counts_filtered <- counts[keep_genes, ]

dim(counts_filtered)

# Check that metadata stays the same
dim(meta)


# 6. Summarize filtering result

genes_removed <- genes_before - genes_after
percent_kept  <- round(100 * genes_after / genes_before, 2)
percent_removed <- round(100 * genes_removed / genes_before, 2)

cat("\n========== LOW-COUNT FILTER SUMMARY ==========\n")
cat("Genes before filtering:", genes_before, "\n")
cat("Genes after filtering:", genes_after, "\n")
cat("Genes removed:", genes_removed, "\n")
cat("Percent kept:", percent_kept, "%\n")
cat("Percent removed:", percent_removed, "%\n")
cat("Filtering rule: keep genes with count >= 10 in at least 3 samples within either HCM or AS group (group-aware filtering)\n")
cat("=============================================\n")


# 7. Save filtered objects

saveRDS(counts_filtered, "data/processed/counts_filtered.rds")
saveRDS(meta, "data/processed/meta_filtered.rds")

cat("\nFiltered objects saved:\n")
cat("- data/processed/counts_filtered.rds\n")
cat("- data/processed/meta_filtered.rds\n")