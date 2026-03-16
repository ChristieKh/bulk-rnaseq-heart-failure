# --------------------------------------------
# 1. Load input files
# --------------------------------------------

counts <- read.delim("data/raw/HCM_vs_stenosis_raw_counts.tsv", row.names = 1, check.names = FALSE)
meta <- read.csv("data/metadata/sample_metadata.csv", stringsAsFactors = FALSE)

# --------------------------------------------
# 2. Quick first look at the objects
# --------------------------------------------

# Show dimensions: number of rows and columns
dim(counts)
dim(meta)
# Show the first few rows
head(counts[, 1:5])
head(meta)

# Show sample names from counts
colnames(counts)
# Show sample names from metadata
meta$sample_id


# --------------------------------------------
# 3. Check the structure of the count matrix
# --------------------------------------------

# Show the internal structure of the counts object
str(counts)

# Check whether all columns are numeric
all_numeric <- all(sapply(counts, is.numeric))
all_numeric

# Check whether all values are integer-like
all_integer_like <- all(counts == round(counts))
all_integer_like

# Check for negative values
has_negative_values <- any(counts < 0)
has_negative_values

# Check for missing values
counts_has_na <- anyNA(counts)
counts_has_na


# --------------------------------------------
# 4. Check gene identifiers
# --------------------------------------------

# Show the first few gene IDs
head(rownames(counts))

# Check for empty or missing gene IDs
has_empty_gene_ids <- any(rownames(counts) == "" | is.na(rownames(counts)))
has_empty_gene_ids

# Check for duplicated gene IDs
duplicated_gene_index <- anyDuplicated(rownames(counts))
duplicated_gene_index

# If duplicates exist, print a few of them
if (duplicated_gene_index != 0) {
  duplicated_genes <- rownames(counts)[duplicated(rownames(counts))]
  head(duplicated_genes)
}


# --------------------------------------------
# 5. Check metadata
# --------------------------------------------

# Show metadata column names
colnames(meta)

# Check that required columns exist
required_columns <- c("sample_id", "condition_short")
missing_required_columns <- setdiff(required_columns, colnames(meta))
missing_required_columns

# Check for missing values in metadata
meta_has_na <- anyNA(meta)
meta_has_na

# Check for empty or missing sample IDs in metadata
has_empty_sample_ids <- any(meta$sample_id == "" | is.na(meta$sample_id))
has_empty_sample_ids

# Check for duplicated sample IDs in metadata
duplicated_sample_index <- anyDuplicated(meta$sample_id)
duplicated_sample_index

# Show the unique condition values before cleaning
unique(meta$condition_short)

# Show counts per group before cleaning
table(meta$condition_short)


# --------------------------------------------
# 6. Clean condition labels
# --------------------------------------------

# Remove accidental spaces
meta$condition_short <- trimws(meta$condition_short)

# Convert to upper case for consistency
meta$condition_short <- toupper(meta$condition_short)

# Show cleaned condition labels
unique(meta$condition_short)
table(meta$condition_short)


# --------------------------------------------
# 7. Check whether sample names match
# --------------------------------------------

# Check whether the set of sample names is the same in both tables
same_sample_set <- setequal(colnames(counts), meta$sample_id)
same_sample_set

# If not the same, show differences
if (!same_sample_set) {
  cat("Samples present in counts but missing in metadata:\n")
  print(setdiff(colnames(counts), meta$sample_id))

  cat("Samples present in metadata but missing in counts:\n")
  print(setdiff(meta$sample_id, colnames(counts)))
}


# --------------------------------------------
# 8. Reorder metadata to match counts columns
# --------------------------------------------

# Reorder metadata rows so that they match the order of columns in counts
meta <- meta[match(colnames(counts), meta$sample_id), ]

# Check that the order now matches exactly
same_sample_order <- all(colnames(counts) == meta$sample_id)
same_sample_order


# --------------------------------------------
# 9. Convert condition to a factor
# --------------------------------------------

# Set AS as reference level and HCM as second level.
# This will make interpretation of contrasts easier later.
meta$condition_short <- factor(meta$condition_short, levels = c("AS", "HCM"))

# Show factor levels
levels(meta$condition_short)

# Show final group sizes
table(meta$condition_short)


# --------------------------------------------
# 10. Final summary of checks
# --------------------------------------------

cat("\n========== FINAL INPUT CHECK SUMMARY ==========\n")
cat("Count matrix dimensions:", dim(counts)[1], "genes x", dim(counts)[2], "samples\n")
cat("Metadata dimensions:", dim(meta)[1], "rows x", dim(meta)[2], "columns\n")
cat("All count columns numeric:", all_numeric, "\n")
cat("All counts integer-like:", all_integer_like, "\n")
cat("Any negative count values:", has_negative_values, "\n")
cat("Any NA in counts:", counts_has_na, "\n")
cat("Any empty/missing gene IDs:", has_empty_gene_ids, "\n")
cat("Duplicated gene IDs index (0 means none):", duplicated_gene_index, "\n")
cat("Any NA in metadata:", meta_has_na, "\n")
cat("Any empty/missing sample IDs:", has_empty_sample_ids, "\n")
cat("Duplicated sample ID index (0 means none):", duplicated_sample_index, "\n")
cat("Same sample set in counts and metadata:", same_sample_set, "\n")
cat("Same sample order after reordering:", same_sample_order, "\n")
cat("Condition levels:", paste(levels(meta$condition_short), collapse = ", "), "\n")
cat("==============================================\n")


# --------------------------------------------
# 11. Stop with an error if critical checks fail
# --------------------------------------------

if (!all_numeric) {
  stop("Not all count columns are numeric.")
}

if (!all_integer_like) {
  stop("Counts are not integer-like. Check whether this is a raw count matrix.")
}

if (has_negative_values) {
  stop("Count matrix contains negative values.")
}

if (counts_has_na) {
  stop("Count matrix contains missing values.")
}

if (has_empty_gene_ids) {
  stop("Some gene IDs are empty or missing.")
}

if (duplicated_gene_index != 0) {
  stop("Duplicated gene IDs found.")
}

if (length(missing_required_columns) > 0) {
  stop("Metadata is missing required columns: ", paste(missing_required_columns, collapse = ", "))
}

if (meta_has_na) {
  stop("Metadata contains missing values.")
}

if (has_empty_sample_ids) {
  stop("Some sample IDs in metadata are empty or missing.")
}

if (duplicated_sample_index != 0) {
  stop("Duplicated sample IDs found in metadata.")
}

if (!same_sample_set) {
  stop("Sample names do not match between count matrix and metadata.")
}

if (!same_sample_order) {
  stop("Sample order still does not match after reordering.")
}


# --------------------------------------------
# 12. Save clean objects for the next step
# --------------------------------------------

# Save cleaned counts and metadata as R objects
saveRDS(counts, "data/raw/counts_checked.rds")
saveRDS(meta, "data/metadata/meta_checked.rds")

cat("\nInput check completed successfully.\n")
cat("Saved cleaned objects:\n")
cat("- data/raw/counts_checked.rds\n")
cat("- data/metadata/meta_checked.rds\n")

