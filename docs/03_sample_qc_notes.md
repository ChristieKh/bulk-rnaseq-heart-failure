# Step 3 — DESeq2 setup and sample-level QC

## Overview

In this step, the filtered count matrix was prepared for downstream differential expression analysis using DESeq2.  
Sample-level quality control (QC) was performed to assess global data structure, detect potential outliers, and evaluate whether samples are comparable after normalization.

---

## DESeq2 dataset construction

A `DESeqDataSet` object was created using:

- filtered count matrix
- curated sample metadata

Design formula:
~ condition_short

Conditions:

- AS (reference)
- HCM

This setup defines the statistical comparison for downstream differential expression analysis.

---

## Size factor estimation

DESeq2 size factors were estimated to normalize for differences in sequencing depth across samples.

Observed range:

- approximately **0.80 – 1.33**

### Interpretation

- moderate variation in library scale
- no extreme outliers

### Conclusion

Library size differences are present but well within expected range and can be handled by DESeq2 normalization.

---

## Library size assessment

A barplot of total raw counts per sample was generated.

### Observation

- all samples have comparable sequencing depth (~30–32 million reads)
- no extreme low-depth or high-depth samples
- no obvious imbalance between AS and HCM groups

### Conclusion

Sequencing depth is consistent across samples, and no sample appears problematic at the raw count level.

---

## VST transformation

Variance stabilizing transformation (VST) was applied to normalized counts.

### Purpose

- reduce dependence of variance on mean expression
- enable meaningful comparison between samples
- prepare data for PCA and distance-based analyses

---

## Distribution of expression values

A boxplot of VST-transformed expression values was generated.

### Observation

- highly similar medians across samples
- comparable interquartile ranges
- similar overall distribution shapes
- no obvious outliers

### Conclusion

Normalization and transformation were successful, and global expression distributions are consistent across samples.

---

## Principal Component Analysis (PCA)

PCA was performed on VST-transformed data.

### Variance explained

- PC1: **22%**
- PC2: **18%**
- total (PC1 + PC2): **40%**

### Observation

- no clear separation between AS and HCM samples
- substantial within-group variability
- especially pronounced dispersion among HCM samples

### Interpretation

- condition is not the dominant source of global variation
- other factors may contribute, such as:
  - biological heterogeneity between individuals
  - age or sex differences
  - disease severity
  - technical variation

### Important note

Lack of clear separation in PCA does **not** imply absence of meaningful biological signal.  
Differential expression may still be present at the gene or pathway level.

---

## Sample-to-sample distance heatmap

A heatmap based on pairwise sample distances (Euclidean distance on VST data) was generated.

### Observation

- no strict clustering of all HCM samples separately from AS
- presence of local clusters of similar samples
- overall heterogeneous structure

### Interpretation

- confirms PCA findings
- indicates substantial variability within both groups

---

## Final QC summary

- Library sizes are consistent across samples
- VST-transformed distributions are highly similar
- No strong sample-level outliers detected
- PCA and distance analysis show substantial within-group heterogeneity
- Condition (AS vs HCM) is not the sole driver of global transcriptomic variation

### Overall conclusion

The dataset appears **technically sound and suitable for downstream differential expression analysis**, despite notable biological heterogeneity between samples.