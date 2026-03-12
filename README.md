# Bulk RNA-seq analysis of Alzheimer’s-associated APP mutation in human iPSC-derived neurons

## Project overview

This project is an analytical reanalysis of a public bulk RNA-seq dataset from human iPSC-derived neurons carrying Alzheimer’s-associated mutations.

For the first version of the project, I focus on one clean and interpretable contrast: **WT vs APP mutant neurons**. The aim is to identify gene-expression changes associated with an Alzheimer’s-related neuronal state and to interpret them in a biologically meaningful way.

## Research question

**How does an Alzheimer’s-associated APP mutation affect gene-expression patterns in human iPSC-derived neurons?**

## Why this dataset

I selected **GSE128343** because it is a compact and biologically clear dataset suitable for a first analytical RNA-seq project:

- human samples
- bulk RNA-seq
- iPSC-derived neurons
- Alzheimer’s-related mutations
- relatively small and manageable sample size
- clean experimental setup with a shared cell-line background

This makes it a good dataset for learning how to connect transcriptomic analysis with biological interpretation.

## Biological background

Induced pluripotent stem cells (iPSCs) are cells that can be reprogrammed into a stem-like state and then differentiated into specific cell types. In this dataset, human iPSCs were differentiated into neurons.

This allows researchers to study disease-relevant molecular changes in a controlled human cellular model. Here, neurons carrying Alzheimer’s-associated mutations are compared with wild-type neurons to investigate transcriptomic changes linked to disease biology.

## Dataset

**Accession:** GSE128343  
**Model:** human iPSC-derived neurons  
**Data type:** bulk RNA-seq gene-level counts  
**Current analysis scope:** WT vs APP

### Samples used in this project

| Group | Sample IDs | n |
|------|------|---|
| WT | 1222, 1223, 1224 | 3 |
| APP | 1225, 1226, 1227 | 3 |

### Samples not yet included

| Group | Sample IDs | n |
|------|------|---|
| PSEN1 | 1228, 1229, 1230 | 3 |
| APP/PSEN1 | 1231, 1232, 1233 | 3 |

These groups are kept for future project extensions.

## Input data

The project uses a public processed count matrix provided on GEO together with sample metadata.

- `data/raw/GSE128343_RNAseq_raw.txt` — gene-by-sample count matrix
- `data/raw/GSE128343_series_matrix.txt` — dataset metadata and sample annotations
- `data/metadata/sample_table_wt_vs_app.csv` — analysis-specific sample table for the first contrast

## Analytical workflow

1. Load count matrix and sample metadata
2. Subset samples for WT vs APP comparison
3. Perform count-level QC and exploratory analysis
4. Run differential expression analysis
5. Visualize transcriptomic differences
6. Perform functional enrichment analysis
7. Interpret results in the context of neuronal dysfunction and Alzheimer’s-related biology

## Main outputs

- sample metadata table
- filtered count matrix
- PCA plot
- sample distance heatmap
- volcano plot
- top differentially expressed genes
- functional enrichment results
- biological interpretation notes

## Repository structure

```
text
data/
scripts/
results/
docs/
```