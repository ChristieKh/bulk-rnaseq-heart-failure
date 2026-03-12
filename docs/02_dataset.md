# Dataset notes

## GEO accession

GSE128343

## Study design

Human iPS cells with different Alzheimer’s-associated mutations were differentiated into neurons and analyzed by RNA-seq.

The full dataset contains 12 samples across 4 groups:

- WT
- APP
- PSEN1
- APP/PSEN1

Each group contains 3 samples.

## Samples used in version 1

For the first version of the project, only the following samples are used:

- WT: 1222, 1223, 1224
- APP: 1225, 1226, 1227

This restriction helps keep the analysis simple, interpretable, and focused on one biologically meaningful comparison.

## Input files

### 1. GSE128343_RNAseq_raw.txt
This file is used as the main count matrix for analysis.

Structure:
- rows: genes
- columns: samples
- values: integer gene-level counts

### 2. GSE128343_series_matrix.txt
This file is used to understand sample annotations and study metadata.

Useful information:
- sample identities
- group labels
- model description
- experimental context

## Why start from processed counts?

The goal of this first project is to focus on transcriptomic reasoning and interpretation rather than on rebuilding the entire RNA-seq preprocessing workflow.

Starting from gene-level counts allows the analysis to begin at the most informative analytical stage:
- sample QC
- exploratory analysis
- differential expression
- pathway interpretation