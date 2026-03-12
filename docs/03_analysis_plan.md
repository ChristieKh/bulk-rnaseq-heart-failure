# Analysis plan

## Main question

How does an Alzheimer’s-associated APP mutation affect gene-expression patterns in human iPSC-derived neurons?

## Contrast

WT vs APP

## Planned steps

### 1. Data loading and preparation
- import count matrix
- import sample table
- subset samples for the selected contrast
- check dimensions and column order

### 2. Count-level QC
- inspect library sizes
- inspect count distributions
- filter low-count genes
- perform variance-stabilizing transformation if needed

### 3. Exploratory analysis
- PCA
- sample distance heatmap
- check whether biological groups separate

### 4. Differential expression
- run DESeq2
- extract log2 fold changes and adjusted p-values
- define significance thresholds

### 5. Visualization
- PCA plot
- volcano plot
- heatmap of top DE genes

### 6. Functional interpretation
- GO enrichment
- pathway enrichment
- interpretation of major functional themes

### 7. Reporting
- summarize major transcriptomic patterns
- discuss biological meaning
- state limitations clearly

## Scope boundaries

This project intentionally excludes:
- raw FASTQ preprocessing
- full replication of the published paper
- multi-contrast analysis in version 1
- advanced batch correction or complex modeling