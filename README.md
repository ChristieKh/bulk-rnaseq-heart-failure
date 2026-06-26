# Bulk RNA-seq analysis of primary versus pressure-overload cardiac hypertrophy in human myocardium

A reproducible bulk RNA-seq reanalysis of human myocardial tissue comparing hypertrophic cardiomyopathy and aortic-stenosis-induced hypertrophy.

## Project overview

This repository contains a focused transcriptomic reanalysis of human cardiac hypertrophy using public bulk RNA-seq data.

The project examines gene-expression differences between two pathological remodeling states in the human heart:

- **hypertrophic cardiomyopathy (HCM)**, representing primary myocardial hypertrophy
- **aortic-stenosis-induced hypertrophy (AS)**, representing secondary hypertrophy caused by chronic pressure overload

The analysis is designed as a compact, interpretable RNA-seq project centered on raw-count-based differential expression, sample-level exploration, and pathway-oriented biological interpretation.

## Research question

**How does hypertrophic cardiomyopathy differ transcriptionally from pressure-overload cardiac hypertrophy in human myocardium?**

## Source study

This project is based on the public dataset [**GSE206978**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206978), generated in the study:

[**Novel Genes Involved in Hypertrophic Cardiomyopathy: Data of Transcriptome and Methylome Profiling**](https://pmc.ncbi.nlm.nih.gov/articles/PMC9739701/)

The original study compared myocardial samples from patients with **HCM** and **AS** and used both **bulk RNA-seq** and **genome-wide DNA methylation profiling**.

This repository focuses only on the **RNA-seq component** and reanalyzes the public raw count matrix as an independent transcriptomic workflow.

## Biological context

Cardiac hypertrophy is not a single biological entity. Similar tissue-level enlargement can arise through different mechanisms.

In **HCM**, hypertrophy reflects a primary disorder of the myocardium, often linked to intrinsic abnormalities of cardiac muscle structure and function. In **AS**, hypertrophy develops as an adaptive response to long-term pressure overload, because the left ventricle must pump against an obstructed aortic valve.

This makes the HCM-versus-AS comparison biologically informative: both groups show hypertrophied myocardium, but the underlying drivers of remodeling differ. The project therefore focuses on transcriptomic features that distinguish **primary myocardial disease** from **pressure-overload remodeling**.

## Why this dataset

I selected **GSE206978** because it is well suited to a first focused analytical RNA-seq project in cardiac disease.

It provides:
- human myocardial tissue
- bulk RNA-seq data
- a gene-level raw count matrix suitable for **DESeq2**
- a clearly defined comparison with direct biological meaning
- a manageable sample size for exploratory and differential-expression analysis

This makes it a practical dataset for building a reproducible workflow while still engaging with a meaningful disease-related question.

## Dataset

**Accession:** GSE206978  
**Tissue:** human myocardium  
**Data type:** bulk RNA-seq raw counts  
**Current comparison:** HCM vs AS

### Samples included in the current analysis

| Group | n |
|------|---:|
| HCM | 8 |
| AS | 5 |

## Input data

The analysis uses public GEO files together with a cleaned project-specific sample sheet.

- `data/raw/HCM_vs_stenosis_raw_counts.tsv` — gene-level raw count matrix (60,723 genes × 13 samples, Ensembl IDs)
- `data/raw/Sample_description.tsv` — sample-level metadata from GEO
- `data/raw/sample_metadata.csv` — cleaned, analysis-ready sample table (sample ID, condition, sex, age)

Validated and intermediate objects produced during the pipeline (checked counts, filtered matrix, DESeq2 objects, results tables) are written to `data/processed/`.

## Analytical workflow

The pipeline is organised as a sequence of numbered R scripts in `scripts/`,
each consuming the output of the previous step:

| Step | Script | What it does |
|---|---|---|
| 01 | `01_load_and_validate_input.R` | Load raw counts + metadata, run integrity checks, align metadata to count columns, set `condition` factor (AS as reference) |
| 02 | `02_filter_low_counts.R` | Group-aware filtering of low-count genes (≥10 counts in ≥3 samples of either group) |
| 03 | `03_deseq2_setup_and_sample_qc.R` | Build DESeq2 dataset, size-factor normalization, VST, sample-level QC (library sizes, PCA, sample-distance heatmap) |
| 04 | `04_differential_expression.R` | Differential expression with **DESeq2** (HCM vs AS), result tables, MA plot |
| 05 | `05_de_visualization.R` | Volcano plot of DE results |
| 06 | `06_heatmap_top_de_genes.R` | Heatmap of the top differentially expressed genes |
| 07 | `07_annotate_top_de_genes.R` | Annotate DE genes (Ensembl → symbol + name), split into up-in-HCM / higher-in-AS |
| 08 | `08_enrichment_ora.R` | Over-representation analysis (GO:BP) of the DE gene sets |
| 09 | `09_gsea.R` | Gene Set Enrichment Analysis (GO:BP) over all genes ranked by log2FC |
| 10 | `10_enrichment_plots.R` | Dot-plot visualisation of the GSEA results |

## Environment and how to run

The analysis runs in R (4.4) with Bioconductor packages, managed through a
conda environment defined in `environment.yml` (key packages: **DESeq2**,
**clusterProfiler**, **enrichplot**, **org.Hs.eg.db**, **pheatmap**, **ggplot2**).

```bash
# create and activate the environment
conda env create -f environment.yml
conda activate hcm-vs-as-rnaseq

# run the pipeline in order
Rscript scripts/01_load_and_validate_input.R
Rscript scripts/02_filter_low_counts.R
# ... through to:
Rscript scripts/10_enrichment_plots.R
```

Each script reads from `data/processed/` (and `data/raw/` for step 01) and
writes its outputs to `results/` and `data/processed/`.

## Main outputs

- cleaned metadata table
- analysis-ready sample sheet
- PCA plot
- sample distance heatmap
- volcano plot
- heatmap of top differentially expressed genes
- differential expression results table
- functional enrichment results
- interpretation notes for key biological patterns

## Repository structure

```text
data/
scripts/
results/
docs/
```

## Project scope

This repository focuses on analytical interpretation of a public bulk RNA-seq dataset using gene-level raw counts.

The current analysis is designed to emphasize:

- study design and contrast definition
- metadata handling
- exploratory transcriptomic analysis
- differential expression analysis with DESeq2
- pathway-level interpretation
- biologically informed discussion of cardiac remodeling

The associated study also included genome-wide DNA methylation analysis, but this repository is currently limited to the RNA-seq component.

## Limitations

- small sample size
- comparison between two pathological states rather than disease versus healthy tissue
- biological conclusions should be interpreted as transcriptomic associations rather than causal mechanisms
- analysis starts from public processed gene-level counts rather than raw FASTQ files
- the differential-expression design models condition only (`~ condition_short`) and does not adjust for sex or age, which are available in the metadata; given the small sample size these potential confounders are not modelled, so some differences may partly reflect demographic imbalance between groups
- exploratory analysis (PCA, sample-distance heatmap) shows no clear global separation between HCM and AS, indicating the transcriptomic differences are subtle and confined to a limited set of genes rather than a genome-wide shift