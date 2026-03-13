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

- `data/raw/GSE206978_HCM_vs_stenosis_raw_counts.tsv.gz` — gene-level raw count matrix
- `data/raw/GSE206978_Sample_description.tsv` — sample-level metadata from GEO
- `data/metadata/sample_table_hcm_vs_as.csv` — cleaned analysis-ready sample table

## Analytical workflow

1. Import raw counts and sample metadata  
2. Build an analysis-ready sample table  
3. Perform sample-level quality checks and exploratory analysis  
4. Apply count transformation for visualization  
5. Run differential expression analysis with **DESeq2**  
6. Generate PCA, volcano plot, and heatmap of top differentially expressed genes  
7. Perform functional enrichment analysis  
8. Interpret results in the context of cardiac hypertrophic remodeling

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