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

## Repository structure

```text
data/
  raw/                 # input count matrix + sample metadata
  processed/           # validated counts, DESeq2 objects, results (.rds)
scripts/               # numbered pipeline scripts 01–10
results/
  figures/             # QC, volcano, heatmap, GSEA dot plot
  tables/              # DE results and enrichment tables
environment.yml        # conda environment definition
```

## Results

### 1. Sample-level QC and exploratory analysis

Library sizes were comparable across all 13 samples, and no sample behaved as
a technical outlier. However, principal component analysis did **not** separate
HCM and AS: the two groups overlap across the plot, and the leading components
explain only a modest fraction of the variance (PC1 ≈ 22%, PC2 ≈ 19%). The
sample-to-sample distance heatmap tells the same story — samples do not cluster
cleanly by condition.

![PCA of HCM vs AS](results/figures/03_sample_qc/pca_plot.png)

This is an informative negative result rather than a failure: both groups are
hypertrophied, diseased myocardium, so their global transcriptomes are similar.
It sets the expectation that any HCM-vs-AS differences will be **subtle and
confined to a limited set of genes**, not a genome-wide shift.

### 2. Differential expression (DESeq2)

Of **16,693** tested genes, **109** were differentially expressed at
`padj < 0.05`, and **35** passed the stricter threshold of
`padj < 0.05 & |log2FoldChange| > 1`. With AS as the reference level, the
results are skewed toward genes that are **higher in AS / lower in HCM**
(67 vs 42 at `padj < 0.05`), consistent with the original study.

Representative genes:

| Direction | Genes |
|---|---|
| Higher in HCM | `CTXND1`, `EIF4EBP3`, `PCDHGC4`, `ATRNL1` |
| Higher in AS (lower in HCM) | `IGF2`, `SPOCK1`, `ITGA11`, `C4B`, `KCNT1` |

The volcano plot shows the overall picture (blue = higher in AS, red = higher
in HCM), and the heatmap of the top genes confirms that these differences are
consistent across samples rather than driven by single outliers.

![Volcano plot](results/figures/05_de_visualization/volcano_plot.png)

![Heatmap of top DE genes](results/figures/06_heatmap_top_genes/heatmap_top_de_genes_symbols_strict.png)

### 3. Functional enrichment

**Over-representation analysis (ORA, GO:BP).** Run on the significant gene sets,
ORA found enrichment only on the AS side: genes higher in AS were
over-represented for **cell-substrate adhesion** (e.g. `ITGA11`, `SPOCK1`,
`CCDC80`, `LAMA5`, `FLNA`, `NOTCH1`). The shorter HCM gene list yielded no
significant terms — expected for ORA on a small, scattered set.

**Gene set enrichment analysis (GSEA, GO:BP).** Using all genes ranked by
log2FoldChange (no cutoff), GSEA was far more sensitive and recovered a clear,
interpretable contrast:

- **Higher in HCM:** mitochondrial and energy-related programs —
  *oxidative phosphorylation*, *aerobic respiration*, *ATP synthesis*,
  *mitochondrial translation* — together with *ribosome biogenesis* and
  *cytoplasmic translation* (energy production + protein-synthesis machinery).
- **Higher in AS:** *extracellular matrix organization*, *cell junction
  assembly*, and a cluster of *neuronal / morphogenesis* terms — a structural,
  remodeling-oriented signature.

![GSEA dot plot](results/figures/10_enrichment/gsea_dotplot.png)

In short: relative to pressure-overload AS, HCM myocardium leans toward an
**energetic / mitochondrial and biosynthetic** program, while AS leans toward
**extracellular-matrix remodeling and structural** programs.

## Comparison with the original study

This reanalysis uses the same discovery cohort (8 HCM / 5 AS) as the source
study, so the results can be compared directly with its published RNA-seq
findings.

| | Original study (RNA-seq) | This reanalysis |
|---|---|---|
| Samples | 8 HCM / 5 AS | 8 HCM / 5 AS |
| DEGs (`padj < 0.05`) | 193 | 109 |
| Direction skew | 149/193 (77%) lower in HCM | 67/109 (61%) higher in AS (lower in HCM) |
| Strict (`|log2FC| > 1`) | 52 (38 down, 14 up) | 35 (20 higher in AS, 15 higher in HCM) |
| GO themes | locomotion, muscle structure development, neuron migration, cytoskeleton | ECM organization + neuronal/morphogenesis (AS); mitochondrial/energetic + translation (HCM) |

**Concordance of key genes.** Several genes highlighted by the original authors
were recovered here with the same direction of change:

| Gene | Original study | This reanalysis |
|---|---|---|
| `IGF2` | lower in HCM | −1.86 (higher in AS) |
| `C4B` | most significant, lower in HCM | among the top AS-side genes |
| `CTXND1` | higher in HCM | +1.49 (higher in HCM) |
| `EIF4EBP3` | higher in HCM | +1.24 (higher in HCM) |

**What agrees.** The direction of the overall signal (most differentially
expressed genes are lower in HCM / higher in AS), several individual marker
genes, and the neuronal/structural GO theme are all reproduced — strong evidence
that the pipeline captures a real biological signal rather than artefacts.

**What differs.** This reanalysis reports fewer DEGs (109 vs 193). This is
expected for an independent reanalysis with different software versions,
annotation, and filtering choices, and reflects a deliberately conservative
pipeline (group-aware low-count filtering plus DESeq2 independent filtering).

**What this analysis adds.** The original study used over-representation
analysis and emphasised the structural/neuronal themes. By additionally running
**GSEA** — which is more sensitive on subtle, coordinated signals — this
reanalysis surfaces a **mitochondrial / energetic and protein-synthesis program
elevated in HCM** that the ORA-only view does not capture. This is consistent
with the well-described energetic remodeling of HCM myocardium.

## Limitations

- small sample size
- comparison between two pathological states rather than disease versus healthy tissue
- biological conclusions should be interpreted as transcriptomic associations rather than causal mechanisms
- analysis starts from public processed gene-level counts rather than raw FASTQ files
- the differential-expression design models condition only (`~ condition_short`) and does not adjust for sex or age, which are available in the metadata; given the small sample size these potential confounders are not modelled, so some differences may partly reflect demographic imbalance between groups
- exploratory analysis (PCA, sample-distance heatmap) shows no clear global separation between HCM and AS, indicating the transcriptomic differences are subtle and confined to a limited set of genes rather than a genome-wide shift

## Conclusion

HCM and AS myocardium are globally similar at the transcriptome level — the
groups do not separate in unsupervised analysis — but they differ in a focused
set of genes and pathways. Relative to pressure-overload AS, HCM is associated
with an **energetic / mitochondrial and biosynthetic** program, while AS is
associated with **extracellular-matrix remodeling and structural** programs.
This independent reanalysis reproduces the key genes and the direction of the
published results, and adds an energetic signal in HCM via a more sensitive
enrichment method. Given the small sample size and the design, these findings
should be read as **transcriptomic associations** rather than causal mechanisms.