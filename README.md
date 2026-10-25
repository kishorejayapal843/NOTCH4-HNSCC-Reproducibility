# NOTCH4-associated pathway analysis in HNSCC

## Overview
This repository contains an R workflow for analysing NOTCH4-associated pathway scores in TCGA-HNSC and two GEO datasets, GSE65858 and GSE39366. The workflow includes differential expression, Hallmark gene-set enrichment, ssGSEA pathway scoring, Spearman correlations, a combined validation table, and an optional mediation analysis.

**Important:** This is a consolidated/reconstructed workflow prepared from the analysis steps and results available during manuscript preparation. It should not be described as an exact reproduction of the original analysis until all input-data provenance, preprocessing, package versions, and results have been verified. The script writes calculated values; it does not replace them with reference values.

## Datasets
- **TCGA-HNSC** — downloaded using `TCGAbiolinks`. The pathway-correlation analysis selects primary tumour samples using sample metadata and expects 520 samples.
- **GSE65858** — retrieved through GEOquery.
- **GSE39366** — uses a processed expression file with a `FeatureNum` column and sample columns. Confirm that this is the same processed file and sample mapping used for the analysis. If needed, place the verified file at `data/GSE39366.txt`.

Use public accession records and their associated publications to verify the source and processing of each dataset.

## Analysis outline
1. Install/load R and Bioconductor dependencies.
2. Retrieve the gene sets from `msigdbr`/MSigDB.
3. Download and preprocess TCGA-HNSC counts; select primary tumours for pathway correlations.
4. Perform NOTCH4 high-versus-low differential expression with DESeq2.
5. Run Hallmark GSEA using the DESeq2 test statistic.
6. Process GSE65858 probe annotations and collapse probes mapped to the same gene.
7. Map GSE39366 feature IDs using GPL9053 annotation and collapse probes mapped to the same gene.
8. Calculate ssGSEA pathway scores for WNT/β-catenin, EMT, and core matrisome.
9. Calculate Spearman correlations between NOTCH4 expression and pathway scores.
10. Export tables, pathway scores, RDS objects, figures, and session logs.

## Run the script
Open RStudio with this repository root as the working directory, then run:

```r
source("scripts/NOTCH4_HNSCC_reproducible_analysis_GITHUB.R")
```

The script creates `data/` and `results/` subfolders if they do not exist. Internet access is required for package installation and public dataset retrieval. For a fully reproducible submission, record the exact R and package versions and consider using `renv` to lock the environment.

## Outputs
Generated output folders include:
- `results/tables/` — differential expression, GSEA, mediation summary, pathway coverage, validation correlations, and reference comparison.
- `results/scores/` — ssGSEA score matrices.
- `results/rds/` — processed expression matrices, gene sets, scores, and model objects.
- `results/logs/` — R session information and run summary.
- `results/` — validation figures.

## Validation status and interpretation
The workflow was run in R 4.5.2 in the current analysis environment. The reported sample counts were TCGA-HNSC primary tumours = 520, GSE65858 = 270, and GSE39366 = 138.

The previously obtained reference correlations were:
- TCGA-HNSC: WNT/β-catenin ρ = 0.4526; EMT ρ = 0.4795; core matrisome ρ = 0.6365 (n = 520).
- GSE65858: WNT/β-catenin ρ = 0.4105; EMT ρ = 0.3657; core matrisome ρ = 0.5193 (n = 270).
- GSE39366: WNT/β-catenin ρ = 0.2507; EMT ρ = 0.3158; core matrisome ρ = 0.4876 (n = 138).

The latest master-script output produced TCGA-HNSC correlations of approximately 0.4572, 0.4867, and 0.6127 for WNT/β-catenin, EMT, and core matrisome, respectively, with n = 520. The GEO correlations matched the reference values to the displayed precision. The TCGA correlations therefore do not exactly reproduce the earlier reference values; this discrepancy should be investigated or clearly documented before describing the workflow as an exact reproduction.

## Mediation analysis
The mediation section is retained as requested. Its output is generated from the model and settings in the script. Do not imply that mediation establishes causality without appropriate design assumptions and caveats.

## Suggested citation and data-use notes
Please cite the original TCGA and GEO data resources and the publications associated with each dataset in any manuscript using these analyses. Confirm each dataset's licence and use conditions before redistribution. Do not commit sensitive or restricted data to this repository.
