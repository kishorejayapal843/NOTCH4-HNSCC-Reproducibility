# Input Data

The analysis uses TCGA-HNSC and two Gene Expression Omnibus (GEO) datasets: GSE65858 and GSE39366.

TCGA-HNSC and GSE65858 are retrieved by the R script from public resources where possible.

For GSE39366, place the verified processed expression file at `data/GSE39366.txt` if the script cannot identify or download the expected supplementary file automatically. The file is expected to contain a `FeatureNum` column and sample-expression columns.

Verify the source, sample mapping, and preprocessing of each dataset before using the results in a manuscript.

Do not upload confidential, identifiable, or restricted data to this repository.
