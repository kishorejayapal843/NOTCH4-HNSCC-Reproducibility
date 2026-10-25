# Input data

TCGA-HNSC and GSE65858 are retrieved by the R script from public resources where possible.

For GSE39366, place the verified processed expression file at `data/GSE39366.txt` if the script cannot identify/download the exact expected supplementary file automatically. The file is expected to contain a `FeatureNum` column plus sample expression columns. Confirm its provenance and sample mapping against the data used for the analysis before submission.

Do not add large or restricted data files to GitHub unnecessarily. Record public accessions, source links, retrieval dates, and any preprocessing details in the main README or a supplementary data manifest.
