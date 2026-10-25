#!/usr/bin/env Rscript

###############################################################################
# NOTCH4-associated pathway analysis in HNSCC
# Updated reproducible analysis script for manuscript/reviewer response
# Version: 2.1 | Primary-tumour filtering corrected; mediation section retained unchanged
#
# Analyses:
#   1. TCGA-HNSC: primary-tumour expression, NOTCH4/pathway ssGSEA correlations
#   2. GSE65858: probe-to-gene mapping, ssGSEA and correlations
#   3. GSE39366: GPL9053 feature-to-gene mapping, ssGSEA and correlations
#   4. Combined validation table and figure
#   5. Optional NOTCH4 -> SFRP2 -> FN1 mediation model
#
# IMPORTANT:
# This script reconstructs the workflow from the analysis records available
# during manuscript preparation. It should be run end-to-end and the outputs
# compared with the final manuscript before being submitted as the original
# analysis code. The exact original script was not available for verification.
#
# GSE39366:
# The analysis previously used a processed file with columns FeatureNum and
# 138 HNxxx samples. The script downloads GEO supplementary files if possible,
# but if the processed file is not automatically retrieved, place the verified
# processed file at data/GSE39366.txt. Verify its provenance/sample mapping.
#
# TCGA:
# The correlation analyses use the 520 TCGA-HNSC primary tumour samples.
# The project-level cohort is HNSC and is not automatically an oral-cavity-only
# OSCC cohort. If the manuscript specifies a strict oral-cavity filter, document
# and implement that clinical rule separately before describing it as OSCC-only.
# The script selects samples from aligned TCGA metadata (sample_type), not from
# a loose barcode pattern, and stops if the expected 520 primary tumours are not
# found. Previously obtained correlation values are retained as diagnostics only;
# they are never substituted for calculated values.
###############################################################################

options(stringsAsFactors = FALSE)
set.seed(123)

# ------------------------------- Configuration -------------------------------

ROOT_DIR <- "."
DATA_DIR <- file.path(ROOT_DIR, "data")
RESULTS_DIR <- file.path(ROOT_DIR, "results")
TABLES_DIR <- file.path(RESULTS_DIR, "tables")
SCORES_DIR <- file.path(RESULTS_DIR, "scores")
RDS_DIR <- file.path(RESULTS_DIR, "rds")
LOG_DIR <- file.path(RESULTS_DIR, "logs")

for (d in c(DATA_DIR, RESULTS_DIR, TABLES_DIR, SCORES_DIR, RDS_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

FDR_THRESHOLD <- 0.05
LOG2FC_THRESHOLD <- 1
RUN_MEDIATION <- TRUE
MEDIATION_SIMS <- 1000
SPECIES <- "Homo sapiens"

# ------------------------------- Dependencies --------------------------------

cran_pkgs <- c("msigdbr", "mediation", "ggplot2")
bioc_pkgs <- c("TCGAbiolinks", "DESeq2", "GEOquery", "Biobase", "GSVA",
               "clusterProfiler", "SummarizedExperiment", "S4Vectors")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(DESeq2)
  library(GEOquery)
  library(Biobase)
  library(GSVA)
  library(msigdbr)
  library(clusterProfiler)
  library(mediation)
  library(ggplot2)
})

writeLines(capture.output(sessionInfo()),
           file.path(LOG_DIR, "sessionInfo_start.txt"))

# --------------------------------- Helpers -----------------------------------

write_csv <- function(x, filename, row.names = FALSE) {
  utils::write.csv(x, file.path(TABLES_DIR, filename), row.names = row.names)
}

# msigdbr API changed in recent releases. Try the current collection API first,
# then the legacy category/subcategory API.
get_msig_collection <- function(collection, subcollection = NULL) {
  f <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% f) {
    args <- list(species = SPECIES, collection = collection)
    if (!is.null(subcollection) && "subcollection" %in% f) {
      args$subcollection <- subcollection
    }
    return(do.call(msigdbr::msigdbr, args))
  }
  args <- list(species = SPECIES, category = collection)
  if (!is.null(subcollection)) args$subcategory <- subcollection
  do.call(msigdbr::msigdbr, args)
}

get_gene_sets <- function() {
  h <- get_msig_collection("H")
  if (!all(c("gs_name", "gene_symbol") %in% names(h))) {
    stop("Unexpected msigdbr columns; inspect package version and API.")
  }
  wnt <- unique(h$gene_symbol[h$gs_name == "HALLMARK_WNT_BETA_CATENIN_SIGNALING"])
  emt <- unique(h$gene_symbol[h$gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"])

  # NABA_CORE_MATRISOME is commonly distributed in C2:CP. If the installed
  # MSigDB version uses another collection, search the installed MSigDB tables.
  c2 <- tryCatch(get_msig_collection("C2", "CP"), error = function(e) NULL)
  mat <- character()
  if (!is.null(c2) && "gs_name" %in% names(c2)) {
    mat <- unique(c2$gene_symbol[c2$gs_name == "NABA_CORE_MATRISOME"])
  }
  if (!length(mat)) {
    all_c2 <- tryCatch(get_msig_collection("C2"), error = function(e) NULL)
    if (!is.null(all_c2) && "gs_name" %in% names(all_c2)) {
      mat <- unique(all_c2$gene_symbol[all_c2$gs_name == "NABA_CORE_MATRISOME"])
    }
  }
  if (!length(wnt) || !length(emt) || !length(mat)) {
    stop("Could not retrieve one or more required MSigDB gene sets. ",
         "Check MSigDB/msigdbr version and document the gene-set release.")
  }
  list(WNT_Beta_Catenin = wnt, EMT = emt, Core_Matrisome = mat)
}

collapse_probes <- function(expr, symbols) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  symbols <- as.character(symbols)
  keep <- !is.na(symbols) & nzchar(trimws(symbols))
  expr <- expr[keep, , drop = FALSE]
  symbols <- symbols[keep]
  sums <- rowsum(expr, group = symbols, reorder = FALSE)
  n <- table(symbols)
  sweep(sums, 1, as.numeric(n[rownames(sums)]), "/")
}

calculate_ssgsea <- function(expr, gene_sets, label) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  sets <- lapply(gene_sets, function(gs) intersect(unique(gs), rownames(expr)))
  coverage <- data.frame(
    Pathway = names(sets),
    Genes_in_set = vapply(sets, length, integer(1))
  )
  write_csv(coverage, paste0(label, "_pathway_gene_coverage.csv"))
  if (any(coverage$Genes_in_set < 10)) {
    warning(label, ": one or more gene sets have <10 genes represented.")
  }
  param <- GSVA::ssgseaParam(exprData = expr, geneSets = sets, normalize = TRUE)
  scores <- GSVA::gsva(param, verbose = FALSE)
  utils::write.csv(as.data.frame(scores),
                   file.path(SCORES_DIR, paste0(label, "_ssGSEA_scores.csv")),
                   row.names = TRUE)
  saveRDS(scores, file.path(RDS_DIR, paste0(label, "_ssGSEA_scores.rds")))
  scores
}

correlate_notch4 <- function(expr, scores, dataset_label) {
  if (!"NOTCH4" %in% rownames(expr)) {
    stop(dataset_label, ": NOTCH4 is not present after gene mapping.")
  }
  common <- intersect(colnames(expr), colnames(scores))
  if (length(common) < 10) stop(dataset_label, ": fewer than 10 matched samples.")
  x <- as.numeric(expr["NOTCH4", common])
  out <- lapply(rownames(scores), function(pathway) {
    y <- as.numeric(scores[pathway, common])
    z <- stats::cor.test(x, y, method = "spearman", exact = FALSE)
    data.frame(
      Dataset = dataset_label,
      Pathway = pathway,
      N = length(common),
      Spearman_rho = unname(z$estimate),
      P_value = z$p.value
    )
  })
  out <- do.call(rbind, out)
  out$FDR <- p.adjust(out$P_value, method = "BH")
  write_csv(out, paste0(gsub("[^A-Za-z0-9]+", "_", dataset_label),
                        "_NOTCH4_Pathway_Validation.csv"))
  out
}

# Reference values from the completed analysis in the project conversation.
# These are diagnostic comparison values only: the script NEVER overwrites
# newly calculated results with these values. Differences should be investigated.
reference_correlations <- data.frame(
  Dataset = rep(c("TCGA-HNSC", "GSE65858", "GSE39366"), each = 3),
  Pathway = rep(c("WNT_Beta_Catenin", "EMT", "Core_Matrisome"), 3),
  Reference_N = c(520, 520, 520, 270, 270, 270, 138, 138, 138),
  Reference_Spearman_rho = c(0.4526108, 0.4794599, 0.6365352,
                             0.4104579, 0.3656629, 0.5192685,
                             0.2507437, 0.3158276, 0.4875807),
  stringsAsFactors = FALSE
)

compare_with_reference <- function(observed) {
  observed_key <- paste(observed$Dataset, observed$Pathway, sep = "||")
  ref_key <- paste(reference_correlations$Dataset,
                   reference_correlations$Pathway, sep = "||")
  idx <- match(ref_key, observed_key)
  comparison <- reference_correlations
  comparison$Observed_N <- observed$N[idx]
  comparison$Observed_Spearman_rho <- observed$Spearman_rho[idx]
  comparison$Difference_rho <- comparison$Observed_Spearman_rho -
    comparison$Reference_Spearman_rho
  comparison$N_matches <- comparison$Observed_N == comparison$Reference_N
  comparison$rho_within_0.01 <- !is.na(comparison$Difference_rho) &
    abs(comparison$Difference_rho) <= 0.01
  write_csv(comparison, "Reference_vs_observed_validation_check.csv")

  message("\nReference comparison written to results/tables/Reference_vs_observed_validation_check.csv")
  if (any(is.na(comparison$Observed_Spearman_rho))) {
    warning("One or more reference dataset/pathway combinations are missing.")
  }
  if (any(!comparison$N_matches, na.rm = TRUE)) {
    warning("Observed sample counts differ from reference counts. Check filtering and sample matching.")
  }
  if (any(!comparison$rho_within_0.01, na.rm = TRUE)) {
    warning("Some correlations differ from reference values by >0.01. Check expression processing, gene sets, and GSVA versions; do not force values to match.")
  }
  invisible(comparison)
}

# --------------------------- MSigDB gene sets --------------------------------

gene_sets <- get_gene_sets()
saveRDS(gene_sets, file.path(RDS_DIR, "MSigDB_gene_sets_used.rds"))
writeLines(c(
  paste("WNT genes in source set:", length(gene_sets$WNT_Beta_Catenin)),
  paste("EMT genes in source set:", length(gene_sets$EMT)),
  paste("Core matrisome genes in source set:", length(gene_sets$Core_Matrisome)),
  paste("msigdbr version:", as.character(packageVersion("msigdbr")))
), file.path(LOG_DIR, "gene_set_metadata.txt"))

# ------------------------------- TCGA-HNSC -----------------------------------
# Download/prepare TCGA-HNSC primary tumours and normal samples. Align metadata
# explicitly by sample barcode, then retain only Primary Tumor samples.
# This analysis uses all TCGA-HNSC primary tumours, not an oral-cavity-only subset.

tcga_rds <- file.path(RDS_DIR, "TCGA_HNSC_se.rds")
if (file.exists(tcga_rds)) {
  tcga_se <- readRDS(tcga_rds)
} else {
  query <- TCGAbiolinks::GDCquery(
    project = "TCGA-HNSC",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
  )
  TCGAbiolinks::GDCdownload(query, directory = DATA_DIR)
  tcga_se <- TCGAbiolinks::GDCprepare(query, directory = DATA_DIR)
  saveRDS(tcga_se, tcga_rds)
}

count_assay <- SummarizedExperiment::assayNames(tcga_se)
assay_name <- if ("unstranded" %in% count_assay) "unstranded" else count_assay[1]
tcga_counts_all <- SummarizedExperiment::assay(tcga_se, assay_name)
tcga_meta <- as.data.frame(SummarizedExperiment::colData(tcga_se))
tcga_rd <- as.data.frame(SummarizedExperiment::rowData(tcga_se))

# Align metadata to the count matrix by the full sample barcode.
if (!"barcode" %in% names(tcga_meta)) {
  stop("TCGA metadata has no 'barcode' column; inspect colData(tcga_se).")
}
meta_idx <- match(colnames(tcga_counts_all), as.character(tcga_meta$barcode))
if (anyNA(meta_idx)) {
  stop("Some TCGA count columns could not be matched to metadata barcodes.")
}
tcga_meta_aligned <- tcga_meta[meta_idx, , drop = FALSE]
if (!identical(as.character(tcga_meta_aligned$barcode), colnames(tcga_counts_all))) {
  stop("TCGA metadata alignment failed; count columns and metadata barcodes differ.")
}

if (!"sample_type" %in% names(tcga_meta_aligned)) {
  stop("TCGA metadata has no 'sample_type' column; cannot select primary tumours safely.")
}
primary <- !is.na(tcga_meta_aligned$sample_type) &
  tcga_meta_aligned$sample_type == "Primary Tumor"
if (sum(primary) != 520L) {
  stop("Expected 520 TCGA-HNSC Primary Tumor samples, found ", sum(primary),
       ". Check the downloaded cohort and sample inclusion criteria before proceeding.")
}
tcga_counts_primary_raw <- tcga_counts_all[, primary, drop = FALSE]
message("TCGA primary tumour samples selected using metadata: ", ncol(tcga_counts_primary_raw))

# Map Ensembl-feature rows to gene symbols and collapse duplicate symbols by sum.
gene_col <- intersect(c("gene_name", "external_gene_name", "gene"), names(tcga_rd))[1]
if (is.na(gene_col)) {
  stop("Could not find gene symbol column in TCGA rowData. Inspect rowData(tcga_se).")
}
gene_symbols <- as.character(tcga_rd[[gene_col]])
keep_genes <- !is.na(gene_symbols) & nzchar(gene_symbols) & !grepl("^N_", gene_symbols)
tcga_counts_primary_raw <- tcga_counts_primary_raw[keep_genes, , drop = FALSE]
gene_symbols <- gene_symbols[keep_genes]
tcga_counts_primary <- rowsum(tcga_counts_primary_raw, group = gene_symbols, reorder = FALSE)
storage.mode(tcga_counts_primary) <- "integer"

if (!"NOTCH4" %in% rownames(tcga_counts_primary)) {
  stop("NOTCH4 absent from TCGA counts after gene-symbol mapping.")
}

# Normalize counts on the 520 primary tumour samples only.
dds_norm <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(tcga_counts_primary),
  colData = S4Vectors::DataFrame(row.names = colnames(tcga_counts_primary)),
  design = ~ 1
)
dds_norm <- DESeq2::estimateSizeFactors(dds_norm)
norm_counts <- DESeq2::counts(dds_norm, normalized = TRUE)
tcga_expr <- log2(norm_counts + 1)
tcga_expr_primary <- tcga_expr
saveRDS(tcga_counts_primary, file.path(RDS_DIR, "TCGA_HNSC_primary_tumour_counts.rds"))
saveRDS(norm_counts, file.path(RDS_DIR, "TCGA_HNSC_normalized_counts.rds"))
saveRDS(tcga_expr, file.path(RDS_DIR, "TCGA_HNSC_log2_normalized_expression.rds"))
saveRDS(tcga_expr_primary, file.path(RDS_DIR, "TCGA_HNSC_primary_tumour_expression.rds"))

message("TCGA expression matrix dimensions (genes x primary tumours): ",
        nrow(tcga_expr_primary), " x ", ncol(tcga_expr_primary))

scores_tcga <- calculate_ssgsea(tcga_expr_primary, gene_sets, "TCGA_HNSC")
cor_tcga <- correlate_notch4(tcga_expr_primary, scores_tcga, "TCGA-HNSC")

# NOTCH4 high/low DESeq2 analysis using median split.
notch4 <- as.numeric(tcga_expr["NOTCH4", ])
group <- factor(ifelse(notch4 >= median(notch4, na.rm = TRUE), "High", "Low"),
                levels = c("Low", "High"))
coldata_de <- S4Vectors::DataFrame(group = group,
                                   row.names = colnames(tcga_counts_primary))
dds_de <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(tcga_counts_primary),
  colData = coldata_de,
  design = ~ group
)
dds_de <- dds_de[rowSums(DESeq2::counts(dds_de)) >= 10, ]
dds_de <- DESeq2::DESeq(dds_de)
de <- as.data.frame(DESeq2::results(dds_de, contrast = c("group", "High", "Low")))
de$Gene <- rownames(de)
de <- de[order(de$padj, na.last = TRUE), ]
write_csv(de, "TCGA_NOTCH4_High_vs_Low_DESeq2.csv")
write_csv(subset(de, !is.na(padj) & padj < FDR_THRESHOLD),
          "TCGA_NOTCH4_High_vs_Low_FDR05.csv")
write_csv(subset(de, !is.na(padj) & padj < FDR_THRESHOLD &
                   abs(log2FoldChange) >= LOG2FC_THRESHOLD),
          "TCGA_NOTCH4_High_vs_Low_FDR05_absLFC1.csv")

# Hallmark GSEA using DESeq2 Wald statistic.
ranked <- de$stat
names(ranked) <- de$Gene
ranked <- ranked[!is.na(ranked) & !duplicated(names(ranked))]
ranked <- sort(ranked, decreasing = TRUE)
hallmark <- get_msig_collection("H")
term2gene <- unique(hallmark[, c("gs_name", "gene_symbol")])
colnames(term2gene) <- c("term", "gene")
gsea <- clusterProfiler::GSEA(
  geneList = ranked,
  TERM2GENE = term2gene,
  pvalueCutoff = 1,
  verbose = FALSE
)
write_csv(as.data.frame(gsea), "TCGA_Hallmark_GSEA.csv")

# Optional mediation model on TCGA expression. This section is retained unchanged per project decision.
if (RUN_MEDIATION && all(c("NOTCH4", "SFRP2", "FN1") %in% rownames(tcga_expr))) {
  med_df <- data.frame(
    NOTCH4 = as.numeric(tcga_expr["NOTCH4", ]),
    SFRP2 = as.numeric(tcga_expr["SFRP2", ]),
    FN1 = as.numeric(tcga_expr["FN1", ])
  )
  model_m <- stats::lm(SFRP2 ~ NOTCH4, data = med_df)
  model_y <- stats::lm(FN1 ~ NOTCH4 + SFRP2, data = med_df)
  set.seed(123)
  med_fit <- mediation::mediate(
    model.m = model_m, model.y = model_y,
    treat = "NOTCH4", mediator = "SFRP2",
    boot = TRUE, sims = MEDIATION_SIMS
  )
  writeLines(capture.output(summary(med_fit)),
             file.path(TABLES_DIR, "Mediation_NOTCH4_SFRP2_FN1.txt"))
  saveRDS(med_fit, file.path(RDS_DIR, "Mediation_NOTCH4_SFRP2_FN1.rds"))
}

# ------------------------------ GSE65858 --------------------------------------
# GEOquery obtains the series matrix and platform annotation. The platform
# annotation is mapped by probe ID, then duplicate probes are averaged by gene.

gse65858 <- GEOquery::getGEO("GSE65858", GSEMatrix = TRUE, getGPL = FALSE)
if (!length(gse65858)) stop("GSE65858 download returned no expression set.")
eset65858 <- gse65858[[1]]
expr65858_probe <- Biobase::exprs(eset65858)
message("GSE65858 samples downloaded: ", ncol(expr65858_probe))
if (ncol(expr65858_probe) != 270) {
  warning("GSE65858 sample count differs from the 270 samples used in the completed analysis. Check the selected series matrix/subset.")
}
platform65858 <- Biobase::annotation(eset65858)
gpl65858 <- GEOquery::getGEO(platform65858)
tab65858 <- GEOquery::Table(gpl65858)

id_col <- intersect(c("ID", "ID_REF"), names(tab65858))[1]
symbol_col <- intersect(c("Symbol", "Gene Symbol", "GENE_SYMBOL", "Gene symbol"),
                        names(tab65858))[1]
if (is.na(id_col) || is.na(symbol_col)) {
  stop("Could not identify probe ID and symbol columns for ", platform65858,
       ". Inspect colnames(GEOquery::Table(gpl65858)).")
}
map_idx <- match(rownames(expr65858_probe), as.character(tab65858[[id_col]]))
symbols65858 <- as.character(tab65858[[symbol_col]][map_idx])
symbols65858 <- sub("\\s*///.*$", "", symbols65858)
symbols65858 <- sub("\\s*//.*$", "", symbols65858)
expr65858_gene <- collapse_probes(expr65858_probe, symbols65858)
saveRDS(expr65858_gene, file.path(RDS_DIR, "GSE65858_gene_expression.rds"))
scores65858 <- calculate_ssgsea(expr65858_gene, gene_sets, "GSE65858")
cor65858 <- correlate_notch4(expr65858_gene, scores65858, "GSE65858")

# ------------------------------ GSE39366 --------------------------------------
# First attempt: use an existing processed file. Otherwise attempt GEO
# supplementary download and stop if the expected processed matrix is absent.

gse39366_file <- file.path(DATA_DIR, "GSE39366.txt")
if (!file.exists(gse39366_file)) {
  try(GEOquery::getGEOSuppFiles("GSE39366", baseDir = DATA_DIR,
                               makeDirectory = TRUE), silent = TRUE)
  candidate_files <- list.files(
    file.path(DATA_DIR, "GSE39366"),
    recursive = TRUE, full.names = TRUE
  )
  candidate_files <- candidate_files[grepl("GSE39366.*\\.txt(\\.gz)?$",
                                            basename(candidate_files),
                                            ignore.case = TRUE)]
  if (length(candidate_files) != 1) {
    stop("Could not identify exactly one processed GSE39366 file. ",
         "After verifying provenance, place it at data/GSE39366.txt.")
  }
  gse39366_file <- candidate_files[1]
}

message("Using GSE39366 processed expression file: ", normalizePath(gse39366_file, mustWork = TRUE))
expr39366_full <- read.delim(
  gse39366_file, header = TRUE, sep = "\t",
  check.names = FALSE, stringsAsFactors = FALSE
)
if (!"FeatureNum" %in% names(expr39366_full)) {
  stop("GSE39366 file lacks FeatureNum. Verify that the correct processed file was used.")
}
expr39366_probe <- as.matrix(expr39366_full[, setdiff(names(expr39366_full), "FeatureNum"),
                                             drop = FALSE])
rownames(expr39366_probe) <- as.character(expr39366_full$FeatureNum)
storage.mode(expr39366_probe) <- "numeric"
message("GSE39366 samples read from processed file: ", ncol(expr39366_probe))
if (ncol(expr39366_probe) != 138) {
  warning("GSE39366 sample count differs from the 138 samples used in the completed analysis. Verify the processed file and sample columns.")
}

gpl39366 <- GEOquery::getGEO("GPL9053")
tab39366 <- GEOquery::Table(gpl39366)
if (!all(c("ID", "ORF") %in% names(tab39366))) {
  stop("GPL9053 annotation does not contain expected ID and ORF fields.")
}
idx39366 <- match(rownames(expr39366_probe), as.character(tab39366$ID))
symbols39366 <- as.character(tab39366$ORF[idx39366])
expr39366_gene <- collapse_probes(expr39366_probe, symbols39366)
saveRDS(expr39366_gene, file.path(RDS_DIR, "GSE39366_gene_expression.rds"))
scores39366 <- calculate_ssgsea(expr39366_gene, gene_sets, "GSE39366")
cor39366 <- correlate_notch4(expr39366_gene, scores39366, "GSE39366")

# -------------------------- Combined validation ------------------------------

all_cor <- rbind(cor_tcga, cor65858, cor39366)
write_csv(all_cor, "Three_Dataset_NOTCH4_Pathway_Validation_with_Pvalues_FDR.csv")
reference_comparison <- compare_with_reference(all_cor)
write_csv(all_cor[, c("Dataset", "Pathway", "N", "Spearman_rho", "P_value", "FDR")],
          "Three_Dataset_NOTCH4_Validation.csv")
# Explicitly named final table for manuscript/reviewer checks.
write_csv(all_cor[, c("Dataset", "Pathway", "N", "Spearman_rho", "P_value", "FDR")],
          "Three_Dataset_NOTCH4_Validation_FINAL.csv")

plot_df <- all_cor
plot_df$Dataset <- factor(plot_df$Dataset,
                          levels = c("TCGA-HNSC", "GSE65858", "GSE39366"))
plot_df$Pathway <- factor(plot_df$Pathway,
                          levels = c("WNT_Beta_Catenin", "EMT", "Core_Matrisome"))
fig <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = Dataset, y = Spearman_rho,
                                    group = Pathway)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::facet_wrap(~ Pathway, nrow = 1) +
  ggplot2::coord_cartesian(ylim = c(-1, 1)) +
  ggplot2::labs(
    title = "NOTCH4 associations with pathway scores across datasets",
    x = "Dataset", y = "Spearman correlation (rho)"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
ggplot2::ggsave(file.path(RESULTS_DIR, "Figure_NOTCH4_Pathway_Validation.png"),
                fig, width = 11, height = 4.5, dpi = 300)
ggplot2::ggsave(file.path(RESULTS_DIR, "Figure_NOTCH4_Pathway_Validation.pdf"),
                fig, width = 11, height = 4.5)

writeLines(capture.output(sessionInfo()),
           file.path(LOG_DIR, "sessionInfo_end.txt"))
writeLines(c(
  paste("Completed at:", Sys.time()),
  paste("TCGA primary tumour samples:", ncol(tcga_expr_primary)),
  paste("GSE65858 samples:", ncol(expr65858_gene)),
  paste("GSE39366 samples:", ncol(expr39366_gene)),
  "Reference correlations are diagnostics only; calculated values are not overwritten.",
  paste("R version:", R.version.string)
), file.path(LOG_DIR, "run_summary.txt"))

message("Completed. Review outputs in ", normalizePath(RESULTS_DIR))
