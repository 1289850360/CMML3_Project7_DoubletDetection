############################################################
# CMML3 Project 7
# Benchmarking doublet detection methods
# Dataset: GSE96583 batch 2, sample 2.1 / 2.2
# This script runs:
# 1. metadata cleaning
# 2. Seurat preprocessing
# 3. DoubletFinder
# 4. exports count matrix for Scrublet
############################################################

############################
# 0. Set working directory
############################

data_dir <- "C:/Users/pjk/Documents/CMML3_Project7/GSE96583"
setwd(data_dir)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

############################
# 1. Install/load packages
############################

cran_pkgs <- c(
  "Seurat", "Matrix", "dplyr", "tidyr", "stringr",
  "ggplot2", "patchwork", "readr", "remotes", "pROC", "PRROC"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
}

library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(readr)
library(DoubletFinder)
library(pROC)
library(PRROC)

set.seed(123)

############################
# 2. Read and clean metadata
############################

meta2_raw <- read_tsv(
  "GSE96583_batch2.total.tsne.df.tsv.gz",
  show_col_types = FALSE
)

# The first column named tsne1 is actually barcode.
# Column names are shifted, so we rename them manually.
meta2 <- meta2_raw %>%
  rename(
    barcode = tsne1,
    tsne1 = tsne2,
    tsne2 = ind,
    ind = stim,
    condition = cluster,
    cluster = cell,
    cell_label_raw = multiplets
  ) %>%
  separate(
    cell_label_raw,
    into = c("celltype", "ground_truth"),
    sep = "\t",
    remove = FALSE,
    fill = "right"
  ) %>%
  mutate(
    barcode = trimws(str_replace_all(as.character(barcode), '"', "")),
    celltype = trimws(str_replace_all(as.character(celltype), '"', "")),
    ground_truth = trimws(str_replace_all(as.character(ground_truth), '"', "")),
    condition = as.character(condition)
  )

cat("\nMetadata dimensions:\n")
print(dim(meta2))

cat("\nFirst metadata barcodes:\n")
print(head(meta2$barcode, 10))

cat("\nGround truth labels:\n")
print(table(meta2$ground_truth, useNA = "ifany"))

cat("\nCondition labels:\n")
print(table(meta2$condition, useNA = "ifany"))

############################
# 3. Read genes
############################

genes_raw <- read_tsv(
  "GSE96583_batch2.genes.tsv.gz",
  col_names = FALSE,
  show_col_types = FALSE
)

if (ncol(genes_raw) >= 2) {
  gene_names <- genes_raw[[2]]
} else {
  gene_names <- genes_raw[[1]]
}

gene_names <- as.character(gene_names)

bad_gene <- is.na(gene_names) | gene_names == "" | gene_names == "NA"
gene_names[bad_gene] <- paste0("Gene_", which(bad_gene))

gene_names <- make.unique(gene_names)

cat("\nNumber of genes:\n")
print(length(gene_names))

cat("\nFirst genes:\n")
print(head(gene_names))

############################
# 4. Read one matrix safely
############################

read_one_sample <- function(matrix_file, barcode_file, sample_name) {
  cat("\nReading sample:", sample_name, "\n")
  
  mat <- readMM(gzfile(matrix_file))
  bc <- readLines(gzfile(barcode_file))
  bc <- trimws(bc)
  bc <- bc[bc != ""]
  
  cat("Raw matrix dim:\n")
  print(dim(mat))
  cat("Number of barcodes:\n")
  print(length(bc))
  
  if (nrow(mat) == length(gene_names) && ncol(mat) == length(bc)) {
    rownames(mat) <- gene_names
    colnames(mat) <- bc
  } else if (ncol(mat) == length(gene_names) && nrow(mat) == length(bc)) {
    mat <- t(mat)
    rownames(mat) <- gene_names
    colnames(mat) <- bc
  } else {
    stop("Matrix dimensions do not match genes/barcodes for sample ", sample_name)
  }
  
  matched <- sum(colnames(mat) %in% meta2$barcode)
  cat("Matched barcodes with metadata:", matched, "\n")
  
  return(list(counts = mat, matched = matched, sample = sample_name))
}

sample_21 <- read_one_sample(
  matrix_file = "GSM2560248_2.1.mtx.gz",
  barcode_file = "GSM2560248_barcodes.tsv.gz",
  sample_name = "2.1"
)

sample_22 <- read_one_sample(
  matrix_file = "GSM2560249_2.2.mtx.gz",
  barcode_file = "GSM2560249_barcodes.tsv.gz",
  sample_name = "2.2"
)

# Select the sample with more matched barcodes.
if (sample_21$matched >= sample_22$matched) {
  counts <- sample_21$counts
  selected_sample <- "2.1"
} else {
  counts <- sample_22$counts
  selected_sample <- "2.2"
}

cat("\nSelected sample:\n")
print(selected_sample)

############################
# 5. Match metadata and counts
############################

meta_current <- meta2[meta2$barcode %in% colnames(counts), ]
meta_current <- as.data.frame(meta_current)

cat("\nMetadata rows matching selected matrix:\n")
print(nrow(meta_current))

cat("\nLabels before removing ambs:\n")
print(table(meta_current$ground_truth, useNA = "ifany"))

# Match metadata to count matrix order
idx <- match(colnames(counts), meta_current$barcode)
keep <- !is.na(idx)

counts <- counts[, keep, drop = FALSE]
meta_current <- meta_current[idx[keep], ]

# Remove ambiguous cells labelled "ambs"; keep only singlet and doublet
keep_label <- meta_current$ground_truth %in% c("singlet", "doublet")

counts <- counts[, keep_label, drop = FALSE]
meta_ctrl <- meta_current[keep_label, ]

rownames(meta_ctrl) <- meta_ctrl$barcode

cat("\nFinal matched labels:\n")
print(table(meta_ctrl$ground_truth, useNA = "ifany"))

cat("\nFinal counts dimension:\n")
print(dim(counts))

cat("\nFirst final count barcodes:\n")
print(head(colnames(counts), 5))

cat("\nFirst final metadata barcodes:\n")
print(head(rownames(meta_ctrl), 5))

if (!all(colnames(counts) == rownames(meta_ctrl))) {
  stop("counts colnames and metadata rownames are not aligned")
}

############################
# 6. Create Seurat object
############################

obj <- CreateSeuratObject(
  counts = counts,
  meta.data = meta_ctrl,
  project = "GSE96583_batch2",
  min.cells = 3,
  min.features = 200
)

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

cat("\nLabels after CreateSeuratObject:\n")
print(table(obj$ground_truth, useNA = "ifany"))

############################
# 7. QC plots and filtering
############################

p_qc_before <- VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.05
)

ggsave(
  "figures/Supp_QC_before_filtering.png",
  p_qc_before,
  width = 10,
  height = 4,
  dpi = 300
)

# Mild QC filtering
obj <- subset(
  obj,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 15
)

cat("\nLabels after QC filtering:\n")
print(table(obj$ground_truth, useNA = "ifany"))

p_qc_after <- VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.05
)

ggsave(
  "figures/Supp_QC_after_filtering.png",
  p_qc_after,
  width = 10,
  height = 4,
  dpi = 300
)

############################
# 8. Standard Seurat workflow
############################

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = rownames(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30)

p_elbow <- ElbowPlot(obj, ndims = 30)
ggsave("figures/Supp_PCA_elbow.png", p_elbow, width = 5, height = 4, dpi = 300)

dims_use <- 1:20

obj <- FindNeighbors(obj, dims = dims_use)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = dims_use)

p_truth <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "ground_truth",
  pt.size = 0.2
) + ggtitle("Ground truth")

ggsave("figures/UMAP_ground_truth.png", p_truth, width = 5, height = 4, dpi = 300)

############################
# 9. Metric function
############################

calc_binary_metrics <- function(truth, pred, score = NULL) {
  truth_bin <- ifelse(truth == "doublet", 1, 0)
  pred_bin <- ifelse(pred == "doublet", 1, 0)
  
  tp <- sum(truth_bin == 1 & pred_bin == 1, na.rm = TRUE)
  fp <- sum(truth_bin == 0 & pred_bin == 1, na.rm = TRUE)
  fn <- sum(truth_bin == 1 & pred_bin == 0, na.rm = TRUE)
  tn <- sum(truth_bin == 0 & pred_bin == 0, na.rm = TRUE)
  
  precision <- ifelse(tp + fp == 0, NA, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, NA, tp / (tp + fn))
  f1 <- ifelse(
    is.na(precision + recall) || precision + recall == 0,
    NA,
    2 * precision * recall / (precision + recall)
  )
  
  auroc <- NA
  auprc <- NA
  
  if (!is.null(score)) {
    score <- as.numeric(score)
    valid <- !is.na(score) & !is.na(truth_bin)
    
    if (length(unique(truth_bin[valid])) == 2) {
      roc_obj <- pROC::roc(truth_bin[valid], score[valid], quiet = TRUE)
      auroc <- as.numeric(pROC::auc(roc_obj))
      
      pr <- PRROC::pr.curve(
        scores.class0 = score[valid][truth_bin[valid] == 1],
        scores.class1 = score[valid][truth_bin[valid] == 0],
        curve = FALSE
      )
      auprc <- pr$auc.integral
    }
  }
  
  return(data.frame(
    TP = tp,
    FP = fp,
    FN = fn,
    TN = tn,
    Precision = precision,
    Recall = recall,
    F1 = f1,
    AUROC = auroc,
    AUPRC = auprc,
    Predicted_doublet_rate = mean(pred == "doublet", na.rm = TRUE)
  ))
}

############################
# 10. Run DoubletFinder
############################

cat("\nRunning DoubletFinder...\n")

expected_rate <- mean(obj$ground_truth == "doublet")
nExp <- round(expected_rate * ncol(obj))

cat("\nExpected doublet rate:\n")
print(expected_rate)

cat("\nExpected doublet number:\n")
print(nExp)

# Param sweep can be slow but helps choose pK
if (exists("paramSweep")) {
  sweep.res <- paramSweep(obj, PCs = dims_use, sct = FALSE)
} else {
  sweep.res <- paramSweep_v3(obj, PCs = dims_use, sct = FALSE)
}

sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

best_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

cat("\nBest pK selected by BCMVN:\n")
print(best_pK)

if (exists("doubletFinder")) {
  obj <- doubletFinder(
    obj,
    PCs = dims_use,
    pN = 0.25,
    pK = best_pK,
    nExp = nExp,
    reuse.pANN = NULL,
    sct = FALSE
  )
} else {
  obj <- doubletFinder_v3(
    obj,
    PCs = dims_use,
    pN = 0.25,
    pK = best_pK,
    nExp = nExp,
    reuse.pANN = NULL,
    sct = FALSE
  )
}

df_class_col <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
df_score_col <- grep("pANN", colnames(obj@meta.data), value = TRUE)

df_class_col <- df_class_col[length(df_class_col)]
df_score_col <- df_score_col[length(df_score_col)]

obj$DoubletFinder_pred <- ifelse(
  obj@meta.data[[df_class_col]] == "Doublet",
  "doublet",
  "singlet"
)

obj$DoubletFinder_score <- obj@meta.data[[df_score_col]]

cat("\nDoubletFinder predictions:\n")
print(table(obj$DoubletFinder_pred, useNA = "ifany"))

p_df <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "DoubletFinder_pred",
  pt.size = 0.2
) + ggtitle("DoubletFinder prediction")

ggsave("figures/UMAP_DoubletFinder.png", p_df, width = 5, height = 4, dpi = 300)

metrics_df <- calc_binary_metrics(
  truth = obj$ground_truth,
  pred = obj$DoubletFinder_pred,
  score = obj$DoubletFinder_score
)

metrics_df$Method <- "DoubletFinder"
metrics_df <- metrics_df %>% select(Method, everything())

write.csv(metrics_df, "results/metrics_DoubletFinder.csv", row.names = FALSE)

cat("\nDoubletFinder metrics:\n")
print(metrics_df)

############################
# 11. Export files for Scrublet
############################

cat("\nSaving files for Scrublet...\n")

counts_for_scrublet <- GetAssayData(obj, assay = "RNA", slot = "counts")

writeMM(counts_for_scrublet, "results/counts_for_scrublet.mtx")

scrublet_meta <- data.frame(
  barcode = colnames(obj),
  ground_truth = obj$ground_truth,
  celltype = obj$celltype,
  condition = obj$condition,
  seurat_cluster = as.character(Idents(obj)),
  UMAP_1 = Embeddings(obj, "umap")[, 1],
  UMAP_2 = Embeddings(obj, "umap")[, 2],
  DoubletFinder_pred = obj$DoubletFinder_pred,
  DoubletFinder_score = obj$DoubletFinder_score
)

write.csv(scrublet_meta, "results/metadata_for_scrublet.csv", row.names = FALSE)

saveRDS(obj, "results/seurat_object_DoubletFinder_only.rds")

cat("\nFinished script 01 successfully.\n")
cat("Generated files:\n")
cat("results/metrics_DoubletFinder.csv\n")
cat("results/metadata_for_scrublet.csv\n")
cat("results/counts_for_scrublet.mtx\n")
cat("results/seurat_object_DoubletFinder_only.rds\n")
cat("figures/UMAP_ground_truth.png\n")
cat("figures/UMAP_DoubletFinder.png\n")
