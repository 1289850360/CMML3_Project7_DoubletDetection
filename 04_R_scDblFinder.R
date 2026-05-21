install.packages("BiocManager")
options(repos = BiocManager::repositories())

BiocManager::install("bluster", ask = FALSE, update = FALSE)
BiocManager::install("scDblFinder", ask = FALSE, update = FALSE)

library(bluster)
library(scDblFinder)
data_dir <- "C:/Users/pjk/Documents/CMML3_Project7/GSE96583"
setwd(data_dir)

library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(dplyr)
library(readr)
library(pROC)
library(PRROC)

obj <- readRDS("results/seurat_object_DoubletFinder_only.rds")

# Use the labelled doublet rate as expected doublet rate
expected_rate <- mean(obj$ground_truth == "doublet")

sce <- as.SingleCellExperiment(obj)

set.seed(123)
sce <- scDblFinder(sce, dbr = expected_rate)

obj$scDblFinder_pred <- as.character(colData(sce)$scDblFinder.class)
obj$scDblFinder_score <- as.numeric(colData(sce)$scDblFinder.score)

calc_binary_metrics <- function(truth, pred, score = NULL) {
  truth_bin <- ifelse(truth == "doublet", 1, 0)
  pred_bin <- ifelse(pred == "doublet", 1, 0)
  
  tp <- sum(truth_bin == 1 & pred_bin == 1, na.rm = TRUE)
  fp <- sum(truth_bin == 0 & pred_bin == 1, na.rm = TRUE)
  fn <- sum(truth_bin == 1 & pred_bin == 0, na.rm = TRUE)
  tn <- sum(truth_bin == 0 & pred_bin == 0, na.rm = TRUE)
  
  precision <- ifelse(tp + fp == 0, NA, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, NA, tp / (tp + fn))
  f1 <- ifelse(precision + recall == 0, NA, 2 * precision * recall / (precision + recall))
  
  auroc <- NA
  auprc <- NA
  
  if (!is.null(score)) {
    valid <- !is.na(score) & !is.na(truth_bin)
    roc_obj <- pROC::roc(truth_bin[valid], score[valid], quiet = TRUE)
    auroc <- as.numeric(pROC::auc(roc_obj))
    
    pr <- PRROC::pr.curve(
      scores.class0 = score[valid][truth_bin[valid] == 1],
      scores.class1 = score[valid][truth_bin[valid] == 0],
      curve = FALSE
    )
    auprc <- pr$auc.integral
  }
  
  data.frame(
    Method = "scDblFinder",
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
  )
}

metrics_scdbl <- calc_binary_metrics(
  truth = obj$ground_truth,
  pred = obj$scDblFinder_pred,
  score = obj$scDblFinder_score
)

write.csv(metrics_scdbl, "results/metrics_scDblFinder.csv", row.names = FALSE)
saveRDS(obj, "results/seurat_object_with_scDblFinder.rds")

print(metrics_scdbl)