############################################################
# Redraw final figures with scDblFinder included
############################################################

data_dir <- "C:/Users/pjk/Documents/CMML3_Project7/GSE96583"
setwd(data_dir)

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

############################
# 1. Read metrics
############################

metrics_df <- read_csv("results/metrics_DoubletFinder.csv", show_col_types = FALSE)
metrics_scr <- read_csv("results/metrics_scrublet.csv", show_col_types = FALSE)
metrics_scdbl <- read_csv("results/metrics_scDblFinder.csv", show_col_types = FALSE)

metrics_all <- bind_rows(metrics_df, metrics_scr, metrics_scdbl) %>%
  select(
    Method, TP, FP, FN, TN,
    Precision, Recall, F1, AUROC, AUPRC,
    Predicted_doublet_rate,
    everything()
  )

write.csv(metrics_all, "results/final_benchmark_metrics_3methods.csv", row.names = FALSE)

print(metrics_all)

############################
# 2. Draw metric bar plot (3 methods)
############################

metrics_long <- metrics_all %>%
  select(Method, Precision, Recall, F1, AUROC, AUPRC) %>%
  pivot_longer(
    cols = c(Precision, Recall, F1, AUROC, AUPRC),
    names_to = "Metric",
    values_to = "Value"
  )

p_metrics <- ggplot(metrics_long, aes(x = Method, y = Value)) +
  geom_col(width = 0.7) +
  facet_wrap(~Metric, nrow = 1) +
  theme_bw() +
  ylim(0, 1) +
  labs(
    title = "Benchmarking doublet detection performance",
    x = NULL,
    y = "Score"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  "figures/Figure2_benchmark_metrics_3methods.png",
  p_metrics,
  width = 11,
  height = 4,
  dpi = 300
)

############################
# 3. If scDblFinder Seurat object exists, draw 4-panel UMAP
############################

if (file.exists("results/seurat_object_with_scDblFinder.rds")) {
  
  obj <- readRDS("results/seurat_object_with_scDblFinder.rds")
  
  # 如果之前没有把 Scrublet prediction 写进这个 object，再从 csv 加进去
  if (!"Scrublet_pred" %in% colnames(obj@meta.data)) {
    scrublet_pred <- read_csv("results/scrublet_predictions.csv", show_col_types = FALSE)
    scrublet_pred <- as.data.frame(scrublet_pred)
    rownames(scrublet_pred) <- scrublet_pred$barcode
    scrublet_pred <- scrublet_pred[colnames(obj), ]
    
    obj$Scrublet_pred <- scrublet_pred$Scrublet_pred
    obj$Scrublet_score <- scrublet_pred$Scrublet_score
  }
  
  p1 <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "ground_truth",
    pt.size = 0.2
  ) + ggtitle("Ground truth")
  
  p2 <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "DoubletFinder_pred",
    pt.size = 0.2
  ) + ggtitle("DoubletFinder")
  
  p3 <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "Scrublet_pred",
    pt.size = 0.2
  ) + ggtitle("Scrublet")
  
  p4 <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "scDblFinder_pred",
    pt.size = 0.2
  ) + ggtitle("scDblFinder")
  
  p_umap <- (p1 | p2) / (p3 | p4)
  
  ggsave(
    "figures/Figure1_UMAP_4methods.png",
    p_umap,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  # 附加：三方法 score distribution
  score_df <- data.frame(
    barcode = colnames(obj),
    ground_truth = obj$ground_truth,
    DoubletFinder = obj$DoubletFinder_score,
    Scrublet = obj$Scrublet_score,
    scDblFinder = obj$scDblFinder_score
  )
  
  score_long <- score_df %>%
    pivot_longer(
      cols = c(DoubletFinder, Scrublet, scDblFinder),
      names_to = "Method",
      values_to = "Score"
    )
  
  p_scores <- ggplot(score_long, aes(x = ground_truth, y = Score)) +
    geom_boxplot(outlier.size = 0.2) +
    facet_wrap(~Method, scales = "free_y") +
    theme_bw() +
    labs(
      title = "Doublet score distributions by ground-truth label",
      x = "Ground truth",
      y = "Doublet score"
    )
  
  ggsave(
    "figures/Supp_score_distributions_3methods.png",
    p_scores,
    width = 10,
    height = 4,
    dpi = 300
  )
  
  saveRDS(obj, "results/final_seurat_object_3methods.rds")
  
  cat("\n4-panel UMAP generated successfully.\n")
  
} else {
  cat("\nNo seurat_object_with_scDblFinder.rds found.\n")
  cat("Metric plot was generated, but 4-panel UMAP could not be created.\n")
}
