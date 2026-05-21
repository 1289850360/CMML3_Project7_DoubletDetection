############################################################
# Final figures and benchmark table
# Combine DoubletFinder + Scrublet results
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
# 1. Load files
############################

obj <- readRDS("results/seurat_object_DoubletFinder_only.rds")

metrics_df <- read_csv("results/metrics_DoubletFinder.csv", show_col_types = FALSE)
metrics_scrublet <- read_csv("results/metrics_scrublet.csv", show_col_types = FALSE)
scrublet_pred <- read_csv("results/scrublet_predictions.csv", show_col_types = FALSE)

############################
# 2. Add Scrublet prediction to Seurat object
############################

scrublet_pred <- as.data.frame(scrublet_pred)
rownames(scrublet_pred) <- scrublet_pred$barcode

# Make sure order matches Seurat object
scrublet_pred <- scrublet_pred[colnames(obj), ]

obj$Scrublet_pred <- scrublet_pred$Scrublet_pred
obj$Scrublet_score <- scrublet_pred$Scrublet_score

############################
# 3. Combine metrics
############################

metrics_all <- bind_rows(metrics_df, metrics_scrublet) %>%
  select(
    Method, TP, FP, FN, TN,
    Precision, Recall, F1, AUROC, AUPRC,
    Predicted_doublet_rate,
    everything()
  )

write.csv(metrics_all, "results/final_benchmark_metrics.csv", row.names = FALSE)

print(metrics_all)

############################
# 4. Main Figure 1: UMAP comparison
############################

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

p_umap <- p1 | p2 | p3

ggsave(
  "figures/Figure1_UMAP_groundtruth_DoubletFinder_Scrublet.png",
  p_umap,
  width = 13,
  height = 4.5,
  dpi = 300
)

############################
# 5. Main Figure 2: metric bar plot
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
  "figures/Figure2_benchmark_metrics.png",
  p_metrics,
  width = 10,
  height = 4,
  dpi = 300
)

############################
# 6. Supplementary: score distributions
############################

score_df <- data.frame(
  barcode = colnames(obj),
  ground_truth = obj$ground_truth,
  DoubletFinder = obj$DoubletFinder_score,
  Scrublet = obj$Scrublet_score
)

score_long <- score_df %>%
  pivot_longer(
    cols = c(DoubletFinder, Scrublet),
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
  "figures/Supp_score_distributions.png",
  p_scores,
  width = 8,
  height = 4,
  dpi = 300
)

############################
# 7. Supplementary: confusion matrix table
############################

confusion_df <- metrics_all %>%
  select(Method, TP, FP, FN, TN)

write.csv(confusion_df, "results/confusion_matrix_summary.csv", row.names = FALSE)

saveRDS(obj, "results/final_seurat_object_with_predictions.rds")

cat("\nFinished final figure script.\n")
cat("Main outputs:\n")
cat("results/final_benchmark_metrics.csv\n")
cat("figures/Figure1_UMAP_groundtruth_DoubletFinder_Scrublet.png\n")
cat("figures/Figure2_benchmark_metrics.png\n")
cat("figures/Supp_score_distributions.png\n")