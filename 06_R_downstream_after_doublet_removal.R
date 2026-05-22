# ============================================================
# 06_R_downstream_after_doublet_removal.R
# Downstream check: reclustering after removing predicted doublets
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# -----------------------------
# 1. Set project paths safely
# -----------------------------

# 手动选择你的 GSE96583 文件夹
# 请选择这个文件夹：文档/CMML3_Project7/GSE96583
project_dir <- choose.dir(caption = "Please select your GSE96583 project folder")

results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")

cat("Project folder is:\n", project_dir, "\n")
cat("Results folder is:\n", results_dir, "\n")
cat("Figures folder is:\n", figures_dir, "\n")

# 检查 results 文件夹是否存在
if (!dir.exists(results_dir)) {
  stop("The results folder was not found. Please make sure you selected the GSE96583 folder, not the results folder itself.")
}

# 检查 figures 文件夹是否存在，不存在才创建
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

# 看一下 results 里有哪些文件
cat("\nFiles in results folder:\n")
print(list.files(results_dir))

# -----------------------------
# 2. Load Seurat object
# -----------------------------

rds_path <- file.path(results_dir, "final_seurat_object_3methods.rds")

if (!file.exists(rds_path)) {
  rds_path <- file.path(results_dir, "final_seurat_object_with_predictions.rds")
}

if (!file.exists(rds_path)) {
  stop("Cannot find final_seurat_object_3methods.rds or final_seurat_object_with_predictions.rds in results folder.")
}

seu <- readRDS(rds_path)

cat("\nLoaded object:\n", rds_path, "\n")
cat("Number of cells:", ncol(seu), "\n")
cat("\nMetadata columns:\n")
print(colnames(seu@meta.data))
# -----------------------------
# 3. Automatically find prediction columns
# -----------------------------
meta_cols <- colnames(seu@meta.data)

find_col <- function(patterns, meta_cols) {
  hits <- c()
  for (p in patterns) {
    hits <- c(hits, grep(p, meta_cols, ignore.case = TRUE, value = TRUE))
  }
  hits <- unique(hits)
  if (length(hits) == 0) return(NA)
  return(hits[1])
}

df_col <- find_col(c("DoubletFinder", "DF.class", "pANN"), meta_cols)
scrub_col <- find_col(c("Scrublet", "scrublet"), meta_cols)
scdbl_col <- find_col(c("scDblFinder", "scDblFinder.class", "scDblFinder_class"), meta_cols)

cat("\nDetected prediction columns:\n")
cat("DoubletFinder:", df_col, "\n")
cat("Scrublet:", scrub_col, "\n")
cat("scDblFinder:", scdbl_col, "\n")

if (is.na(df_col) | is.na(scrub_col) | is.na(scdbl_col)) {
  stop("Could not automatically find all three prediction columns. Check colnames(seu@meta.data) printed above.")
}

# -----------------------------
# 4. Helper function: convert labels to singlet/doublet
# -----------------------------
standardise_label <- function(x) {
  x <- as.character(x)
  x_low <- tolower(x)
  
  out <- rep(NA, length(x_low))
  
  out[grepl("singlet|single", x_low)] <- "singlet"
  out[grepl("doublet|double", x_low)] <- "doublet"
  
  # Some tools may use TRUE/FALSE or 1/0
  out[x_low %in% c("true", "1")] <- "doublet"
  out[x_low %in% c("false", "0")] <- "singlet"
  
  return(out)
}

seu$DF_label_std <- standardise_label(seu@meta.data[[df_col]])
seu$Scrublet_label_std <- standardise_label(seu@meta.data[[scrub_col]])
seu$scDblFinder_label_std <- standardise_label(seu@meta.data[[scdbl_col]])

cat("\nLabel summaries:\n")
cat("\nDoubletFinder:\n")
print(table(seu$DF_label_std, useNA = "ifany"))
cat("\nScrublet:\n")
print(table(seu$Scrublet_label_std, useNA = "ifany"))
cat("\nscDblFinder:\n")
print(table(seu$scDblFinder_label_std, useNA = "ifany"))

# -----------------------------
# 5. Function: remove doublets and rerun clustering/UMAP
# -----------------------------
rerun_after_filtering <- function(obj, label_col = NULL, method_name = "Original",
                                  dims_use = 1:20, resolution_use = 0.5) {
  
  if (!is.null(label_col)) {
    keep_cells <- rownames(obj@meta.data)[obj@meta.data[[label_col]] == "singlet"]
    obj <- subset(obj, cells = keep_cells)
  }
  
  # rerun standard Seurat workflow
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution_use, verbose = FALSE)
  obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)
  
  n_cells <- ncol(obj)
  n_clusters <- length(unique(obj$seurat_clusters))
  
  p <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE,
    pt.size = 0.25
  ) +
    ggtitle(paste0(method_name, "\n", n_cells, " cells; ", n_clusters, " clusters")) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.position = "none"
    )
  
  return(list(object = obj, plot = p, cells = n_cells, clusters = n_clusters))
}

# -----------------------------
# 6. Run downstream comparison
# -----------------------------
original_res <- rerun_after_filtering(
  seu,
  label_col = NULL,
  method_name = "Original"
)

df_res <- rerun_after_filtering(
  seu,
  label_col = "DF_label_std",
  method_name = "After DoubletFinder removal"
)

scrub_res <- rerun_after_filtering(
  seu,
  label_col = "Scrublet_label_std",
  method_name = "After Scrublet removal"
)

scdbl_res <- rerun_after_filtering(
  seu,
  label_col = "scDblFinder_label_std",
  method_name = "After scDblFinder removal"
)

# -----------------------------
# 7. Save UMAP comparison figure
# -----------------------------
fig_s3 <- (original_res$plot | df_res$plot) /
  (scrub_res$plot | scdbl_res$plot)

fig_s3 <- fig_s3 +
  plot_annotation(
    title = "Downstream clustering after removing predicted doublets",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

ggsave(
  filename = file.path(figures_dir, "FigureS3_downstream_umap_after_doublet_removal.png"),
  plot = fig_s3,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  filename = file.path(figures_dir, "FigureS3_downstream_umap_after_doublet_removal.pdf"),
  plot = fig_s3,
  width = 10,
  height = 8
)

cat("\nSaved Figure S3 to figures folder.\n")

# -----------------------------
# 8. Save summary table
# -----------------------------
summary_table <- data.frame(
  Dataset = c(
    "Original",
    "After DoubletFinder removal",
    "After Scrublet removal",
    "After scDblFinder removal"
  ),
  Cells_retained = c(
    original_res$cells,
    df_res$cells,
    scrub_res$cells,
    scdbl_res$cells
  ),
  Cells_removed = c(
    0,
    original_res$cells - df_res$cells,
    original_res$cells - scrub_res$cells,
    original_res$cells - scdbl_res$cells
  ),
  Percent_removed = round(c(
    0,
    (original_res$cells - df_res$cells) / original_res$cells * 100,
    (original_res$cells - scrub_res$cells) / original_res$cells * 100,
    (original_res$cells - scdbl_res$cells) / original_res$cells * 100
  ), 2),
  Number_of_clusters = c(
    original_res$clusters,
    df_res$clusters,
    scrub_res$clusters,
    scdbl_res$clusters
  )
)

write.csv(
  summary_table,
  file = file.path(results_dir, "downstream_doublet_removal_summary.csv"),
  row.names = FALSE
)

cat("\nDownstream summary table:\n")
print(summary_table)

# -----------------------------
# 9. Optional marker gene DotPlot
# -----------------------------
markers <- c(
  "CD3D", "CD3E",
  "MS4A1", "CD79A",
  "NKG7", "GNLY",
  "LST1", "S100A8",
  "PPBP", "PF4"
)

markers_present <- markers[markers %in% rownames(scdbl_res$object)]

cat("\nMarkers present in dataset:\n")
print(markers_present)

if (length(markers_present) >= 3) {
  p_marker <- DotPlot(
    scdbl_res$object,
    features = markers_present,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis() +
    ggtitle("Marker gene expression after scDblFinder doublet removal") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(
    filename = file.path(figures_dir, "FigureS4_marker_dotplot_after_scDblFinder_removal.png"),
    plot = p_marker,
    width = 9,
    height = 5,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(figures_dir, "FigureS4_marker_dotplot_after_scDblFinder_removal.pdf"),
    plot = p_marker,
    width = 9,
    height = 5
  )
  
  cat("\nSaved optional marker DotPlot as Figure S4.\n")
} else {
  cat("\nNot enough marker genes found. Marker DotPlot skipped.\n")
}

# -----------------------------
# 10. Save filtered Seurat objects for record
# -----------------------------
saveRDS(df_res$object, file.path(results_dir, "seurat_after_DoubletFinder_removal.rds"))
saveRDS(scrub_res$object, file.path(results_dir, "seurat_after_Scrublet_removal.rds"))
saveRDS(scdbl_res$object, file.path(results_dir, "seurat_after_scDblFinder_removal.rds"))

cat("\nAll done.\n")