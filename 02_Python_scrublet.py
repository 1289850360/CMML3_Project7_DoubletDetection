import os
import time
import numpy as np
import pandas as pd
import scipy.io
import scrublet as scr
from sklearn.metrics import (
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    average_precision_score,
    confusion_matrix
)

data_dir = r"C:/Users/pjk/Documents/CMML3_Project7/GSE96583"
os.chdir(data_dir)

os.makedirs("results", exist_ok=True)

counts = scipy.io.mmread("results/counts_for_scrublet.mtx").tocsr()

# R 导出的是 genes x cells，Scrublet 需要 cells x genes
counts = counts.T.tocsr()

meta = pd.read_csv("results/metadata_for_scrublet.csv")

truth = (meta["ground_truth"] == "doublet").astype(int).values
expected_rate = truth.mean()

print("Counts shape, cells x genes:", counts.shape)
print("Ground truth doublet rate:", expected_rate)

start = time.time()

scrub = scr.Scrublet(
    counts,
    expected_doublet_rate=expected_rate,
    random_state=123
)

doublet_scores, predicted_doublets = scrub.scrub_doublets(
    min_counts=2,
    min_cells=3,
    min_gene_variability_pctl=85,
    n_prin_comps=20
)

runtime = time.time() - start

pred = predicted_doublets.astype(int)

precision = precision_score(truth, pred, zero_division=0)
recall = recall_score(truth, pred, zero_division=0)
f1 = f1_score(truth, pred, zero_division=0)

try:
    auroc = roc_auc_score(truth, doublet_scores)
except Exception:
    auroc = np.nan

try:
    auprc = average_precision_score(truth, doublet_scores)
except Exception:
    auprc = np.nan

tn, fp, fn, tp = confusion_matrix(truth, pred).ravel()

metrics = pd.DataFrame([{
    "Method": "Scrublet",
    "TP": tp,
    "FP": fp,
    "FN": fn,
    "TN": tn,
    "Precision": precision,
    "Recall": recall,
    "F1": f1,
    "AUROC": auroc,
    "AUPRC": auprc,
    "Predicted_doublet_rate": pred.mean(),
    "Runtime_seconds": runtime
}])

scrublet_result = pd.DataFrame({
    "barcode": meta["barcode"],
    "ground_truth": meta["ground_truth"],
    "Scrublet_score": doublet_scores,
    "Scrublet_pred": np.where(predicted_doublets, "doublet", "singlet")
})

scrublet_result.to_csv("results/scrublet_predictions.csv", index=False)
metrics.to_csv("results/metrics_scrublet.csv", index=False)

print(metrics)
print("Finished Scrublet successfully.")