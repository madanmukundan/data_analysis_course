# Import necessary libraries
import argparse
from pathlib import Path

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Parse arguments
parser = argparse.ArgumentParser(description='Filter cells based on QC metrics')
parser.add_argument("--input_file", type=Path, required=True, help="Input file")
parser.add_argument("--output_file", type=Path, required=True, help="Output file with good-quality cells")
parser.add_argument("--qc_violin_plot", type=Path, help="Path to save violin plot of QC metrics")
# TODO: Add the scatterplot path from snakemake rule so python can use it while executing the script
# TODO: Add the parameter from snakemake so python can use it while executing the script
parser.add_argument("--scatter_plot", type=Path, help="Scatter plot")
parser.add_argument("--min_num_genes", type=int, help="Minimum numer of genes")
args = parser.parse_args()

# Create directory for figures if one does not exist
# NOTE: Unless you save the scatterplot elsewhere, you do not need to repeat for scatterplot
args.qc_violin_plot.parent.mkdir(parents=True, exist_ok=True)

# Load input data
adata = sc.read_h5ad(args.input_file)

# Add sample information from the file path
adata.obs["sample"] = args.input_file.stem

# # # Calculating QC metrics and thresholds
# Find mitochondrial genes
var_names = adata.var_names.astype("str")
mito_genes = var_names.str.startswith("mt-")
# Annotate mitochondrial genes in adata.var
adata.var["mito"] = mito_genes

adata.obs["fraction_mito"] = (
        1 + np.sum(adata[:, mito_genes].X, axis=1).A1) / (1 + np.sum(adata.X, axis=1).A1)
adata.obs["n_counts"] = adata.X.sum(axis=1).A1
adata.obs["n_genes"] = (adata.X > 0).sum(axis=1).A1

# NOTE: Edit the min_n_genes so that it uses your snakemake parameter instead of a hardcoded value
min_n_genes = args.min_num_genes
min_n_counts = 1000
# NOTE: Edit the fraction_mito threshold so it filters everything above mean + 2*sd ~0.079
# max_fraction_mito = np.mean(adata.obs["fraction_mito"]) + 2 * np.std(adata.obs["fraction_mito"])
max_fraction_mito = adata.obs["fraction_mito"].mean() + 2 * adata.obs["fraction_mito"].std()
# max_fraction_mito = 0.05

# # # Filtering

# Label cells in unfiltered data based on whether they pass each threshold
adata.obs["pass_n_genes"] = adata.obs["n_genes"] >= min_n_genes
adata.obs["pass_n_counts"] = adata.obs["n_counts"] >= min_n_counts
adata.obs["pass_fraction_mito"] = adata.obs["fraction_mito"] <= max_fraction_mito
adata.obs["pass_all_qc"] = adata.obs["pass_n_genes"] & adata.obs["pass_n_counts"] & adata.obs["pass_fraction_mito"]

print(f"Filtering cells based on QC metrics: n_genes >= {min_n_genes}, n_counts >= {min_n_counts}, fraction_mito <= {max_fraction_mito}")
print(f"Number of cells before filtering: {adata.n_obs}")
print(f"Number of cells passing all QC thresholds: {np.sum(adata.obs['pass_all_qc'])}")
print(f"Fraction of cells passing all QC thresholds: {np.mean(adata.obs['pass_all_qc'])}")

# Save filtered data
adata_filtered = adata[adata.obs["pass_all_qc"]].copy()
adata_filtered.write_h5ad(args.output_file)

# # # Plotting
# Plot QC violin plots with seaborn
# Set the max y-value for each axis manually to allow comparisons between samples
# Note for students: The max value was decided post-hoc based on the data in this case
fig, axes = plt.subplots(1, 3, figsize=(10, 5))
sns.set_theme(style="whitegrid")
sns.violinplot(y=adata.obs["n_genes"], ax=axes[0])
sns.violinplot(y=adata.obs["n_counts"], ax=axes[1])
sns.violinplot(y=adata.obs["fraction_mito"], ax=axes[2])
axes[0].set_ylim(-100, 8000)
axes[1].set_ylim(-1000, 48000)
axes[2].set_ylim(-0.05, 0.55)
fig.suptitle(f"QC metrics for {adata.obs['sample'].unique()[0]}")
fig.tight_layout()
plt.savefig(args.qc_violin_plot)

# NOTE: Add scatterplot with seaborn and save to location specified in snakemake rule
fig, axes = plt.subplots(1, 3, figsize=(10, 5))
# plt.style.use('ggplot')
sns.set_theme(style="whitegrid")
sns.scatterplot(x=adata.obs["n_genes"],y=adata.obs["n_counts"], ax=axes[0])
# axes[0].scatter(x=adata.obs["n_genes"],y=adata.obs["n_counts"])
axes[0].set_xlim(-100, 8000) # genes
axes[0].set_ylim(-1000, 48000) # counts
axes[0].hlines(min_n_counts, xmin=0, xmax=7000, color='r')
axes[0].vlines(min_n_genes, ymin=-1000, ymax=40000,  color='r')

sns.scatterplot(x=adata.obs["n_genes"],y=adata.obs["fraction_mito"], ax=axes[1])
# axes[1].scatter(x=adata.obs["n_genes"],y=adata.obs["fraction_mito"])
axes[1].set_xlim(-100, 8000) # genes
axes[1].set_ylim(-0.05, 0.55) # frac mito
axes[1].hlines(max_fraction_mito, xmin=0, xmax=7000, color='r')
axes[1].vlines(min_n_genes, ymin=0, ymax=0.3, color='r')

sns.scatterplot(x=adata.obs["n_counts"], y=adata.obs["fraction_mito"], ax=axes[2])
# axes[2].scatter(x=adata.obs["n_counts"], y=adata.obs["fraction_mito"])
axes[2].set_ylim(-1000, 48000) # counts
axes[2].set_ylim(-0.05, 0.55) # mito
axes[2].hlines(max_fraction_mito,xmin=-1000, xmax=48000,color='r')
axes[2].vlines(min_n_counts, ymin=0, ymax=0.3,color='r')
# NOTE: Add the thresholds as horizontal and vertical lines to the plot
fig.suptitle(f"QC metrics for {adata.obs['sample'].unique()[0]}")
fig.tight_layout()
plt.savefig(args.scatter_plot)
