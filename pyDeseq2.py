#!/usr/bin/env python3
# coding: utf-8

"""
Sample script for analyzing RNA-seq data using PyDESeq2, Scanpy, and other libraries.
Code skeleton from: https://github.com/linkangit/RNAseq-analysis-PyDESeq2/tree/main?tab=readme-ov-file

Usage:
    python alt_rna_analysis.py
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
from pathlib import Path
import random  # used later for picking random sets of genes
import sys

# PyDESeq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# GSEA
import gseapy as gp
from gseapy.plot import gseaplot

# For label adjustment
from adjustText import adjust_text

#######################################################
# Misc Notes
#######################################################
# P. putida KT2440
# gene_annotation_csv_path = 'C:/Users/gebeb/Desktop/GitHub/GESA/data_P_putida/gene_annotations.csv'
# gene_annotation_gname_col = 'PGD Gene ID'  # Column that matches gene names in 'read_counts.csv'
# gene_annotation_gene_symbol_col = 'Locus Tag'  # For P. putida - PP_4718, PP_0025, etc.

# P. aeruginosa PAO1
# sample_info_file_path = 'C:/Users/gebeb/Desktop/Research/RNAseq/PAO1_PRJNA739737/sample_info.csv'
# read_counts_csv_path = 'C:/Users/gebeb/Desktop/Research/RNAseq/PAO1_PRJNA739737/read_counts.csv'
# gene_annotation_csv_path = 'C:/Users/gebeb/Desktop/GitHub/GESA/data_P_aeruginosa/gene_annotations.csv'
# gene_annotation_gname_col = 'PGD Gene ID'
# gene_annotation_gene_symbol_col = 'Locus Tag'

# E. coli K12
# ...


#######################################################
# 1. Prepare sample information (metadata)
#######################################################
sample_info_file_path = 'C:/Users/gebeb/Desktop/Research/RNAseq/PAO1_PRJNA892976/sample_info.csv'
read_counts_csv_path = 'C:/Users/gebeb/Desktop/Research/RNAseq/PAO1_PRJNA892976/read_counts.csv'
gene_annotation_csv_path = 'C:/Users/gebeb/Desktop/GitHub/GESA/data_P_aeruginosa/gene_annotations.csv'
gene_annotation_gname_col = 'PGD Gene ID'  # Column that matches gene names in 'read_counts.csv'
gene_annotation_gene_symbol_col = 'Locus Tag'  # For P. putida - PP_4718, PP_0025, etc.
control_condition = 'Native'
variant_condition = 'Geraniol'
n_labels_volcano = 50 # Number of genes to label on the volcano plot

# Create folder to save plots
results_folder = "DGE_results"
Path(results_folder).mkdir(parents=True, exist_ok=True)

# sample_info_dict = {
#     "Sample": ["Kelley_17", "Kelley_18", "Kelley_19", "Kelley_20",
#                "Kelley_21", "Kelley_22", "Kelley_23", "Kelley_24"],
#     "Condition": ["C", "C", "C", "C",
#                   "mut", "mut", "mut", "mut"]
# }


sample_info = pd.read_csv(sample_info_file_path)
print("SAMPLE INFO:\n", sample_info, "\n")

# Identify control vs. mutant groups
ctrl_list = sample_info.loc[sample_info['Condition'] == control_condition, 'Sample'].tolist()
variant_list = sample_info.loc[sample_info['Condition'] == variant_condition, 'Sample'].tolist()
print("Control IDs:", ctrl_list)
print("Mutant IDs:", variant_list, "\n")

#######################################################
# 2. Load and clean up read-count data
#######################################################

read_df = pd.read_csv(read_counts_csv_path)
print("Initial shape of read-count table:", read_df.shape)
print("Column headers in read-count table:\n", list(read_df.columns))

# Make 'gname' column the index
read_df = read_df.set_index('gname')

# Exclude genes with zero total counts
read_df = read_df[read_df.sum(axis=1) > 0]

#######################################################
# 3. Filter columns based on metadata info
#######################################################

# Keep copy of full data tables
full_read_df = read_df.copy()
full_sample_info = sample_info.copy()
full_sample_info = full_sample_info.set_index("Sample")
# Match index order of sample info to read counts table
full_sample_info = full_sample_info.loc[full_read_df.T.index]

# Keep only samples that are in filtered read_df
sample_info = sample_info.loc[(sample_info["Condition"] == control_condition) | (sample_info["Condition"] == variant_condition)].copy()

select_samples = ctrl_list + variant_list
read_df = read_df[select_samples]
print("\nFiltered read-count table dimensions:", read_df.shape)
print("Filtered columns:", list(read_df.columns))

# Adjust sample_info indexing for DESeq2
sample_info.set_index('Sample', inplace=True)
sample_info["Condition"] = sample_info["Condition"].astype("category")
sample_info["Condition"] = sample_info["Condition"].cat.reorder_categories([control_condition, variant_condition], ordered=True)

#######################################################
# 4. Create DESeqDataSet and run differential analysis
#######################################################

# PyDESeq2 expects samples in rows and genes in columns
transposed_df = read_df.T

# Ensure indices match exactly
sample_info = sample_info.loc[transposed_df.index]

my_dds = DeseqDataSet(
    counts=transposed_df,
    metadata=sample_info,
    design_factors=["Condition"]
)

print("\nPyDESeq2 dataset before DESeq2 run:\n", my_dds)

my_dds.deseq2()
print("\nPyDESeq2 dataset after DESeq2 run:\n", my_dds)

# Create and save Variance Stabilizing Transformation (VST) with full data
full_dds = DeseqDataSet(
    counts=full_read_df.T,
    metadata=full_sample_info,
    design_factors=["Condition"]
)

full_dds.vst_fit(use_design=False)
vst_counts = full_dds.vst_transform()
vst_df = pd.DataFrame(vst_counts.T, index=full_read_df.index, columns=full_read_df.columns)
vst_df.to_csv(f"{results_folder}/VST_normalized_counts.csv", index=True)
print("\nSaved VST-normalized counts table.")

#######################################################
# 5. Extract DE outputs
#######################################################

my_stats = DeseqStats(my_dds, contrast=["Condition", variant_condition, control_condition])
stats_summary = my_stats.summary()
results_df = my_stats.results_df

print("\nHead of DE outputs:\n", results_df.head())

#######################################################
# 6. Merge annotation details
#######################################################

gene_info = pd.read_csv(gene_annotation_csv_path)
gene_info.set_index(gene_annotation_gname_col, inplace=True)
gene_info = gene_info.drop_duplicates(subset=gene_annotation_gene_symbol_col, keep="last") # Remove rows with duplicate gene ID labels

# Example: filter out rows without symbols
gene_info['check_symbol'] = gene_info[gene_annotation_gene_symbol_col].index.isnull()
gene_info = gene_info.loc[gene_info['check_symbol'] == False]

complete_df = pd.merge(
    left=results_df,
    right=gene_info,
    how="inner",
    left_index=True,
    right_index=True
)

# Incorporate Symbol into results_df
results_df[gene_annotation_gene_symbol_col] = complete_df[gene_annotation_gene_symbol_col]

# Filter out low-expression genes
results_df = results_df[results_df.baseMean >= 10]

# Find significant genes
significant_hits = results_df[(results_df.padj < 0.1) & (abs(results_df.log2FoldChange) > 0.5)]
print("\nSignificant Genes:\n", significant_hits)

#######################################################
# 7. Exploratory data analysis with Scanpy
#######################################################

try:
    norm_data = my_dds.layers["normed_counts"]
except AttributeError:
    print("\n[Note] 'my_dds.layers[\"normed_counts\"]' not found. Check PyDESeq2 docs.")
    norm_data = None

if norm_data is not None:
    # Build AnnData object
    adata = anndata.AnnData(X=norm_data, obs=my_dds.obs, var=pd.DataFrame(index=my_dds.var.index))
    adata.layers["counts"] = norm_data
    adata.layers["log1p"] = np.log1p(norm_data)

    # PCA Plot TODO: Edit plot to include ALL sample types
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata)
    sc.pl.pca(adata, color='Condition', size=150, show=False)
    plt.title("PCA of Samples")
    plt.savefig(f"{results_folder}/PCA_{control_condition}_vs_{variant_condition}.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Distance Heatmap
    corr_mat = np.corrcoef(norm_data)
    sns.clustermap(
        corr_mat,
        row_cluster=True, col_cluster=True,
        cmap="vlag",
        xticklabels=adata.obs.index, yticklabels=adata.obs.index
    )
    plt.title("Sample Correlation Heatmap")
    plt.savefig(f"{results_folder}/Correlation_Heatmap_{control_condition}_vs_{variant_condition}.png", dpi=300, bbox_inches="tight")
    plt.close()

    # MA Plot
    plt.figure(figsize=(6, 6))
    plt.scatter(
        x=results_df['baseMean'],
        y=results_df['log2FoldChange'],
        c=results_df['log2FoldChange'],
        alpha=0.5,
        cmap='viridis'
    )
    plt.xscale('log')
    plt.axhline(y=0, color='red', linestyle='--')
    plt.xlabel("Mean Normalized Counts (log scale)")
    plt.ylabel("Log2 FC")
    plt.title("MA Plot")
    plt.savefig(f"{results_folder}/MA_plot_{control_condition}_vs_{variant_condition}.png", dpi=300, bbox_inches="tight")
    plt.close()

#######################################################
# Volcano Plot Section
#######################################################

results_df = results_df.copy()
results_df.rename(columns={gene_annotation_gene_symbol_col: "symbol"}, inplace=True)
results_df = results_df.dropna(subset=["padj", "symbol", "log2FoldChange", "baseMean"])

# Replace zeros with a small epsilon (for edge cases where log10(0) occurs)
epsilon = 1e-300
results_df["padj"] = results_df["padj"].clip(lower=epsilon)

results_df["nlog10"] = -np.log10(results_df["padj"])

# Random picks for color categories
# picked_set1 = random.choices(results_df.symbol.tolist(), weights=results_df.nlog10.tolist(), k=250)
# picked_set2 = random.choices(results_df.symbol.tolist(), weights=results_df.nlog10.tolist(), k=300)
# picked_set2 = [x for x in picked_set2 if x not in picked_set1]

def assign_color(row):
    fold_change, gene_symbol, minuslog10 = row
    if abs(fold_change) < 1 or minuslog10 < 2:
        return "not_significant"
    if fold_change > 1:
        return "upregulated"
    if fold_change < -1:
        return "downregulated"

results_df["Color"] = results_df[["log2FoldChange", "symbol", "nlog10"]].apply(assign_color, axis=1)

# picked_set3 = random.choices(results_df.symbol.tolist(), weights=results_df.nlog10.tolist(), k=250)
# picked_set4 = random.choices(results_df.symbol.tolist(), weights=results_df.nlog10.tolist(), k=300)
# picked_set4 = [x for x in picked_set4 if x not in picked_set3]

# def assign_shape(gsym):
#     if gsym in picked_set3:
#         return "shapeA"
#     if gsym in picked_set4:
#         return "shapeB"
#     return "shapeNeutral"

# results_df["shape"] = results_df.symbol.map(assign_shape)

# Filter out the super highly expressed genes
results_df = results_df.sort_values(by=["nlog10"], ascending=False).iloc[5:,:]

# Plotting
plt.figure(figsize=(6, 6))
ax = sns.scatterplot(
    data=results_df,
    x="log2FoldChange",
    y="nlog10",
    hue="Color",
    hue_order=["not_significant", "upregulated", "downregulated"],
    palette=["lightgrey", "orange", "purple"],
    # style="shape",
    # style_order=["shapeA", "shapeB", "shapeNeutral"],
    # markers=["^", "s", "o"],
    size="baseMean",
    sizes=(40, 400)
)

ax.axhline(2, zorder=0, c="k", lw=2, ls="--")
ax.axvline(1, zorder=0, c="k", lw=2, ls="--")
ax.axvline(-1, zorder=0, c="k", lw=2, ls="--")

### Label highly significant points ###
# label_texts = []
# for idx in range(len(results_df)):
#     myrow = results_df.iloc[idx]
#     if myrow.nlog10 > 5 and abs(myrow.log2FoldChange) > 2:
#         label_texts.append(
#             plt.text(
#                 x=myrow.log2FoldChange,
#                 y=myrow.nlog10,
#                 s=myrow.symbol,
#                 fontsize=12, weight="bold"
#             )
#         )

### Label genes with -1 > log2FoldChange > 1 ###
subset = results_df[(abs(results_df.log2FoldChange) > 1) & (results_df.nlog10 > 2)]
# keep top N by magnitude
subset = subset.nlargest(n_labels_volcano, "nlog10")

label_texts = [
    plt.text(row.log2FoldChange, row.nlog10, row.symbol, fontsize=10, weight="bold")
    for _, row in subset.iterrows()
]

adjust_text(label_texts, arrowprops=dict(arrowstyle="-", color="k"))

plt.legend(loc=1, bbox_to_anchor=(1.4, 1), frameon=False, prop={"weight": "bold"})

for border in ["bottom", "left"]:
    ax.spines[border].set_linewidth(2)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(width=2)

plt.xticks(size=12, weight="bold")
plt.yticks(size=12, weight="bold")
plt.xlabel("$\\log_{2}$FoldChange", size=15)
plt.ylabel("$-\\log_{10}$FDR", size=15)

plt.savefig(f"{results_folder}/Volcano_n{n_labels_volcano}_{control_condition}_vs_{variant_condition}.png", dpi=300, bbox_inches="tight", facecolor="white")
# plt.show()
print(f"\nVolcano plot saved as 'Volcano_{control_condition}_vs_{variant_condition}.png'.")

##########################################################
# Final steps: Summaries, heatmaps, and saving results
##########################################################

# A) Identify significant hits again
significant_hits = results_df[(results_df["padj"] < 0.1) & (abs(results_df["log2FoldChange"]) > 0.5)].copy()
if len(significant_hits) == 0:
    print("No genes match significance criteria.")
else:
    print(f"Total significant genes found: {len(significant_hits)}")

# B) Take top 10 up and top 10 down from significant group
top_up10 = significant_hits.sort_values("log2FoldChange", ascending=False).head(10)
top_down10 = significant_hits.sort_values("log2FoldChange", ascending=True).head(10)
highlighted_genes = list(top_up10.index) + list(top_down10.index)

# C) Extract normalized counts for heatmap
try:
    norm_data = my_dds.layers["normed_counts"]
except AttributeError:
    raise ValueError("Could not locate 'my_dds.layers[\"normed_counts\"]'. "
                     "Verify your PyDESeq2 version includes normalized counts.")

norm_data_df = pd.DataFrame(
    data=norm_data,
    index=my_dds.obs_names,
    columns=my_dds.var.index
)

########################################################## TODO: Replace gene IDs with Locus Tags
# HEATMAP A: All significantly altered genes
##########################################################
if len(significant_hits) > 0:
    sign_data_df = norm_data_df[significant_hits.index]
    sns.clustermap(sign_data_df.T, z_score=0, cmap="viridis", figsize=(8, 10))
    plt.title("Cluster Heatmap - Significant Genes")
    plt.savefig(f"{results_folder}/Heatmap_All_Sig_Genes_{control_condition}_vs_{variant_condition}.png", dpi=300, bbox_inches="tight")
    # plt.show()

########################################################## TODO: Replace gene IDs with Locus Tags
# HEATMAP B: Top 20 (10 up, 10 down)
##########################################################
if len(highlighted_genes) > 0:
    top_data_df = norm_data_df[highlighted_genes]
    sns.clustermap(top_data_df.T, z_score=0, cmap="magma", figsize=(8, 10))
    plt.title("Cluster Heatmap - Top 10 Up & Top 10 Down")
    plt.savefig(f"{results_folder}/Heatmap_Top20_Up_Down_{control_condition}_vs_{variant_condition}.png", dpi=300, bbox_inches="tight")
    # plt.show()

##########################################################
# Save the final DE results
##########################################################
results_df.to_csv(f"{results_folder}/DE_results_{control_condition}_vs_{variant_condition}.csv")
print(f"DE results saved to 'DE_results_{control_condition}_vs_{variant_condition}.csv'.")

##########################################################
# Gene regulatory networks from significantly DE genes
##########################################################
# TODO: REFERENCE CODE FROM 'GESA/networks.py'
# 1. Co-expression heatmap - ALL GENES (see Junier, I. et al. 2016. PLOS One.)
# 2. Co-expression network - DEGs ONLY (see 'GESA/networks.py')