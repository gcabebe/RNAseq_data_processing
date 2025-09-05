#!/usr/bin/env python3

import re
import sys
import os
import warnings
import itertools
import subprocess
import numpy as np
import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

# Ignore all warnings
warnings.filterwarnings("ignore")

def download_sra_metadata(project_id, output_path):
    """
    Downloads the SraRunTable.csv using pysradb and saves it to output_path.
    """
    try:
        subprocess.run([
            "pysradb", "metadata", project_id,
            "--desc", "--expand", "-o", output_path
        ], check=True)
        print(f"Successfully downloaded SRA metadata to {output_path}")
    except subprocess.CalledProcessError:
        raise RuntimeError(f"Failed to download metadata for {project_id} using pysradb.")


def infer_condition_column(metadata_df):
    """
    Attempt to automatically detect the correct condition/treatment column
    based on number of unique non-null values and common names.
    """
    exclude_cols = {"Run", "Sample Name", "LibraryName", "avgLength", "Bases", "spots", "Experiment", "BioProject"}

    # Drop excluded and all-null columns
    candidates = metadata_df.drop(columns=[col for col in exclude_cols if col in metadata_df.columns])
    candidates = candidates.dropna(axis=1, how='all')

    # Filter to columns with low cardinality (2-10 unique values), ignoring 'NA'
    candidate_scores = {}
    for col in candidates.columns:
        non_na = candidates[col].dropna()
        n_unique = non_na.nunique()
        if 2 <= n_unique <= 10:
            candidate_scores[col] = n_unique

    if not candidate_scores:
        raise ValueError("No suitable condition column found in metadata!")

    # Prefer columns named like 'condition', 'treatment', etc.
    priority_cols = ['condition', 'treatment', 'group', 'factor']
    for col in priority_cols:
        if col in candidate_scores:
            print(f"Automatically selected condition column: '{col}'")
            return col

    # Otherwise pick the one with fewest unique values (simplest grouping)
    best_col = min(candidate_scores, key=candidate_scores.get)
    print(f"⚠️ Using best guess for condition column: '{best_col}'")
    return best_col


def get_readcounts(project_id, metadata_path):
    """
    Load and filter read counts and metadata.
    """
    # Read counts
    readcounts_path = f"{root_path}/read_counts/featureCounts_results_{project_id}.txt"
    readcounts = pd.read_csv(readcounts_path, sep="\t", skiprows=1).set_index("Geneid")
    readcounts = readcounts.iloc[:, 6:]  # Drop annotation cols
    readcounts = readcounts[[col for col in readcounts.columns if not re.search(r"paired|unpaired", col)]]

    # Metadata
    metadata = pd.read_csv(metadata_path, dtype=str)
    condition_col = infer_condition_column(metadata)

    metadata = metadata.dropna(subset=[condition_col])
    run_ids = metadata['Run'].astype(str)
    run_to_condition = metadata.set_index('Run')[condition_col].to_dict()

    # Standardize readcount colnames to SRR IDs
    readcounts.columns = [re.search(r"(SRR\d+)", col).group(1) if re.search(r"(SRR\d+)", col) else col for col in readcounts.columns]
    filtered_readcounts = readcounts.loc[:, readcounts.columns.isin(run_ids)]
    filtered_readcounts.columns = [run_to_condition.get(srr, srr) for srr in filtered_readcounts.columns]

    return filtered_readcounts


def tpm(counts, lengths):
    rate = counts / lengths
    return rate / sum(rate) * 1e6


def generate_tpm_table(tpm_gene_ref_path, readcounts, project_id):
    genes = pd.read_csv(tpm_gene_ref_path).drop(columns=["Unnamed: 0"]).set_index("gene_id")
    tpms = readcounts.apply(lambda x: tpm(x, genes['length']), axis=0)
    tpms.to_csv(f"{root_path}/{project_id}_tpm_counts.csv", index=True)
    print("✅ TPM table saved.")


def preprocess_for_deseq(readcounts):
    readcounts = readcounts.loc[:, (readcounts != 0).any(axis=0)]
    averaged_readcounts = pd.DataFrame(index=readcounts.index)
    unique_samples = list(dict.fromkeys(readcounts.columns))

    for sample in unique_samples:
        sample_cols = readcounts.loc[:, readcounts.columns == sample]
        averaged_cols = []
        for i in range(0, len(sample_cols.columns), 4):
            if i + 1 < len(sample_cols.columns):
                avg_col = sample_cols.iloc[:, i:i+4].mean(axis=1).round(0).astype(int)
            else:
                avg_col = sample_cols.iloc[:, i].round(0).astype(int)
            averaged_cols.append(avg_col)

        averaged_sample_df = pd.concat(averaged_cols, axis=1)
        averaged_sample_df.columns = [sample] * len(averaged_sample_df.columns)
        averaged_readcounts = pd.concat([averaged_readcounts, averaged_sample_df], axis=1)

    condition_time = [col.rsplit('_', 1)[0] for col in averaged_readcounts.columns]
    sample_info = pd.DataFrame({'condition': condition_time}, index=averaged_readcounts.columns)
    return averaged_readcounts, sample_info


def do_deseq_analys(readcounts):
    readcounts, sample_info = preprocess_for_deseq(readcounts)

    genes_to_keep = readcounts.columns[readcounts.sum(axis=0) >= 10]
    readcounts = readcounts[genes_to_keep]

    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=readcounts.T,
        metadata=sample_info,
        design_factors="condition",
        refit_cooks=True,
        inference=inference
    )

    dds.vst_fit(use_design=False)
    vst_counts = dds.vst_transform()
    vst_df = pd.DataFrame(vst_counts.T, index=readcounts.index, columns=readcounts.columns)
    vst_df.to_csv(f"{root_path}/{project_id}_rlog_normalized_counts.csv", index=True)
    print("✅ Saved rlog/VST-normalized counts.")

    dds.deseq2()

    print("📊 Generating DESeq2 pairwise comparisons...")
    for A, B in itertools.combinations(dds.obs["condition"].unique(), 2):
        print(f"  ➤ Comparing {B} vs {A}")
        ds = DeseqStats(dds, contrast=["condition", B, A], inference=inference)
        ds.summary()
        ds.results_df.to_csv(f"{root_path}/DGE_{B}_vs_{A}_{project_id}.csv", index=True)

    print("🎉 All DESeq2 results saved.")


# === MAIN ===

project_id = sys.argv[1]
tpm_gene_ref_path = sys.argv[2]
root_path = sys.argv[3]

metadata_path = f"{project_id}/SraRunTable.csv"
download_sra_metadata(project_id, metadata_path)

readcounts = get_readcounts(project_id, metadata_path)
generate_tpm_table(tpm_gene_ref_path, readcounts, project_id)
do_deseq_analys(readcounts)
