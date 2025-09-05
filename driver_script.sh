#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4

set -e  # Exit on error

# === Load config ===
source /home/gcabebe/rnaseq_mining_scripts/config.sh

# === Create folders ===
mkdir -p metadata runlists

# === Summary log file ===
SUMMARY_LOG="rnaseq_mining_output_summary.txt"
> "$SUMMARY_LOG"  # Clear if exists

# === Step 1: Download runinfo from NCBI ===
echo "Fetching metadata for: $ESEARCH_QUERY"
esearch -db sra -query "$ESEARCH_QUERY" | efetch -format runinfo > "$EFETCH_FILENAME"

# === Step 2: Extract unique BioProjects ===
cut -d',' -f1 --complement "$EFETCH_FILENAME" | cut -d',' -f1 > tmp_bioprojects.txt  # crude but safe
grep -oE 'PRJ[EDN][A-Z0-9]+' "$EFETCH_FILENAME" | sort -u > bioprojects.txt

# === Step 3: Loop over each BioProject and extract runs ===
while read -r PROJECT_ID; do
    echo "Processing project: $PROJECT_ID"
    METADATA_CSV="metadata/${PROJECT_ID}_metadata.csv"
    RUNLIST_FILE="runlists/${PROJECT_ID}.runs.txt"

    # Fetch run-level metadata
    pysradb metadata "$PROJECT_ID" --desc > "$METADATA_CSV"

    # VERSION 1: Extract Run column (skip header)
    #awk -F',' 'NR==1 {for(i=1;i<=NF;i++) if($i=="Run") col=i} NR>1 && $col != "" {print $col}' "$METADATA_CSV" > "$RUNLIST_FILE"
	# VERSION 2: Extract only SRR/ERR/DRR run accessions from the metadata file
	grep -oE '\b(SRR|ERR|DRR)[0-9]+\b' "$METADATA_CSV" | sort -u > "$RUNLIST_FILE"

	
    echo "Runlist for $PROJECT_ID:"
    cat "$RUNLIST_FILE"

    # Submit SLURM job
    sbatch pipeline.sh "$PROJECT_ID" "$RUNLIST_FILE"

done < bioprojects.txt

# === Step 4: Monitor SLURM jobs and update summary ===
echo "Waiting for pipeline jobs to finish..."
while true; do
    echo "=== Output Summary Update $(date) ===" >> "$SUMMARY_LOG"
    find . -maxdepth 1 -name "*.out" -exec bash -c 'echo "=== {} ==="; head -1 "{}"; tail -5 "{}"; echo' \; >> "$SUMMARY_LOG"

    # Exit if no pipeline jobs are running
    if ! squeue -u "$USER" | grep -q pipeline; then
        echo "All pipeline jobs finished." >> "$SUMMARY_LOG"
        break
    fi
    sleep 300  # Wait 5 min
done

echo "Driver script complete."
