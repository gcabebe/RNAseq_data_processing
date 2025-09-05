#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p long
#SBATCH -t 120:00:00

set -e

# Required inputs
PROJECT_ID=$1
RUNLIST_FILE=$2

# === Set root and source config ===
ROOT_PATH=$(pwd)/projects/$PROJECT_ID
mkdir -p "$ROOT_PATH"

source /home/gcabebe/rnaseq_mining_scripts/config.sh

FASTQ_RAW_DIR="$ROOT_PATH/fastq_raw"
FASTQ_TRIM_DIR="$ROOT_PATH/fastq_trim"
FASTQC_RESULTS="$ROOT_PATH/fastqc"
ALIGN_DIR="$ROOT_PATH/align_star"
ALIGN_HQ_DIR="$ROOT_PATH/indexed_HQ"
READ_COUNTS_DIR="$ROOT_PATH/read_counts"

# Reference files
TRIMMOMATIC_DIR=/home/gcabebe/trx_tools/Trimmomatic-0.39
ADAPTER_PE=$TRIMMOMATIC_DIR/adapters/TruSeq3-PE.fa
ADAPTER_SE=$TRIMMOMATIC_DIR/adapters/TruSeq3-SE.fa

mkdir -p "$FASTQ_RAW_DIR" "$FASTQ_TRIM_DIR" "$FASTQC_RESULTS" "$ALIGN_DIR" "$ALIGN_HQ_DIR" "$READ_COUNTS_DIR"

##########################################################
##################### DOWNLOAD DATA ######################
##########################################################
# Debug #1
echo "Contents of $RUNLIST_FILE:"
cat "$RUNLIST_FILE"

while read SRR; do
    [[ -z "$SRR" ]] && continue  # Skip empty lines
    [[ ! "$SRR" =~ ^(SRR|DRR|ERR) ]] && continue  # Skip malformed entries, only take samples with SRR/DRR/ERR accessions

    echo "Fetching: $SRR"
    prefetch "$SRR" --output-directory "$FASTQ_RAW_DIR" || { echo "Prefetch failed for $SRR"; exit 1; }

    echo "Running fasterq-dump for $SRR..."
    fasterq-dump "$FASTQ_RAW_DIR/"$SRR"/$SRR.sra" --outdir "$FASTQ_RAW_DIR" --split-files || { echo "fasterq-dump failed for $SRR"; exit 1; }

done < "$RUNLIST_FILE"



##########################################################
##################### TRIMMING READS #####################
##########################################################

cd "$FASTQ_RAW_DIR"

# Abort job if fastq files not found
if ! ls *.fastq 1> /dev/null 2>&1; then
    echo "No FASTQ files found — aborting Trimmomatic step"
    exit 1
fi

PAIRED_BAMS=()
SINGLE_BAMS=()


for filename in *.fastq; do
    base=$(basename "$filename" .fastq)
    sample_id=${base%_1}  # Remove _1 or _2 suffix if present

    if [[ "$filename" == *_1.fastq ]]; then
        R1="$filename"
        R2="${base/_1/_2}.fastq"

        # Run Trimmomatic for PE
        java -jar $TRIMMOMATIC_DIR/trimmomatic-0.39.jar PE \
            "$R1" "$R2" \
            "$FASTQ_TRIM_DIR/${sample_id}_paired_1.fastq.gz" "$FASTQ_TRIM_DIR/${sample_id}_unpaired_1.fastq.gz" \
            "$FASTQ_TRIM_DIR/${sample_id}_paired_2.fastq.gz" "$FASTQ_TRIM_DIR/${sample_id}_unpaired_2.fastq.gz" \
            ILLUMINACLIP:$ADAPTER_PE:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

        # Track BAM file after alignment
        PAIRED_BAMS+=("$ALIGN_HQ_DIR/${sample_id}_1_high_mapq_reads.bam")

    elif [[ "$filename" != *_2.fastq ]]; then
        # Run Trimmomatic for SE
        java -jar $TRIMMOMATIC_DIR/trimmomatic-0.39.jar SE -phred33 "$filename" \
            "$FASTQ_TRIM_DIR/${base}.trim.fastq.gz" \
            ILLUMINACLIP:$ADAPTER_SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

        # Track BAM file after alignment
        SINGLE_BAMS+=("$ALIGN_HQ_DIR/${base}_high_mapq_reads.bam")
    fi
done


##########################################################
##################### FASTQC #############################
##########################################################

for file in "$FASTQ_TRIM_DIR"/*.fastq.gz; do
    fastqc "$file" --extract -o "$FASTQC_RESULTS"
done

multiqc "$FASTQC_RESULTS" -o "$FASTQC_RESULTS"

##########################################################
##################### ALIGNMENT (STAR) ###################
##########################################################

cd "$ALIGN_DIR"

for file in "$FASTQ_TRIM_DIR"/*_paired_1.fastq.gz; do
    base=$(basename "$file" _paired_1.fastq.gz)
    STAR --genomeDir $STAR_INDEX \
        --readFilesIn "$FASTQ_TRIM_DIR/${base}_paired_1.fastq.gz" "$FASTQ_TRIM_DIR/${base}_paired_2.fastq.gz" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${base}_" \
        --outFilterMultimapNmax 1 \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --runThreadN 4 \
        --limitBAMsortRAM 1700000000
done

##########################################################
##################### FILTER HIGH MAPQ ###################
##########################################################

cd "$ALIGN_HQ_DIR"

for bam in "$ALIGN_DIR"/*Aligned.sortedByCoord.out.bam; do
    base=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
    samtools view -h -b -q 20 "$bam" > "${base}_high_mapq_reads.bam"
done

##########################################################
##################### FEATURE COUNTS #####################
##########################################################

cd "$READ_COUNTS_DIR"

# Run featureCounts for paired-end reads
if [ ${#PAIRED_BAMS[@]} -gt 0 ]; then
    featureCounts -a "$GTF_PATH" -o "featureCounts_results_${PROJECT_ID}_paired.txt" \
        -t CDS -p "${PAIRED_BAMS[@]}"
fi

# Run featureCounts for single-end reads
if [ ${#SINGLE_BAMS[@]} -gt 0 ]; then
    featureCounts -a "$GTF_PATH" -o "featureCounts_results_${PROJECT_ID}_single.txt" \
        -t CDS "${SINGLE_BAMS[@]}"
fi


#featureCounts -a $GTF_PATH -o "featureCounts_results_${PROJECT_ID}.txt" \
#    "$ALIGN_HQ_DIR"/*.bam -t CDS

##########################################################
########## RUN PYTHON NORMALIZATION SCRIPT ################
##########################################################

python RawRCs_to_NormTables.py "$PROJECT_ID" "$GENE_LENGTHS_REF_PATH" "$ROOT_PATH"

