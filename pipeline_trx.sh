#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p long
#SBATCH -t 100:00:00

##########################################################
########################## SETUP #########################
##########################################################

# Before running, may need to convert sh scripts to Unix with one of the following:
	# [On CLI] dos2unix <file_name.sh>
	# [On Notepad++] Edit > EOL Conversion > Unix (LF)

# Exit script upon encountering an error
set -e

# Source configuration
source local_config.sh

# Absolute paths to reference files
index_KT2440='/home/gcabebe/rnaseq/strain_KT2440/index_STAR'
gtf_KT2440='/home/gcabebe/rnaseq/strain_KT2440/refseqKT2440/Pseudomonas_putida_KT2440_110.gtf'

index_ecoli_k12='/home/gcabebe/rnaseq/strain_Ecoli_K12_MG1655/index_STAR'
gtf_ecoli_k12='/home/gcabebe/rnaseq/strain_Ecoli_K12_MG1655/ref/GTF_E_coli_str_k_12.gtf'

index_paerug_pao1='/home/gcabebe/rnaseq/strain_Paeruginosa_PAO1/index_STAR'
gtf_paerug_pao1='/home/gcabebe/rnaseq/strain_Paeruginosa_PAO1/ref_genome/GTF_Paeruginosa_PAO1_107.gtf'

index_pglucan='/home/gcabebe/rnaseq/strain_Pglucanolyticus_J12/index_STAR'
gtf_pglucan='/home/gcabebe/rnaseq/strain_Pglucanolyticus_J12/genomic.gtf'

index_saureus_hg003='/home/gcabebe/rnaseq/strain_Saureus/HG003/index_STAR'
gtf_saureus_hg003='/home/gcabebe/rnaseq/strain_Saureus/HG003/genomic_HG003.gtf'

index_saureus_hg001='/home/gcabebe/rnaseq/strain_Saureus/HG001/index_STAR'
gtf_saureus_hg001='/home/gcabebe/rnaseq/strain_Saureus/HG001/genomic_HG001.gtf'

index_saureus_atcc_51811='/home/gcabebe/rnaseq/strain_Saureus/ATCC_51811/index_STAR'
gtf_saureus_atcc_51811='/home/gcabebe/rnaseq/strain_Saureus/ATCC_51811/genomic_ATCC_51811.gtf'

index_saureus_newman='/home/gcabebe/rnaseq/strain_Saureus/Newman/index_STAR'
gtf_saureus_newman='/home/gcabebe/rnaseq/strain_Saureus/Newman/genomic_Newman.gtf'

index_saureus_rn4220='/home/gcabebe/rnaseq/strain_Saureus/RN4220/index_STAR'
gtf_saureus_rn4220='/home/gcabebe/rnaseq/strain_Saureus/RN4220/genomic_RN4220.gtf'

index_saureus_usa300='/home/gcabebe/rnaseq/strain_Saureus/USA300/index_STAR'
gtf_saureus_usa300='/home/gcabebe/rnaseq/strain_Saureus/USA300/genomic_USA300_TCH1516.gtf'

index_saureus_t145='/home/gcabebe/rnaseq/strain_Saureus/t145/index_STAR'
gtf_saureus_t145='/home/gcabebe/rnaseq/strain_Saureus/t145/genomic_t145.gtf'

index_saureus_jkd6008='/home/gcabebe/rnaseq/strain_Saureus/JKD6008/index_STAR'
gtf_saureus_jkd6008='/home/gcabebe/rnaseq/strain_Saureus/JKD6008/genomic_JKD6008.gtf'

# Determine STAR index and GTF based on organism
case "$ORGANISM_NAME" in
    "E coli K12")
		echo "Using E. coli reference files..."
        STAR_INDEX="$index_ecoli_k12"
        GTF_PATH="$gtf_ecoli_k12"
        ;;
    "Pseudomonas putida KT2440")
		echo "Using P. putida KT2440 reference files..."
        STAR_INDEX="$index_KT2440"
        GTF_PATH="$gtf_KT2440"
        ;;
    "Pseudomonas aeruginosa PAO1")
		echo "Using P. aeruginosa PAO1 reference files..."
        STAR_INDEX="$index_paerug_pao1"
        GTF_PATH="$gtf_paerug_pao1"
        ;;
	"Paenibacillus glucanolyticus J12")
		echo "Using P. glucanolyticus J12 reference files..."
		STAR_INDEX="$index_pglucan"
        GTF_PATH="$gtf_pglucan"
		;;
	"S aureus HG003")
		echo "Using S. aureus HG003 reference files..."
        STAR_INDEX="$index_saureus_hg003"
        GTF_PATH="$gtf_saureus_hg003"
		;;
	"S aureus HG001")
		echo "Using S. aureus HG001 reference files..."
        STAR_INDEX="$index_saureus_hg001"
        GTF_PATH="$gtf_saureus_hg001"
		;;
	"S aureus ATCC_51811")
		echo "Using S. aureus ATCC_51811 reference files..."
        STAR_INDEX="$index_saureus_atcc_51811"
        GTF_PATH="$gtf_saureus_atcc_51811"
		;;
	"S aureus Newman")
		echo "Using S. aureus Newman reference files..."
        STAR_INDEX="$index_saureus_newman"
        GTF_PATH="$gtf_saureus_newman"
		;;
	"S aureus RN4220")
		echo "Using S. aureus RN4220 reference files..."
        STAR_INDEX="$index_saureus_rn4220"
        GTF_PATH="$gtf_saureus_rn4220"
		;;
	"S aureus USA300")
		echo "Using S. aureus USA300 reference files..."
        STAR_INDEX="$index_saureus_usa300"
        GTF_PATH="$gtf_saureus_usa300"
		;;
	"S aureus t145")
		echo "Using S. aureus t145 reference files..."
        STAR_INDEX="$index_saureus_t145"
        GTF_PATH="$gtf_saureus_t145"
		;;
	"S aureus JKD6008")
		echo "Using S. aureus JKD6008 reference files..."
        STAR_INDEX="$index_saureus_jkd6008"
        GTF_PATH="$gtf_saureus_jkd6008"
		;;
    *)
        echo "Error: Unknown organism '$ORGANISM_NAME'"
        exit 1
        ;;
esac


##########################################################
###################### ENA DOWNLOAD ######################
##########################################################

cd "$ROOT_PATH"

FASTQ_RAW_DIR='fastq_raw'
mkdir -p $FASTQ_RAW_DIR
cd $FASTQ_RAW_DIR

# Run the download script within the same shell session
source "$ROOT_PATH/$DOWNLOAD_SCRIPT_NAME"

cd ..

##########################################################
##################### TRIMMING READS #####################
##########################################################

FASTQ_TRIM_DIR="fastq_trim"
mkdir -p "$FASTQ_TRIM_DIR"
cd "$FASTQ_TRIM_DIR"

declare -A seen

for fq in "$ROOT_PATH/$FASTQ_RAW_DIR"/*.{fastq,fq}.gz; do
    [[ -e "$fq" ]] || continue
    [[ -n "${seen[$fq]}" ]] && continue

    filename=$(basename "$fq")

    # Remove .gz first, then .fastq or .fq
    base_no_gz="${filename%.gz}"
    base_no_ext="${base_no_gz%.fastq}"
    base_no_ext="${base_no_ext%.fq}"

    # Paired-end
    if [[ "$base_no_ext" =~ _1$ ]]; then
        base="${base_no_ext%_1}"

		if [[ -f "$ROOT_PATH/$FASTQ_RAW_DIR/${base}_1.fastq.gz" ]]; then
            ext="fastq.gz"
        else
            ext="fq.gz"
        fi

        R1="$ROOT_PATH/$FASTQ_RAW_DIR/${base}_1.$ext"
        R2="$ROOT_PATH/$FASTQ_RAW_DIR/${base}_2.$ext"

        if [[ -f "$R2" ]]; then
            out1="${base}_paired_trimm1.fastq.gz"
            out2="${base}_paired_trimm2.fastq.gz"

            if [[ -f "$out1" && -f "$out2" ]]; then
                echo "[Skipping] PE already trimmed: $base"
                seen["$R1"]=1
                seen["$R2"]=1
                continue
            fi

            echo "[Trimming] PE detected: $base"

            java -jar /home/gcabebe/trx_tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
                -threads 4 -phred33 \
                "$R1" "$R2" \
                "$out1" \
                "${base}_unpaired_trimm1.fastq.gz" \
                "$out2" \
                "${base}_unpaired_trimm2.fastq.gz" \
                ILLUMINACLIP:/home/gcabebe/trx_tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True \
                LEADING:3 TRAILING:3 MINLEN:36

            seen["$R1"]=1
            seen["$R2"]=1
            continue
        fi
    fi

    # Single-end
    base="$base_no_ext"
    out="${base}_trim.fastq.gz"

    if [[ -f "$out" ]]; then
        echo "[Skipping] SE already trimmed: $base"
        seen["$fq"]=1
        continue
    fi

    echo "[Trimming] SE detected: $base"

    java -jar /home/gcabebe/trx_tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
        -threads 4 -phred33 \
        "$fq" \
        "$out" \
        ILLUMINACLIP:/home/gcabebe/trx_tools/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10:2:true \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    seen["$fq"]=1
done

cd ..


##########################################################
##################### FASTQC #####################
##########################################################

FASTQC_RESULTS='fastqc'
mkdir -p $FASTQC_RESULTS

if [ -f "$FASTQC_RESULTS/multiqc_report.html" ] || [ -f "$FASTQC_RESULTS/multiqc_report.zip" ]; then
    echo "[FastQC] MultiQC report present. Skipping FastQC + MultiQC."
else
	# Run FastQC on all fastq.gz files in the specified directory
	for filename in "$FASTQ_TRIM_DIR"/*.{fastq,fq}.gz; do
		# Run FASTQC on current file and extract to specified folder
		fastqc "$filename" --extract -o "$(pwd)/$FASTQC_RESULTS"
	done

	# run MultiQC in the directory containing fastqc results
	cd $FASTQC_RESULTS
	multiqc .

	cd ..
fi


##########################################################
##################### ALIGN READS WITH STAR #####################
##########################################################

ALIGN_DIR="align_star"
mkdir -p "$ALIGN_DIR"

echo "[Alignment] Running STAR..."

cd "$ALIGN_DIR"

declare -A seen

for fq in "$ROOT_PATH/$FASTQ_TRIM_DIR"/*.fastq.gz "$ROOT_PATH/$FASTQ_TRIM_DIR"/*.fq.gz; do
    [[ -e "$fq" ]] || continue
    [[ -n "${seen[$fq]}" ]] && continue

    filename=$(basename "$fq")

    # Trimmomatic's PE-mode orphaned/unpaired leftover reads are not real
    # samples - skip aligning them (and any already-existing SE-trimmed
    # variant of them from older runs)
    if [[ "$filename" =~ _unpaired_trimm[12] ]]; then
        seen["$fq"]=1
        continue
    fi

    #################################################
    # Paired-end detection
    #################################################

    if [[ "$filename" =~ _1\.fastq\.gz$ ]]; then

        base="${filename%_1.fastq.gz}"
        R1="$ROOT_PATH/$FASTQ_TRIM_DIR/${base}_1.fastq.gz"
        R2="$ROOT_PATH/$FASTQ_TRIM_DIR/${base}_2.fastq.gz"

    elif [[ "$filename" =~ _paired_trimm1\.fastq\.gz$ ]]; then

        base="${filename%_paired_trimm1.fastq.gz}"
        R1="$ROOT_PATH/$FASTQ_TRIM_DIR/${base}_paired_trimm1.fastq.gz"
        R2="$ROOT_PATH/$FASTQ_TRIM_DIR/${base}_paired_trimm2.fastq.gz"

    else
        R1=""
        R2=""
    fi

    #################################################
    # Run paired-end alignment
    #################################################

    if [[ -n "$R1" && -f "$R2" ]]; then

        outbam="${base}_Aligned.sortedByCoord.out.bam"

        if [[ -f "$outbam" ]]; then
            echo "[Skipping] PE already aligned: $base"

            seen["$R1"]=1
            seen["$R2"]=1
            continue
        fi

        echo "[Alignment] PE detected: $base"

        if ! STAR \
            --genomeDir "$STAR_INDEX" \
            --readFilesIn "$R1" "$R2" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${base}_" \
            --outFilterMultimapNmax 20 \
            --outReadsUnmapped Fastx \
            --outSAMtype BAM SortedByCoordinate \
            --twopassMode Basic \
            --runThreadN 2 \
            --limitBAMsortRAM 2793708443
        then
            echo "[Alignment] ERROR: STAR failed for $base - skipping, will retry on next run"
            echo "$(date '+%F %T') PE $base" >> "$ROOT_PATH/$ALIGN_DIR/failed_alignments.log"
            rm -f "$outbam"

            seen["$R1"]=1
            seen["$R2"]=1
            continue
        fi

        seen["$R1"]=1
        seen["$R2"]=1

        continue
    fi

    #################################################
    # Single-end alignment
    #################################################

    base="${filename%.fastq.gz}"
    base="${base%.fq.gz}"

    outbam="${base}_Aligned.sortedByCoord.out.bam"

    if [[ -f "$outbam" ]]; then
        echo "[Skipping] SE already aligned: $base"

        seen["$fq"]=1
        continue
    fi

    echo "[Alignment] SE detected: $base"

    if ! STAR \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$fq" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${base}_" \
        --outFilterMultimapNmax 20 \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --runThreadN 2 \
        --limitBAMsortRAM 2793708443
    then
        echo "[Alignment] ERROR: STAR failed for $base - skipping, will retry on next run"
        echo "$(date '+%F %T') SE $base" >> "$ROOT_PATH/$ALIGN_DIR/failed_alignments.log"
        rm -f "$outbam"
    fi

    seen["$fq"]=1

done

cd ..

##########################################################
############# FILTER OUT LOW QUALITY READS ###############
##########################################################

ALIGN_HQ_DIR='indexed_HQ'

mkdir -p $ALIGN_HQ_DIR
cd $ALIGN_HQ_DIR

for SAMPLE in "$ROOT_PATH"/"$ALIGN_DIR"/*Aligned.sortedByCoord.out.bam
do
	# Skip stray alignments of Trimmomatic's orphaned/unpaired leftover reads
	# from older runs - these aren't real samples
	if [[ "$(basename "$SAMPLE")" =~ _unpaired_trimm[12] ]]; then
		continue
	fi

    # shave down filename to just its SRR # (eg. 'SRR6012666_2.fastq.gz' to just 'SRR6012666')
	base_file_name=$(basename ${SAMPLE%%.*})
	sample_base=${base_file_name%_Aligned}

	hq_bam="${base_file_name}_high_mapq_reads.bam"
	fwd_bam="${base_file_name}_high_mapq_reads_forward.sorted.bam"
	rev_bam="${base_file_name}_high_mapq_reads_reverse.sorted.bam"

	# PE dUTP/TruSeq stranded samples get split into forward/reverse strand BAMs;
	# SE samples (no _paired_trimm1 fastq) just get the quality-filtered BAM.
	if [[ -f "$ROOT_PATH/$FASTQ_TRIM_DIR/${sample_base}_paired_trimm1.fastq.gz" ]]; then

		if [[ -f "$fwd_bam" && -f "$rev_bam" ]]; then
			echo "[Skipping] PE stranded already filtered: ${sample_base}"
			continue
		fi

		echo "[Filtering Aligned Reads - PE stranded] Currently on ${sample_base}"

		tmp1="${base_file_name}_tmp1.bam"
		tmp2="${base_file_name}_tmp2.bam"
		tmp3="${base_file_name}_tmp3.bam"
		tmp4="${base_file_name}_tmp4.bam"

		if ! (
			set -e
			samtools view -h -b -q 20 "$SAMPLE" > "$hq_bam"

			# Forward strand reads
			samtools view -b -f 128 -F 16 "$hq_bam" > "$tmp1"
			samtools view -b -f 80 "$hq_bam" > "$tmp2"
			samtools merge -f "${base_file_name}_high_mapq_reads_forward.bam" "$tmp1" "$tmp2"
			samtools sort -o "$fwd_bam" "${base_file_name}_high_mapq_reads_forward.bam"

			# Reverse strand reads
			samtools view -b -f 144 "$hq_bam" > "$tmp3"
			samtools view -b -f 64 -F 16 "$hq_bam" > "$tmp4"
			samtools merge -f "${base_file_name}_high_mapq_reads_reverse.bam" "$tmp3" "$tmp4"
			samtools sort -o "$rev_bam" "${base_file_name}_high_mapq_reads_reverse.bam"
		); then
			echo "[Filtering] ERROR: PE stranded filtering failed for ${sample_base} - skipping, will retry on next run"
			echo "$(date '+%F %T') PE_FILTER ${sample_base}" >> "$ROOT_PATH/$ALIGN_HQ_DIR/failed_filtering.log"
			rm -f "$hq_bam" "$fwd_bam" "$rev_bam" \
				"${base_file_name}_high_mapq_reads_forward.bam" \
				"${base_file_name}_high_mapq_reads_reverse.bam"
			rm -f "$tmp1" "$tmp2" "$tmp3" "$tmp4"
			continue
		fi

		# tmp files and the unsorted merged BAMs were only intermediates for the
		# split; hq_bam (the unstranded total) is kept for featureCounts below
		rm -f "$tmp1" "$tmp2" "$tmp3" "$tmp4" \
			"${base_file_name}_high_mapq_reads_forward.bam" \
			"${base_file_name}_high_mapq_reads_reverse.bam"

	else

		if [[ -f "$hq_bam" ]]; then
			echo "[Skipping] SE already filtered: ${sample_base}"
			continue
		fi

		echo "[Filtering Aligned Reads - SE] Currently on ${sample_base}"

		if ! samtools view -h -b -q 20 "$SAMPLE" > "$hq_bam"; then
			echo "[Filtering] ERROR: SE filtering failed for ${sample_base} - skipping, will retry on next run"
			echo "$(date '+%F %T') SE_FILTER ${sample_base}" >> "$ROOT_PATH/$ALIGN_HQ_DIR/failed_filtering.log"
			rm -f "$hq_bam"
			continue
		fi

	fi
done

cd ..

##########################################################
################ GENERATE COVERAGE FILES ##################
##########################################################

COVERAGE_DIR='genomecov'

mkdir -p "$COVERAGE_DIR"
cd "$COVERAGE_DIR"

echo "[Bedtools Coverage] Generating coverage files..."

# Stranded coverage for PE samples (forward/reverse split BAMs)
for fwd_bam in "$ROOT_PATH/$ALIGN_HQ_DIR"/*_high_mapq_reads_forward.sorted.bam
do
	[[ -e "$fwd_bam" ]] || continue

	base=$(basename "$fwd_bam" "_forward.sorted.bam")
	rev_bam="$ROOT_PATH/$ALIGN_HQ_DIR/${base}_reverse.sorted.bam"

	fwd_out="${base}.forward.coverage"
	rev_out="${base}.reverse.coverage"

	if [[ -f "$fwd_out" && -f "$rev_out" ]]; then
		echo "[Skipping] Stranded coverage already exists: ${base}"
		continue
	fi

	echo "[Bedtools Coverage] Processing ${base} (stranded)"
	bedtools genomecov -d -ibam "$fwd_bam" > "$fwd_out"
	bedtools genomecov -d -ibam "$rev_bam" > "$rev_out"
done

# Unstranded coverage for SE samples (PE samples already got stranded coverage above)
for bam in "$ROOT_PATH/$ALIGN_HQ_DIR"/*_high_mapq_reads.bam
do
	[[ -e "$bam" ]] || continue

	base=$(basename "$bam" .bam)

	if [[ -f "$ROOT_PATH/$ALIGN_HQ_DIR/${base}_forward.sorted.bam" ]]; then
		continue
	fi

	out="${base}.coverage"

	if [[ -f "$out" ]]; then
		echo "[Skipping] Coverage already exists: ${base}"
		continue
	fi

	echo "[Bedtools Coverage] Processing ${base} (unstranded)"
	bedtools genomecov -d -ibam "$bam" > "$out"
done

cd ..

##########################################################
############## GENERATE READ COUNT TABLE #################
##########################################################

READ_COUNTS_DIR='read_counts'

mkdir -p $READ_COUNTS_DIR
cd $READ_COUNTS_DIR

echo "[Readcounts table] Running featureCounts..."

# bam files are of several strains but we will use only KT2440
# Only the SE quality-filtered BAMs remain named *_high_mapq_reads.bam directly in
# $ALIGN_HQ_DIR; PE samples were split into forward/reverse strand BAMs above, so
# count those instead to avoid mixing totals with per-strand splits.
featureCounts -a "$GTF_PATH" -o "./featureCounts_results_${STUDY_ID}.txt" \
	"${ROOT_PATH}/${ALIGN_HQ_DIR}"/*_high_mapq_reads.bam \
	"${ROOT_PATH}/${ALIGN_HQ_DIR}"/*_high_mapq_reads_forward.sorted.bam \
	"${ROOT_PATH}/${ALIGN_HQ_DIR}"/*_high_mapq_reads_reverse.sorted.bam \
	-t CDS
