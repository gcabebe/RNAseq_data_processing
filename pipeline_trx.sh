#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p long
#SBATCH -t 100:00:00

##########################################################
########################## SETUP #########################
##########################################################

# Before running, may need to convert sh scripts to Unix with - dos2unix <file_name.sh>

# Exit script upon encountering an error
set -e

# Source configuration
source local_config.sh

# Absolute paths to reference files
index_KT2440='/home/gcabebe/rnaseq/strain_KT2440/index_STAR'
gtf_KT2440='/home/gcabebe/rnaseq/strain_KT2440/refseqKT2440/Pputida_gtf_COPY.gtf'

index_ecoli_k12='/home/gcabebe/rnaseq/strain_Ecoli_K12_MG1655/STAR_index'
gtf_ecoli_k12='/home/gcabebe/rnaseq/strain_Ecoli_K12_MG1655/ref/GTF_E_coli_str_k_12.gtf'

index_paerug_pao1='/home/gcabebe/rnaseq/strain_Paeruginosa_PAO1/index_STAR'
gtf_paerug_pao1='/home/gcabebe/rnaseq/strain_Paeruginosa_PAO1/ref_genome/GTF_Paeruginosa_PAO1_107.gtf'

# Determine STAR index and GTF based on organism
case "$ORGANISM_NAME" in
    "E coli K12")
        STAR_INDEX="$index_ecoli_k12"
        GTF_PATH="$gtf_ecoli_k12"
        ;;
    "Pseudomonas putida KT2440")
        STAR_INDEX="$index_KT2440"
        GTF_PATH="$gtf_KT2440"
        ;;
    "Pseudomonas aeruginosa PAO1")
        STAR_INDEX="$index_paerug_pao1"
        GTF_PATH="$gtf_paerug_pao1"
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

##########################################################
##################### TRIMMING READS #####################
##########################################################

# run Trimmomatic on all fastq.gz files in the directory
# first recursively go into subdirectories to find fastq files
# https://stackoverflow.com/questions/107995/how-do-you-recursively-unzip-archives-in-a-directory-and-its-subdirectories-from
find . -name "*.fastq.gz" | while read filename
do 
	echo $(basename ${filename%_*})
	# if paired-end reads, run the Trimmomatic command for PE 
	# https://stackoverflow.com/questions/50321291/find-if-filename-contains-a-certain-string-in-bash-script
	if [[ "$filename" == *"_1.fastq.gz"* ]];then
		java -jar /home/gcabebe/trx_tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE $filename "${filename%_*}_2.fastq.gz" "$(basename ${filename%_*}_paired_trimm1.fastq.gz)" "$(basename ${filename%_*}_unpaired_trimm1.fastq.gz)" "$(basename ${filename%_*}_paired_trimm2.fastq.gz)" "$(basename ${filename%_*}_unpaired_trimm2.fastq.gz)"  ILLUMINACLIP:/home/gcabebe/trx_tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
	# skipped the PE read ending with "_2.fastq.gz" because already processed with "_1.fq.gz"
	elif [[ "$filename" == *"_2.fastq.gz"* ]];then
		echo "Skipping $file ..."
	# run Trimmomatic command for SE reads (file naming is a little messed up, fix later)
	else
		java -jar /home/gcabebe/trx_tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 $filename "`basename $filename`.trim.fastq.gz" ILLUMINACLIP:/home/gcabebe/trx_tools/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  	fi
done

cd ..

# store all reads into new directory 'fastq_trim'
FASTQ_TRIM_DIR='fastq_trim'
mkdir -p $FASTQ_TRIM_DIR
cd $FASTQ_RAW_DIR
find . -maxdepth 1 \( -name "*_paired_*" -o -name "*.trim.*" \) -exec mv -t "$ROOT_PATH/$FASTQ_TRIM_DIR" {} +

cd ..

##########################################################
##################### FASTQC #####################
##########################################################

FASTQC_RESULTS='fastqc'
mkdir -p $FASTQC_RESULTS

# Run FastQC on all fastq.gz files in the specified directory
for filename in "$FASTQ_TRIM_DIR"/*.fastq.gz; do
    # Run FASTQC on current file and extract to specified folder
    fastqc "$filename" --extract -o "$(pwd)/$FASTQC_RESULTS"
done

# run MultiQC in the directory containing fastqc results
cd $FASTQC_RESULTS
multiqc .

cd ..

##########################################################
##################### ALIGN READS WITH STAR #####################
##########################################################

ALIGN_DIR='align_star'

mkdir -p $ALIGN_DIR
cd $ALIGN_DIR

for SAMPLE in "$ROOT_PATH/$FASTQ_TRIM_DIR"/*.fastq.gz
do
    base_file_name=$(basename ${SAMPLE%%.*})
  
    # Check if output BAM file already exists
    if [ ! -e "${ROOT_PATH}/${ALIGN_DIR}/${base_file_name}_Aligned.sortedByCoord.out.bam" ]; then
        echo "Calculating alignment for $SAMPLE..."
        STAR --genomeDir "$STAR_INDEX" --readFilesIn $SAMPLE --readFilesCommand zcat --outFileNamePrefix ${base_file_name}_ --outFilterMultimapNmax 1 --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --twopassMode Basic --runThreadN 2 --limitBAMsortRAM 1700000000
    else
        echo "Alignment already exists for $SAMPLE. Skipping..."
    fi
done

cd ..

##########################################################
##################### FILTER OUT LOW QUALITY READS #####################
##########################################################

ALIGN_HQ_DIR='indexed_HQ'

mkdir -p $ALIGN_HQ_DIR
cd $ALIGN_HQ_DIR

# BAM_FILES=`ls -m ${ROOT_PATH}/${ALIGN_DIR}/*Aligned.sortedByCoord.out.bam | sed 's/,//g'`

for SAMPLE in "$ROOT_PATH"/"$ALIGN_DIR"/*Aligned.sortedByCoord.out.bam
do
    # shave down filename to just its SRR # (eg. 'SRR6012666_2.fastq.gz' to just 'SRR6012666')
	base_file_name=$(basename ${SAMPLE%%.*})
	project_name=${base_file_name%%_*}
	echo "${project_name}"
    samtools view -h -b -q 20 $SAMPLE > ${base_file_name}_high_mapq_reads.bam
done

cd ..

##########################################################
##################### GENERATE READ COUNT TABLE #####################
##########################################################

READ_COUNTS_DIR='read_counts'

mkdir -p $READ_COUNTS_DIR
cd $READ_COUNTS_DIR

# bam files are of several strains but we will use only KT2440
featureCounts -a "$GTF_PATH" -o "./featureCounts_results_${STUDY_ID}.txt" "${ROOT_PATH}/${ALIGN_HQ_DIR}"/*.bam -t CDS	

####################################################################################
##################### GENERATE TPM, VST, and DESeq2 NORMALIZED TABLES ##############
####################################################################################

#python RawRCs_to_NormTables.py $STUDY_ID $TPM_GENE_REF_PATH $ROOT_PATH "$ROOT_PATH/SraRunTable.csv"