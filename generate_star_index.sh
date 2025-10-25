#!/bin/bash
#SBATCH -N 2
#SBATCH -n 4

script=$0
# RefSeq directory
REF_DIR=$1 # /home/gcabebe/rnaseq/refseqKT2440/refseq/Pseudomonas_putida_KT2440_110.fna
# GTF file directory
GTF_DIR=$2 # /home/gcabebe/rnaseq/refseqKT2440/gtf/Pseudomonas_putida_KT2440_110.gtf
# directory where STAR index will be placed
INDEX_DIR=$3 # /home/gcabebe/rnaseq/STARindex

# run STAR in"genomeGenerate" mode
STAR --runMode genomeGenerate --runThreadN 1 --genomeDir ${INDEX_DIR} --genomeFastaFiles ${REF_DIR} --sjdbGTFfile ${GTF_DIR} --sjdbOverhang 49 --quantMode GeneCounts --genomeSAindexNbases 10
