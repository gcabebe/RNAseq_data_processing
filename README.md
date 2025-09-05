# TODO!!
- Check files under - /home/gcabebe/rnaseq_mining_scripts
- Check/edit environment.yml file

# Summary
This repository provides an executable script for pre-processing individual RNA-seq studies into read counts per gene tables. Codes assume the use of a high power computing (HPC) cluster and Simple Linux Utility for Resource Management (SLURM), but can be modified for running on a local computer.

## Example overview
The following RNA-seq project comes from ... (links to ENA and NCBI SRA page). In brief, they tested the gene expression of _P. putida_ in (control, ?? and ???). From the raw FASTQ files provided, this pipeline will process the data to create the following:

Raw read counts table:

(example readcounts table here)

DESeq2-generated differential gene expression (DGE) table:

(example DGE table here)

# Setup

## Files to Download
For each organism, you will need the files listed below. For _P. putida_, _P. aeruginosa_, and _E. coli_, these files are already provided within the folder ```FOLDER NAME HERE!!!```.
- Genome sequence (FASTA)
- Gene lengths (CSV)
  - Columns: gene_id, length
- GFF
- GTF

## HPC Requirements
Anaconda is essential for downloading packages and creating environments. Follow the [Anaconda documentation](/docs/getting-started/anaconda/install#verify-your-install) to install the Linux version of Anaconda.

Create a new conda environment by running the following:
```
conda env create -f environment.yml 
INCLUDE THESE LIBRARIES - numpy pandas pydeseq2 pysradb
(Below line removed since we are using environment.yml instead):
conda install -c conda-forge trimmomatic fastqc multiqc samtools qualimap subread
```

This will install the Linux and Python libraries that are used in this pipeline automatically, instead of having to install each package individually.

The following are downloaded from the source and need to be updated in the PATH environment (ie. on the .bashrc script):
- [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
- [STAR](https://github.com/alexdobin/STAR)


## Modifying the Config File


# Running the Script
Once you have generated your STAR index, modified