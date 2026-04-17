> [!NOTE]
> This repository is still under construction and is being updated week-to-week.

# Summary
This repository provides an executable script for pre-processing individual RNA-seq studies into read counts per gene tables. Codes assume the use of a high performance computing (HPC) cluster and Simple Linux Utility for Resource Management (SLURM), but can be modified for running on a local computer.

## Example dataset
The following RNA-seq project comes from a lab at the Technical University of Denmark. In brief, they tested the gene expression of _P. putida_ on non-trivial carbon sources (ie. citrate, ferulic acid, serine, and glucose) using transcriptomics and genome-scale modeling ([D'Arrigo et al 2019](https://doi.org/10.1111/1758-2229.12704)). Their data is publicly available and can be accessed on [NCBI](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP157937&o=acc_s%3Aa) and the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/view/srp157937).

---
# Setup

## Files to download
For each organism, you will need the files listed below. For [_P. putida_](https://www.pseudomonas.com/strain/show/110), [_P. aeruginosa_](https://www.pseudomonas.com/strain/show?id=107), and [_E. coli_](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/), these files are already provided within this repository for _P. putida_ under the ```data_P_putida``` folder.
- Genome sequence (FASTA) - nucleotide sequence of an organism
- Gene annotations file (CSV) - describes genes and other features of DNA, RNA and protein sequences
- Gene-finding format (GFF) - similar to the gene annotations file, but used mainly in the pre-procecssing steps for RNA-seq
- Gene transfer format (GTF) - similar to GFF in that it's the main file used in pre-processing, different table formatting and conventions
- Gene lengths (CSV) - made from a gene annotations file + Python processing script
  - Columns: gene_id, length

## HPC requirements
This repository assumes work on the [Turing/Ace HPC cluster](https://docs.turing.wpi.edu/), but code can be easily modified for other Linux-run platforms. Refer to the 'Getting Started' page of the Turing/Ace documentation for setup.

Anaconda is essential for downloading packages and creating environments. Follow the [Anaconda documentation](/docs/getting-started/anaconda/install#verify-your-install) to install the Linux version of Anaconda.

Create a new conda environment by running the following:
```
conda env create -f environment.yml 
```

This will install the Linux and Python libraries that are used in this pipeline automatically, instead of having to install each package individually.

If for some reason not all packages were able to be installed within the environment, activate the conda environment and install each package manually. Below is an example:

```
conda install bioconda::samtools
```

When in doubt, search packages on the [Anaconda website](https://anaconda.org/) to see how to install a package.   

[comment]: <> (The following are downloaded from the source and need to be updated in the PATH environment &#40;ie. on the .bashrc script&#41;:)

[comment]: <> (- [Entrez Direct]&#40;https://www.ncbi.nlm.nih.gov/books/NBK179288/&#41;)

[comment]: <> (- [STAR]&#40;https://github.com/alexdobin/STAR&#41;)

# Preparing files
### FASTQ download script
Once you are on the ENA page for the RNA-seq project you want to analyze ([Example RNA-seq project with data deposited on ENA](https://www.ebi.ac.uk/ena/browser/view/srp157937)), scroll to the bottom of the page and click on the "Download All" button. The download file will be in this format (You will need to add the header manually, as it is specifically required for Turing/Slurm bash files):

```
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
wget -nc <URL_TO_SAMPLE_1.fastq.gz>
wget -nc <URL_TO_SAMPLE_2.fastq.gz>
wget -nc <URL_TO_SAMPLE_3.fastq.gz>
wget -nc <URL_TO_SAMPLE_4.fastq.gz>
```

### Create ```local_config.sh```
The ```local_config.sh``` Will contain variables and file paths that will be read by the main pipeline script. It will have the following format:

```
ORGANISM_NAME='Pseudomonas aeruginosa PAO1'
ROOT_PATH="/path/to/main/analysis/folder"
STUDY_ID="PROJECT_ACCESSION_ID"
DOWNLOAD_SCRIPT_NAME="ena_download.sh"
```
- ORGANISM_NAME: Choices (as of now) are `E coli K12`, ```P putida KT2440```, and ```P aeruginosa PAO1```. To include another organism, the code from ```pipeline_trx.sh``` will need to be modified.
- ROOT_PATH: The absolute path to the folder that will hold all files resulting from using this pipeline.
- DOWNLOAD_SCRIPT_NAME: Make sure to change this name to match the script you downloaded from ENA.

### Generate a STAR index (if needed)

In order to quantify gene expression, you will need to know how many sequencing reads are mapped to each gene. And in order to know how many sequencing read fragments are mapped to each gene, you will need to have an indexed genome. Indexing a genome can be thought of as an index at the beginning of a book. If you want to go to a specific chapter, you can find the page number from the index instead of flipping through every page.

There are multiple tools that can be used for genome indexing and alignment. This pipeline will use the STAR aligner, but [BWA](https://github.com/lh3/bwa) is also sufficient for analysis of bacterial organisms.

Before running the full RNA-seq analysis pipeline, generate a STAR index by running the following:

```
bash generate_star_index.sh <GENOME_FASTA_PATH> <GTF_PATH> <OUTPUT_FOLDER_PATH>
```

Maker sure to replace variables in <brackets> with the real, absolute paths to your files or folders.

### Modify `pipeline_trx.sh`

To include another organism for analysis in this pipeline, you will need to modify the setup part of the code.

Add variables that have the paths to the STAR index and GTF of your desired organism and modify the file paths for _P. putida_, _E. coli_, and _P. aeruginosa_ if you plan on using those organisms:

```
# Absolute paths to reference files
index_KT2440=<PATH_TO_STAR_INDEX>
gtf_KT2440=<PATH_TO_GTF>

index_ecoli_k12=<PATH_TO_STAR_INDEX>
gtf_ecoli_k12=<PATH_TO_GTF>

index_paerug_pao1=<PATH_TO_STAR_INDEX>
gtf_paerug_pao1=<PATH_TO_GTF>

index_pstutz=<PATH_TO_STAR_INDEX>
gtf_pstutz=<PATH_TO_GTF>
```
Within the case/esac coding chunk, add a condition for your new organism. An example is annotated below:

```
# Determine STAR index and GTF based on organism
case "$ORGANISM_NAME" in
    "E coli K12")
        STAR_INDEX="$index_ecoli_k12"
        GTF_PATH="$gtf_ecoli_k12"
        ;;
    "P putida KT2440")
        STAR_INDEX="$index_KT2440"
        GTF_PATH="$gtf_KT2440"
        ;;
    "P aeruginosa PAO1")
        STAR_INDEX="$index_paerug_pao1"
        GTF_PATH="$gtf_paerug_pao1"
        ;;
    ##### START OF EXAMPLE - adding new organism #####
    "P stutzeri STRAIN1")
        STAR_INDEX="$gtf_pstutz"
        GTF_PATH="$gtf_pstutz"
        ;;
    ##### END OF EXAMPLE #####
    *)
        echo "Error: Unknown organism '$ORGANISM_NAME'"
        exit 1
        ;;
esac
```

# Running the Script

Once you have generated your STAR index, modified `local_config.sh` and `pipeline_trx.sh`, you are now ready to run the script!

To begin the RNA-seq preprocessing, run the following:

```
bash pipeline_trx.sh
```

Depending on how many FASTQ files you start with and how big each file is, the piepline may take from 2-80+ hours. If you're using an HPC to run this, you can logout and put your computer to sleep while the pipeline is running. If you're running this on your local computer, make sure your computer stays awake for the whole process.

After successfully running `pipeline_trx.sh`, you should now have the following files and folders in your `ROOT_PATH`:

```
ROOT_PATH
  fastq_raw
  fastq_trim
  fastqc
  align_star
  indexed_HQ
  read_counts
  ena_download.sh
```

By the end, you should have a folder called `read_counts` which holds a .txt table. The first column should have the gene IDs for your organism. Columns 2-6 should have information on the gene type (chromosome, plasmid, etc.), gene start, end, strand, and length. The rest of the columns should be the read counts associated with each starting FASTQ sample you downloaded.
