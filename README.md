> [!NOTE]
> This repository is still under construction and is being updated week-to-week.

## TODO!!
- Check files under - /home/gcabebe/rnaseq_mining_scripts
- Check/edit environment.yml file

# Summary
This repository provides an executable script for pre-processing individual RNA-seq studies into read counts per gene tables. Codes assume the use of a high power computing (HPC) cluster and Simple Linux Utility for Resource Management (SLURM), but can be modified for running on a local computer.

## Example dataset
The following RNA-seq project comes from a lab at the Technical University of Denmark. In brief, they tested the gene expression of _P. putida_ on non-trivial carbon sources (ie. citrate, ferulic acid, serine, and glucose) using transcriptomics and genome-scale modeling ([D'Arrigo et al 2019](https://doi.org/10.1111/1758-2229.12704)). Their data is publicly available and can be accessed on [NCBI](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP157937&o=acc_s%3Aa) and the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/view/srp157937).

---
# Setup

## Files to Download
For each organism, you will need the files listed below. For [_P. putida_](https://www.pseudomonas.com/strain/show/110), [_P. aeruginosa_](https://www.pseudomonas.com/strain/show?id=107), and [_E. coli_](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/), these files are already provided within this repository for _P. putida_ under the ```data_P_putida``` folder.
- Genome sequence (FASTA) - nucleotide sequence of an organism
- Gene annotations file (CSV) - describes genes and other features of DNA, RNA and protein sequences
- Gene-finding format (GFF) - similar to the gene annotations file, but used mainly in the pre-procecssing steps for RNA-seq
- Gene transfer format (GTF) - similar to GFF in that it's the main file used in pre-processing, different table formatting and conventions
- Gene lengths (CSV) - made from a gene annotations file + Python processing script
  - Columns: gene_id, length

## HPC Requirements
This repository assumes work on the [Turing/Ace HPC cluster](https://docs.turing.wpi.edu/), but code can be easily modified for other Linux-run platforms. Refer to the 'Getting Started' page of the Turing/Ace documentation for setup.

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