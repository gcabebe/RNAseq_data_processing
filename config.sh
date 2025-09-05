RNASEQ_SCRIPTS_PATH="/home/gcabebe/rnaseq_mining_scripts"

# FOR PSEUDOMONAS PUTIDA KT2440
ESEARCH_QUERY="Pseudomonas putida KT2440[ORGN] AND RNA-Seq"
EFETCH_FILENAME="pputida_KT2440_runinfo.csv"
STAR_INDEX="/home/gcabebe/rnaseq/strain_KT2440/index_STAR"
GTF_PATH="/home/gcabebe/rnaseq/strain_KT2440/refseqKT2440/KT2440_110.gtf"
GENE_LENGTHS_REF_PATH="/home/gcabebe/rnaseq/strain_KT2440/tpm_normalize_generef_KT2440.csv"
#SRA_TABLE="$ROOT_PATH/SraRunTable.csv"

# FOR E COLI K-12
# ESEARCH_QUERY="Escherichia coli K-12[ORGN] AND RNA-Seq"
# EFETCH_FILENAME="ecoli_K12_runinfo.csv"
#STAR_INDEX="/home/gcabebe/rnaseq/strain_Ecoli_K12_MG1655/STAR_index"
#GTF_PATH="/home/gcabebe/rnaseq/strain_Ecoli_K12_MG1655/ref/GTF_E_coli_str_k_12.gtf"
#GENE_LENGTHS_REF_PATH="/home/gcabebe/rnaseq/strain_KT2440/gene_lengths_E_coli_str_k_12.csv"

# FOR PSEUDOMONAS AERUGINOSA
# ESEARCH_QUERY="Pseudomonas aeruginosa PAO1[ORGN] AND RNA-Seq"
# EFETCH_FILENAME="paeruginosa_pao1_runinfo.csv"
#STAR_INDEX="/home/gcabebe/rnaseq/strain_Paeruginosa_PAO1/index_STAR"
#GTF_PATH="/home/gcabebe/rnaseq/strain_Paeruginosa_PAO1/ref_genome/GTF_Paeruginosa_PAO1_107.gtf"
#GENE_LENGTHS_REF_PATH="/home/gcabebe/rnaseq/strain_Paeruginosa_PAO1/ref_genome/gene_lengths_Paeruginosa_PAO1.csv"