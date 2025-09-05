# MAKE SURE TO CHANGE 'gff_file', 'csv_file' AND 'gene_lengths_file'!!!

import csv
import pandas as pd

# CHANGE THESE VARIABLES!!!
gff_file = 'references_Ecoli/Ecoli_ASM584v2.gff'
csv_file = 'references_Ecoli/gene_annotation_E_coli_str_k_12.csv'
gene_lengths_file = 'references_Ecoli/gene_lengths_E_coli_str_k_12.csv'

def parse_attributes(attr_str):
    attrs = {}
    for pair in attr_str.strip().split(';'):
        if '=' in pair:
            key, value = pair.split('=', 1)
            attrs[key] = value
    return attrs

rows = []
fieldnames = set(['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase'])

# First pass: parse and collect all fieldnames
with open(gff_file, 'r') as infile:
    for line in infile:
        if line.startswith('#'):
            continue

        fields = line.strip().split('\t')
        if len(fields) != 9:
            continue

        seqid, source, ftype, start, end, score, strand, phase, attr_str = fields
        base_row = {
            'seqid': seqid,
            'source': source,
            'type': ftype,
            'start': start,
            'end': end,
            'score': score,
            'strand': strand,
            'phase': phase
        }

        attr_dict = parse_attributes(attr_str)
        full_row = {**base_row, **attr_dict}
        rows.append(full_row)
        fieldnames.update(attr_dict.keys())

# Second pass: write to CSV with all fields
fieldnames = list(fieldnames)

with open(csv_file, 'w', newline='') as outfile:
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
        writer.writerow(row)

df = pd.read_csv(csv_file)
df = df[df['product'].notna()]

df.to_csv(csv_file, index=False)

# Create gene lengths annotation file
df['length'] = abs(df['end'] - df['start'])
gene_lengths_df = df[['protein_id',  'length']]
gene_lengths_df.to_csv(gene_lengths_file, index=False)