# Check for premature stop codons in nt sequences of paralogs
import pandas as pd
import os
from textwrap import wrap
from Bio import SeqIO


# create dictionary out of codon table
# squeeze = T - important for non-nested dict
dna_codon_table = pd.read_csv('DNA_codon_table.csv', header=None, index_col=0, squeeze=True).to_dict()

stop_codons = []
# path = # add path to dir ("BUSCO_gene_sets/gene_fam_expansions/expansions_all_sp")
for file_name in os.listdir(path):
    for record in SeqIO.parse(path + file_name, "fasta"):
        s = str(record.seq)
        dna_codons = wrap(s, 3)
        for p, c in enumerate(dna_codons):
            if dna_codon_table.get(c) == "Stop" and p + 1 != len(dna_codons):
                print(f"{file_name}:{record.id}: {p + 1}/{len(dna_codons)}, {c}")

# output:
# file name:record ID: codon position/number of codons, stop codon