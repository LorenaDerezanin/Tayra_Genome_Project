#!/usr/bin/python

import pandas as pd
from textwrap import wrap


# create dictionary out of codon table
# squeeze = T - important for non-nested dict
dna_codon_table = pd.read_csv('DNA_codon_table.csv', header=None, index_col=0, squeeze=True).to_dict()

for record in SeqIO.parse('', "fasta"):
    