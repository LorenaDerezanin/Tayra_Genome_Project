
import pandas as pd
import os
from textwrap import wrap
from Bio import SeqIO


# create dictionary out of codon table
# squeeze = T - important for non-nested dict
dna_codon_table = pd.read_csv('DNA_codon_table.csv', header=None, index_col=0, squeeze=True).to_dict()

s = []
path = '/run/user/1000/gvfs/sftp:host=62.141.164.6,port=22220,user=derezanin/home/derezanin/dorina_dir/BUSCO_gene_sets/gene_fam_expansions/expansions_all_sp/'
for i in os.listdir(path):
    for record in SeqIO.parse(path+i, "fasta"):
        s = str(record.seq)
        dna_codons = wrap(s, 3)
        
