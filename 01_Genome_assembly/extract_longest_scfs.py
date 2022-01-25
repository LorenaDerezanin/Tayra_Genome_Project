# done in Supernova_asm_tayra conda env
# python 3.7.3
# run as time python ../08_circos_plots/extract_scfs_by_length.py 


from Bio import SeqIO

longest_scaffolds = [] # empty list

for record in SeqIO.parse("tayra_asm2_haplo.1.fasta", "fasta"):
    if len(record.seq) >= 100000 : # e.g. 1 Mbp long
        longest_scaffolds.append(record) # add record to the list

print("Found %i long scaffolds" % len(longest_scaffolds))
# e.g. found 58 long scfs (10 Mbp or longer)

SeqIO.write(longest_scaffolds, "tayra_longest_scfs_100kb.fasta", "fasta")


# create a file with scf IDs and lengths

for record in SeqIO.parse("tayra_longest_scfs_100kb.fasta", "fasta"):
    print("scf_%s %i" % (record.id, len(record.seq)))
   
# SeqIO.write(scf_list, "tayra_scf_list_lengths_100kb.tab", "tab")


# record_dict = SeqIO.to_dict(SeqIO.parse("longest_scfs_10M.fasta", "fasta"))
# print(record_dict["105"]) 



##### create args for easier usage

# shorten the names of scaffolds 
# sed '/^>/ s/ .*//' tayra_longest_scfs_100kb.fasta >> tayra_longest_scfs_100kb_short_names.fa
