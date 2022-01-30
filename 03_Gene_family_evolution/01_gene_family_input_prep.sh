
#### GENE FAMILY ANALYSIS ####

# extract ortologs for each species, rm dups, convert other statuses into counts
for f in *_busco4/full_table*
do
  dir_name=$(dirname $f)
  base=$(basename $dir_name "_busco4")
  grep -v '^#' $f | cut -f 1,2  | sed "s/Complete/1/" - | sed "s/Missing/0/" - | \
  sed "s/Fragmented/NA/" - | grep -v 'Duplicated' > ${base}_OGs_extracted.txt
done


# count duplicates and append to previous output

for f in *_busco4/full_table*
do
  dir_name=$(dirname $f)
  base=$(basename $dir_name "_busco4")
  grep -v '^#' $f | cut -f 1,2 | grep 'Duplicated' | sort | cut -f 1 | \
  uniq -c | awk 'BEGIN {FS=" "; OFS="\t";}{print $2,$1;}' > ${base}_OGs_duplicated.txt
done


#  add headers to coulmns, sort
for f in *_OGs_duplicated.txt
do
  base=$(basename $f "_OGs_duplicated.txt")
  echo -e "ID\t${base}" > ${base}_OGs_extracted_srt.txt
  cat $f ${base}_OGs_extracted.txt | sort >> ${base}_OGs_extracted_srt.txt
done


# get IDs of all orthos
cut -f 1 can_fam_OGs_extracted_srt.txt > IDs_all_OGs.txt


# get gene counts for each species
for f in *_OGs_extracted_srt.txt
do
  base=$(basename $f "_OGs_extracted_srt.txt")
  cut -f 2 $f > ${base}_f2.txt
done


# paste back the first column with gene IDs to each sp. file
paste IDs_all_OGs.txt *_f2.txt > gene_counts_all_sp.txt



# remove lines with fragmented genes ('NA')
grep -v 'NA' gene_counts_all_sp.txt > gene_counts_all_sp_noNA.txt


## switch to R, count genes per species

                     



