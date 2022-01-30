
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

                     
# count species which have at least 1 copy of the gene

library(dplyr)

gene_counts_all_sp_noNA <- read.table("gene_counts_all_sp_noNA.txt", header = T) %>%
                           mutate(Description = rowSums(.[2:10] > 0)) %>%  # sum number of elements in each row that are > 0
                           select(Description, everything()) %>%
                           filter(., Description > 4 & mus_musculus > 0) %>% 
                           # filter rows where at least 5 sp. have the gene (outgroup included)
                           rename(canfam = can_fam, eirabarbara = eira_barbara, feliscatus = felis_catus, gulogulo = gulo_gulo,
                           marteszibellina = martes_zibellina, miroungaangustirostris = mirounga_angustirostris,
                           musmusculus = mus_musculus, mustelaputoriusfuro = mustela_putorius_furo,
                           odobenusrosmarus = odobenus_rosmarus)

write.table(gene_counts_all_sp_noNA, file = "gene_counts_all_sp_noNA_sp_counted.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)



