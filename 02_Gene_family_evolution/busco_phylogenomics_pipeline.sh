
#### GENE FAMILY ANALYSIS ####

# extract ortologs for each sp., remove dups, convert other statuses into counts
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


## switch to R

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
                     



## CAFE - run on Allegro, conda cafe


# estimate one or more birth-death (λ) parameters for the given ultrametric sp. tree and gene family counts
# λ parameter describes the probability that a gene will be gained or lost


# limit memory usage and CPUs, super greedy tool
# no underscores in species names - cafe crashes, fix names


# set up config file for cafe - run_cafe_params.txt:

#!shell
date

# load gene counts
load -i gene_counts_all_sp_noNA_sp_counted.txt -p 0.01 -t 10 -l log.txt

# load ultrametric tree
tree (((canfam:40.75,(((eirabarbara:8.19,(gulogulo:7.11,marteszibellina:7.11):1.09):3.97,mustelaputoriusfuro:12.16):19.66,(miroungaangustirostris:16.65,odobenusrosmarus:16.65):15.18):8.92):13.25,feliscatus:54):42,musmusculus:96)

# david's tree
tree ((((odobenusrosmarus:1.0,miroungaangustirostris:1.0):3.0,(mustelaputoriusfuro:3.0,(eirabarbara:2.0,(marteszibellina:1.0,gulogulo:1.0):1.0):1.0):1.0):1.0,canfam:5.0):1.0,feliscatus:6.0)

# search for 1-parameter model
lambda -s -t (((1,(((1,(1,1)1)1,1)1,(1,1)1)1)1,1)1,1)

# david's tree
lambda -s -t ((((1,1)1,(1,(1,(1,1)1)1)1)1,1)1,1)

report report_run

date


#extract gene IDs from cafe report
sed -E "s/([0-9]+at[0-9]+)([^,]+),?/\1\t\2\n/g"


# warning: not an ultrametric tree, but report created 

# get results, sort by p-value
grep 'at40674' report_run.cafe | sort -nk3



# switch to R

# filter CAFE output

cafe_report <- read.table("report_run.cafe", header = F, sep = "\t", skip = 11) %>%  # skip header and info lines
               select(-V5) %>%
               rename(ID = V1, Newick = V2, Family_wide_p_value = V3, Viterbi_p_value = V4) %>%
               filter(., Family_wide_p_value < 0.01) %>%   # keep results with p-value below cut-off
               arrange(., Family_wide_p_value)  # sort in ascending order 


# summary of the cafe report - all gains/losses

python cafe_scripts/cafetutorial_report_analysis.py -i report_run2.cafe -r 0 -o carnivorans_all_g_fams

# split per mustelid species/node



# plot expansions/contractions on the tree

ggtree
treeio


python3.7 cafe_scripts/CAFE_fig/CAFE_fig.py report_run2.cafe -pf 0.01 --dump test_tree/ -g .svg 



# convert phy to fasta:
  # insert ">" in front of sp. name
  # remove 1st row
  # rename from phy to fa 

for f in *.phy
do
  base=$(basename $f ".phy")
  awk '/_/{gsub (/^/,">")}1' $f | \
  sed '1d' > ${base}.fa
done


# run t-coffee evaluation on wolverine alignments
for f in *.fa
do
  base=$(basename $f ".fa")
  t_coffee -infile $f -evaluate  -outfile ${base}.txt \
  -output=score_ascii,aln,score_html 2> ${base}.log
done


# grep sequence lengths in aln
for f in *.log
do
  base=$(basename $f ".log")
  grep "Length" $f > ${base}_length.txt
done


### Cafe Error Correction ###

# prepare cafe_config.sh:

#!/home/derezanin/miniconda3/envs/cafe/bin/cafe
load -i gene_counts_all_sp_noNA_no_mouse.txt -t 4 -l cafferror.txt
tree ((((odobenusrosmarus:1.0,miroungaangustirostris:1.0):3.0,(mustelaputoriusfuro:3.0,(eirabarbara:2.0,(marteszibellina:1.0,gulogulo:1.0):1.0):1.0):1.0):1.0,canfam:5.0):1.0,feliscatus:6.0)
lambda -s -t ((((1,1)1,(1,(1,(1,1)1)1)1)1,1)1,1)
report caferror_1/cafe_final_report


# run error model estimation
python caferror.py -i cafe_config.sh -v 0 -f 0



# run cafe with the best error model specified

#!cafe
load -i gene_counts_all_sp_noNA_no_mouse.txt -t 4 -l cafferror_model2.txt
tree ((((odobenusrosmarus:1.0,miroungaangustirostris:1.0):3.0,(mustelaputoriusfuro:3.0,(eirabarbara:2.0,(marteszibellina:1.0,gulogulo:1.0):1.0):1.0):1.0):1.0,canfam:5.0):1.0,feliscatus:6.0)
errormodel -model cafe_scripts/caferror_1/cafe_errormodel_0.03173828125.txt -all
lambda -s -t ((((1,1)1,(1,(1,(1,1)1)1)1)1,1)1,1)
report report_run2_caferror_model


# summary of report
python cafe_scripts/cafetutorial_report_analysis.py -i report_run2_caferror_model.cafe -r 0 -o carnivorans_all_g_fams_error_model



# run cafe with separate birth and death rate estimates for test

#!cafe
#version
#date
load -i gene_counts_all_sp_noNA_no_mouse.txt -t 4 -l caffe_sep_estimates.txt
tree ((((odobenusrosmarus:1.0,miroungaangustirostris:1.0):3.0,(mustelaputoriusfuro:3.0,(eirabarbara:2.0,(marteszibellina:1.0,gulogulo:1.0):1.0):1.0):1.0):1.0,canfam:5.0):1.0,feliscatus:6.0)
errormodel -model cafe_scripts/caferror_1/cafe_errormodel_0.03173828125.txt -all
lambdamu -s -t ((((1,1)1,(1,(1,(1,1)1)1)1)1,1)1,1)
report report_run3_err_sep_estimates
 

