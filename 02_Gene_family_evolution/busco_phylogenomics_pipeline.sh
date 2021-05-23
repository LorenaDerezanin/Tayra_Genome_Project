
### BUSCO gene sets for estimating gene and species trees ###


# 1. Get list of BUSCO IDs for protein/nucleotide sequences for each species

for f in *_busco4
do
  base=$(basename $f "_busco4")
  ls *.fna $f/run_mammalia_odb10/busco_sequences/single_copy_busco_sequences/ > ${base}_sc_fna.txt
done

# concatenate all sp. gene lists in one
cat *.txt > 03_merged_gene_set_fna.txt 


# 2. Get single-copy busco IDs shared in all species

# changed grep to 9 species (mus musculus added)
cat 01_merged_gene_set.txt | sort | uniq -c | sort -nr | grep -P " 9 " | cut -f 8 -d " " > 02_overlap_gene_set.txt

# Get species names from subdir names
sed 's/\(.*\)busco4/\1/' sp_list.txt > species_list_clean.txt




# Add species name to sequence ID line, save files in new dir
for sp in $(cat species_list_clean.txt)
do
  mkdir -p ${sp}_busco_species_added_faa
  for f in ${sp}_busco4/run_mammalia_odb10/busco_sequences/single_copy_busco_sequences/*.faa
  do
    base=$(basename $f ".faa")
    sed -e "s/^>\(.*\) .*/>${sp}/" $f > ${sp}_busco_species_added_faa/${base}.faa
  done
done


# add sp. names to nucleotide sequences, save in new dir, use short sequence names

for sp in $(cat species_list_clean.txt)
do
  mkdir -p ${sp}_busco_species_added_fna
  for f in ${sp}_busco4/run_mammalia_odb10/busco_sequences/single_copy_busco_sequences/*.fna
  do
    base=$(basename $f ".fna")
    sed -e "s/^>\(.*\) .*/>${sp}/" $f > ${sp}_busco_species_added_fna/${base}.fna
  done
done




# Extract and concatenate faa sequences from all species for each busco ID

mkdir -p 01_busco4_genes_concatenated

for line in $(cat 00_busco4_genes_overlap/02_overlap_gene_set.txt)
do
  cat busco4/*_busco_species_added_faa/${line} > 01_busco4_genes_concatenated/protein_sequences/${line}
done



# extract and concat fna seq. from all species for each BUSCO ID
for line in $(cat 00_busco4_genes_overlap/04_overlap_gene_set_fna.txt)
do
  cat busco4/*_busco_species_added_fna/${line} > 01_busco4_genes_concatenated/nucleotide_sequences/${line}
done



# Multiple sequence alignment


# ran on Allegro (conda env mafft_trimal)
# mafft v7.471 
# trimAl v. 1.4.1


BUSCOS=/home/derezanin/NO_BACKUP/01_bff_mustela_nigripes/02_busco_runs/01_busco4_genes_overlap
ALN=/home/derezanin/NO_BACKUP/01_bff_mustela_nigripes/02_busco_runs/02_busco4_alignments
TRM=/home/derezanin/NO_BACKUP/01_bff_mustela_nigripes/02_busco_runs/03_busco4_trimmed_alignments

for f in $BUSCOS/*.faa
do
  faa_base=$(basename $f ".faa")
  mafft --thread 20 $f > $ALN/${faa_base}.aln
  alignment="$ALN/${faa_base}.aln"
  aln_base=$(basename $alignment ".aln")
  trimal -in $alignment -out $TRM/${aln_base}.trm -strictplus
done


# 1500 files processed in ~13 min
# whole set (6020 buscos) ~50 min



#### Codon MSA ####


# conda prank
# prank v.170427


# run as array job on Allegro

FNA_BUSCOS=/home/derezanin/NO_BACKUP/01_bff_mustela_nigripes/02_busco_runs/01_busco4_genes_concatenated/nucleotide_sequences
FNA_ALN=/home/derezanin/NO_BACKUP/01_bff_mustela_nigripes/02_busco_runs/02_busco4_alignments/nucleotide_alignments
FNA_TRM=/home/derezanin/NO_BACKUP/01_bff_mustela_nigripes/02_busco_runs/03_busco4_trimmed_alignments/nucleotide_alignments_trm


LIST=$(ls $FNA_BUSCOS/*.fna)
files_arr=($LIST) # transform list in an array
# indexes in files_arr are 0-based

# get array element which is indexed as slurm array task id 
# slurm array task id is the ordinal number of the file


for (( i=$SLURM_ARRAY_TASK_ID; i<6019; i+=100 ))
do
  f=${files_arr[${i}]}
  echo "Processing $i $f" >> $FNA_ALN/prank_array.log
  fna_base=$(basename $f ".fna")
  prank -d=$f -o=$FNA_ALN/${fna_base}_codon.aln -codon -F -verbose
done



# conda mafft_trimal as array
for f in $FNA_ALN/*.aln
do
  fna_aln_base=$(basename $f ".aln.best.fas")
  trimal -in $fna_alignment -out $FNA_TRM/${fna_aln_base}.trm -strictplus
done



### GENE TREES ###


# ran locally on Allegro (conda env raxml), sbatch script was significantly slower with 20 cpu(1GB per cpu) - SSE3 issue?
# raxml v.8.2.12

GENE_TREES=/home/derezanin/NO_BACKUP/01_bff_mustela_nigripes/02_busco_runs/04_busco4_gene_trees_local

for f in $TRM/*.trm
do
  trm_base=$(basename $f ".trm")
  raxmlHPC-PTHREADS-SSE3 -T 10 -s $f -n ${trm_base}.gtree -m PROTGAMMAJTT -p 12345 -w $GENE_TREES
done

# ran locally on 10 cpus in ~1,2 h



# compare gene trees with newick utils to get most common trees

# draw and check one gene tree
nw_display RAxML_bestTree.100879at40674.gtree

# remove branch lengths, preserving only the topology, view the tree
nw_topology RAxML_bestTree.100879at40674.gtree | nw_display -

# reroot the tree with mus musculus as an outgroup
nw_topology RAxML_bestTree.100879at40674.gtree | nw_reroot - mus_musculus_NC_000069.7 | nw_display -

# sort alphabetically
nw_topology RAxML_bestTree.100879at40674.gtree | nw_reroot - mus_musculus_NC_000069.7 | nw_order - | nw_display -


# do it for all trees
cat RAxML_bestTree.* > best_trees.txt


# reroot and sort trees
nw_topology best_trees_renamed.txt | nw_reroot - mus_musculus | nw_order - > best_trees_rerooted_srt.txt

# sort output numerically (n) by the first column (k1)
 sort best_trees_rerooted_srt.txt | uniq -c | sort -nrk1 

# too many trees occurring in low freq due to diff. scaffold IDs, rename all species to short names
# prepare uniq old names mapped to new ones
for sp in $(cat species_list.txt)
do
  nw_labels best_trees.txt | grep "$sp" | sort | uniq | awk '{print $0,"\t",sp}' sp="$sp" > ${sp}_long_name.txt
done

# cat all in 1 map 
cat *_long_name.txt > all_sp_long2short_names.map

# rename species in all trees
nw_rename best_trees.txt all_sp_long2short_names.map > best_trees_renamed.txt

# reroot, sort, and count again 
nw_topology best_trees_renamed.txt | nw_reroot - mus_musculus | nw_order - > best_trees_rerooted_srt.txt
sort best_trees_rerooted_srt.txt | uniq -c | sort -nrk1 

# top 10 trees:
    # 630 
    # 384 
    # 170 
    # 155 
    # 128 
    # 106 
    #  88 
    #  85 
    #  76 
    #  72 

# check first 3 with nw_display
echo '(((can_fam,(((eira_barbara,(gulo_gulo,martes_zibellina)),mustela_putorius_furo),(mirounga_angustirostris,odobenus_rosmarus))),felis_catus),mus_musculus);' \
> topology1.txt
echo '(((can_fam,((((eira_barbara,martes_zibellina),gulo_gulo),mustela_putorius_furo),(mirounga_angustirostris,odobenus_rosmarus))),felis_catus),mus_musculus);' \
> topology2.txt
echo '(((can_fam,((((eira_barbara,gulo_gulo),martes_zibellina),mustela_putorius_furo),(mirounga_angustirostris,odobenus_rosmarus))),felis_catus),mus_musculus);' \
> topology3.txt



### SPECIES TREE ###


# Concatenate busco sequences in busco ID order for each species

# for each .trm file
# save/append basename of each .trm file in a new file  
# nested loop - for each sp. - extract sp. name after ^> and before 2nd "_" and create new file with each sp. name as file name
# locate sequence under seq id line, strip it of \n 
# save/append sequence to previously created sp. file


mkdir -p 05_busco4_merged_trm_alignments

for f in $TRM/*.trm
do
  trm_base=$(basename $f ".trm")
  # keep busco ids to confirm order
  echo $trm_base >> busco_ids.txt
  cat $f | while read line 
  do
    if [[ $line =~ ^\> ]]
    then 
      # [^_]* - any char except the "_"
      species_name=$(echo $line | sed -e 's/^>\([^_]*_[^_]*\)_.*/\1/')   
      # check if file doesn't exist
      if [ ! -f 05_busco4_merged_trm_alignments/${species_name}_merged_buscos.aln ]
      then 
        # create file and write sp name
        echo ">$species_name" > 05_busco4_merged_trm_alignments/${species_name}_merged_buscos.aln
      fi
    else
      echo -n $line >> 05_busco4_merged_trm_alignments/${species_name}_merged_buscos.aln 
    fi
  done
done

# concatenate each sp merged busco aln into one file
for f in 05_busco4_merged_trm_alignments/*_merged_buscos.aln
do
  # add new line before new species name
  (cat $f; echo) >> concatenated_buscos_all_sp.aln
done



# trim concat aln
trimal -in concatenated_buscos_all_sp.aln -out concatenated_buscos_all_sp.trm -strictplus

# remove gaps from trm aln
trimal -in concatenated_buscos_all_sp.trm -out concatenated_buscos_all_sp.trm.nogaps -nogaps

# build species tree
raxmlHPC-PTHREADS-SSE3 -T 12 -s concatenated_buscos_all_sp.trm.nogaps \
-n speciestree1000bs -f a -N 1000 -x 12345 -p 12345 -m PROTGAMMAJTT >& raxml_sp_tree1000bs.log

# reroot and sort sp. tree
nw_reroot RAxML_bipartitions.speciestree1000bs mus_musculus | nw_order - | nw_display -

# count topology ccurrences
nw_reroot RAxML_bootstrap.speciestree1000bs mus_musculus | nw_order - | sort | uniq -c | sort -nrk1


### Ultrametric tree - equal branch lengths ###

nw_reroot RAxML_bestTree.speciestree1000bs mus_musculus | nw_order - > RAxML_bestTree.speciestree1000_reroot_srt


# switch to R

library(ape)
moltree<-read.tree("RAxML_bestTree.speciestree1000_reroot_srt")
calib<-makeChronosCalib(moltree,node="root",age.min=96)
# enter interactive mode and check node numbers 
# makeChronosCalib(moltree,interactive = TRUE)

# pass multiple divergence times
calib_new <- makeChronosCalib(moltree,node = c(10,11),age.min = c(96,54))


timtree<-chronos(moltree,calibration=calib,lambda=1,model="discrete")
timtree<-chronos(moltree,calibration=calib_new,lambda=1,model="discrete")

is.ultrametric(timtree)
write.tree(timtree)
timtree$edge.length<-round(as.numeric(timtree$edge.length), digits=2)
is.ultrametric(timtree)
write.tree(timtree)

# save tree as svg
nw_display -t -u 'MY' -s -w 1000 sp_time_tree2 > my_time_tree.svg


nw_reroot RAxML_bipartitions.speciestree1000bs mus_musculus | nw_order - | nw_display - -s -w 1000 > species_tree_branch_support.svg


# canidae and felidea divergence time and carnivora vs mouse divergence time





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
 

