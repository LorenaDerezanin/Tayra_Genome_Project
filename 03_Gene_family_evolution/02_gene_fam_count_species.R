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
