
### PREPARATION of 10XGENOMICS LINKED-READS for the SUPERNOVA ASSEMBLY ###

# Supernova(2.1.1)

# Convert BCL to FASTQ #

# samplesheet.csv needs to be provided, create one at 
# https://support.10xgenomics.com/de-novo-assembly/software/pipelines/latest/using/bcl2fastq-direct
# PATH= abs.path to mustelids/tayra/MDC/191108_A00643_0023_AHHCJVDRXX 

supernova mkfastq --run $PATH \
--id tayra --csv samplesheet.csv --qc --ignore-dual-index --localmem=50 \
--jobmode=local -r 12 -p 10 -w 12

# --qc flag not supported for NovaSeq flowcells, but outputs qc for samples/reads

# run mkfastq in conda env with bcl2fastq v2.20 (NovaSeq requires v.2.20)


# TEST ASSEMBLY RUN with RAW DATA #

# WD= abs.path to mustelids/tayra/MDC/tayra/outs/fastq_path/HHCJVDRXX/tayra

### FIRST ASM RUN ###

supernova run --id=tayra --sample=tayra --lanes=1 --fastqs=$WD --localcores=32 --localmem=258 --maxreads=900000000
# start: 15:20 (3.12.19)																			
# end: 15:50 (8.12.19)
# keep fastq files in the Supernova WD, no subdirs

# run mkoutput to create asm fasta files from binary data
supernova mkoutput \
        --style=pseudohap2 \
        --asmdir= # abs. path to tayra/outs/assembly \
        --outprefix=tayra_asm_haplo 
# pseudohap2 calls 2 separate haplotypes
# time: 18 min (8.12.19)

# run BUSCO on both haplotypes - genome completness assessment

# busco doesn't run, install blast 2.2.31, augustus 3.3.3, busco 3.0.2 through conda
# change bashrc path to busco config, source bashrc in conda env, run in conda
# runs now, run it on 2nd haplo
time python /data/fg2/derezanin/envs/Supernova_asm_tayra/bin/run_busco -i tayra_asm_haplo.2.fasta -o tayra_haplo2_busco \
-l /home/derezanin/species_comp/caniformes/mustelids/wolverine/busco/mammalia_odb9 \
 -m geno -c 10 2> busco_tayra2.log

# run QUAST for more genome metrics 
/data/bioinf/quast/quast_py2.7.py tayra_asm_haplo.1.fasta -o tayra1_quast_out --large -t 6 2> quast_tayra1.log  
# /data/bioinf/quast/quast_py2.7.py tayra_asm_haplo.2.fasta -o tayra2_quast_out --large -t 6 2> quast_tayra2.log  

# contig report, run again with -s (scaffolds option)
/data/bioinf/quast/quast_py2.7.py tayra_asm_haplo.1.fasta -s -o tayra1scf_quast_out --large -t 6 2> quast_tayra1scf.log  


### SECOND ASM RUN ###

# run supernova again with 1,3 bil. reads (80x raw cov) to get ~70x effective cov
supernova run --id=tayra2 --sample=tayra --lanes=1 --fastqs=$WD --localcores=26 --localmem=340 --maxreads=1300000000
# start: at 10:40 (10.12.19) 																									     
# crashed on 12.12.19, assign more cores 
supernova run --id=tayra2 --sample=tayra --lanes=1 --fastqs=$WD --localcores=30 --localmem=340 --maxreads=1300000000
# restarted at 12:15 (12.12.19)
# found orphaned local stage, resumed
# end: 19:00 (16.12.19)
# optimal asm

# run mkoutput to create asm fasta files from binary data
supernova mkoutput \
        --style=pseudohap2 \
        --asmdir=/home/derezanin/supernova_tayra/tayra2/outs/assembly \
        --outprefix=tayra_asm2_haplo 
# time: 23 min (18.12.19)

/data/bioinf/quast/quast_py2.7.py tayra_asm2_haplo.1.fasta -o tayra2_quast_out --large -t 2 2> quast_tayra2.log  

### THIRD ASM RUN ###

# start supernova run on allegretto 
# 1.7 bil. reads as input - use extreme cov parameter
supernova run --id=tayra3 --sample=tayra --lanes=1 --fastqs=$WD --localcores=32 --localmem=320 \
--maxreads=1700000000 --accept-extreme-coverage
# time:  ~5 days with some crashes, suboptimal asm


# run BUSCO and QUAST to inspect genome completness and get assembly metrics
time python /data/fg2/derezanin/envs/Supernova_asm_tayra/bin/run_busco -i tayra_asm3.1.fasta -o tayra_asm3.1_busco \
-l /home/derezanin/species_comp/caniformes/mustelids/wolverine/busco/mammalia_odb9 \
 -m geno -c 10 2> busco_tayra_asm2.log

/data/bioinf/quast/quast_py2.7.py tayra_asm3.1.fasta -o tayra3_quast_out --large -t 2 2> quast_tayra3.log  




