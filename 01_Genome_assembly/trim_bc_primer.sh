
# trim 10x barcodes(16bp) and primer sequence (7bp)

READS=/home/derezanin/supernova_tayra/04_tayra_reads

python2.7 process_10xReads.py --bctrim 16 --trim 7 -a -o $READS/eira_bc_trmd -1 $READS/barcoded_reads/tayra_S1_L001_R1_001.fastq.gz -2 $READS/barcoded_reads/tayra_S1_L001_R2_001.fastq.gz

trim_galore --paired -q 30 --length 80 -j 8 --basename eira_bc_qual_trmd -o $READS/bc_qual_trimmed_reads/ $READS/eira_bc_trmd_R1_001.fastq.gz $READS/eira_bc_trmd_R2_001.fastq.gz
