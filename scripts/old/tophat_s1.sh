#!/bin/bash
#PBS -N sample1_3
#PBS -q batch
#PBS -l nodes=1:ppn=4:AMD
#PBS -l walltime=96:00:00
#PBS -l pmem=20gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/fastq
module load tophat/2.1.1

time  tophat -p 4 -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/tophat/SC001_raw_tophat_out \
-i 10 -I 1000 \
--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/transcriptome_data/known \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/genome \
./SC001_raw.fastq
