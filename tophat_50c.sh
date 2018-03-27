#!/bin/bash
#PBS -N sample1_3
#PBS -q batch
#PBS -l nodes=1:ppn=4:HIGHMEM
#PBS -l walltime=96:00:00
#PBS -l pmem=200gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_new/fastq
module load tophat/2.1.1

time  tophat -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed/tophat/Holly_50C_S18_R1_001_trimmed_tophat_out \
-i 10 -I 1000 \
--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
./Holly_50C_S18_R1_001_trimmed.fq
