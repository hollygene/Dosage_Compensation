#!/bin/bash
#PBS -N tophat_allMAand66AGC
#PBS -q batch
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l walltime=96:00:00
#PBS -l pmem=100gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_new/fastq
module load tophat/2.1.1

for file in ./*.fq

do

FBASE=$(basename $file .fq)
BASE=${FBASE%.fq}
tophat -p $THREADS -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/${BASE}_tophat_out \
-i 10 -I 1000 \
--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
./${BASE}.fq

done

tophat -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_GC_66_S12_R1_001_trimmed_tophat_out \
-i 10 -I 1000 \
--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/Holly_GC_66_S12_R1_001_trimmed.fq
