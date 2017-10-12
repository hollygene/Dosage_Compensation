#!/bin/bash 
#PBS -N test.sample1
#PBS -q batch
#PBS -l nodes=1:ppn=4:AMD
#PBS -l walltime=96:00:00
#PBS -l pmem=20gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae
#PBS -o /lustre1/hcm14449/SC_RNAseq/sort_all.o
#PBS -e /lustre1/hcm14449/SC_RNAseq/sort_all.e

cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat
module load tophat/2.1.1

time  tophat -p 4 -o ./SC001_1_raw.fastq_tophat_out -i 10 -I 1000 \
--./transcriptome-index=transcriptome_data/known genome SC001_1_raw.fastq
