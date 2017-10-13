#!/bin/bash

#PBS -q batch
#PBS -l nodes=1:ppn=4:AMD
#PBS -l walltime=96:00:00
#PBS -l pmem=20gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae
#PBS -o /lustre1/hcm14449/SC_RNAseq/sort_all.o
#PBS -e /lustre1/hcm14449/SC_RNAseq/sort_all.e

module load tophat/2.1.1

cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat

for file in ./*.fastq

do

FBASE=$(basename $file .fastq)

BASE=${FBASE%.fastq}

tophat -p 4 -o ./${BASE}_tophat_out -i 10 -I 1000 --transcriptome-index=transcriptome_data/known genome ${BASE}.fastq

done
