#!/bin/bash
#PBS -N assembly
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=96:00:00
#PBS -l pmem=5gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

basedir="/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly"
cd $basedir
THREADS=4

# put bowtie 2.2.9 executables in $PATH
module load bowtie2/2.2.9

# build forward and backward indices of reference genome for mapping reads with bowtie
bowtie2-build -f genome.fa genome

#unload bowtie2
module unload bowtie2/2.2.9

#load in tophat
module load tophat/2.1.1

#index transcriptome file
tophat -G genes.gtf --transcriptome-index=transcriptome_data/known genome

#run tophat on all samples
#do GC ones first
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC
mkdir "/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test/"

for file in ./*.fastq

do

FBASE=$(basename $file .fastq)
BASE=${FBASE%.fastq}
tophat -p $THREADS -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test/${BASE}_tophat_out \
-i 10 -I 1000 \
--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
./${BASE}.fastq

done
