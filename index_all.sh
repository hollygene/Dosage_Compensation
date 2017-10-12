#!/bin/bash

#PBS -N index_all
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=100gb
#PBS -l walltime=480:00:00
#PBS -M hmcqueary@uga.edu
#PBS -m ae
#PBS -o /lustre1/hcm14449/SC_RNAseq/index_all.o
#PBS -e /lustre1/hcm14449/SC_RNAseq/index_all.e

module load samtools/1.3.1
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat

while read SampleName
do
time  samtools index \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/SORT${SampleName}/${SampleName}.sorted.bam

done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/samples.txt
