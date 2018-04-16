#!/bin/bash

#PBS -N sort_all
#PBS -q batch
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l mem=200gb
#PBS -l walltime=480:00:00
#PBS -M hmcqueary@uga.edu
#PBS -m ae
#PBS -o /lustre1/hcm14449/SC_RNAseq/sort_all.o
#PBS -e /lustre1/hcm14449/SC_RNAseq/sort_all.e

module load samtools/1.3.1
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat

mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/newSort
while read SampleName
do
time samtools sort -m 100G -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/newSort/${SampleName}/${SampleName}.sorted.bam \
-T /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/${SampleName}/accepted_hits \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/${SampleName}/accepted_hits.bam

done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/samples.txt
