#!/bin/bash

#PBS -N sort_allMA
#PBS -q batch
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l mem=200gb
#PBS -l walltime=480:00:00
#PBS -M hmcqueary@uga.edu
#PBS -m ae
#PBS -o /lustre1/hcm14449/SC_RNAseq/sort_all.o
#PBS -e /lustre1/hcm14449/SC_RNAseq/sort_all.e

module load samtools/1.3.1
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat

while read SampleName
do
time samtools sort -m 200G \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/${SampleName}/${SampleName}.sorted.bam \
   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/${SampleName}/accepted_hits.bam

done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/samples.txt
