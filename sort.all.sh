#!/bin/bash

#PBS -N sort_all
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=100gb
#PBS -l walltime=480:00:00

module load samtools/1.3.1
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat

while read SampleName
do
  mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/SORT${SampleName}

time samtools sort -m 100G -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/SORT${SampleName}/${SampleName}.sorted.bam \
-T /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/SORT${SampleName}/accepted_hits \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/SORT${SampleName}/accepted_hits.bam

done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/samples.txt
