#!/bin/bash

#PBS -N sort_all
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=100gb
#PBS -l walltime=480:00:00

while read SampleName
do
  mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/SORT${SampleName}
cd $PBS_O_WORKDIR
module load samtools/1.3.1
time  samtools sort \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/${SampleName}/accepted_hits.bam

done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/samples.txt
