#!/bin/bash

#PBS -N j_s_samtools
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=100gb
#PBS -l walltime=480:00:00

while read SampleName
do
  mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/INDEX${SampleName}
cd $PBS_O_WORKDIR
module load samtools/1.3.1
time  samtools index \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/${SampleName}/accepted_hits.bam

done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/samples.txt
