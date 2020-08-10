#!/bin/bash

#PBS -q batch
#PBS -l nodes=2:ppn=3:AMD
#PBS -l walltime=96:00:00
#PBS -l pmem=10gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/Cufflinks
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/Cufflinks

module load cufflinks/2.2.1

while read SampleName
do
mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/Cufflinks/cufflinks${SampleName}

time cufflinks \
-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/genes.gtf \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/Cufflinks/cufflinks${SampleName} \
-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/genome.fa \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/tophat/${SampleName}/accepted_hits.bam

done < /home/hcm14449/Github/Dosage_Compensation/MAsamples.txt
