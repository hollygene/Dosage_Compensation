#!/bin/bash

#PBS -q batch
#PBS -l nodes=1:ppn=4:AMD
#PBS -l walltime=96:00:00
#PBS -l pmem=20gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/Cufflinks
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/Cufflinks

module load cufflinks/2.2.1

while read SampleName
do
mkdir /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/Cufflinks/cufflinks${SampleName}

time cufflinks -p 4
-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/genes.gtf \
-o /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/Cufflinks/cufflinks${SampleName} \
-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophatgenome.fa \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/tophat${SampleName}/accepted_hits.bam

done < /home/hcm14449/Github/Dosage_Compensation/MAsamples.txt
