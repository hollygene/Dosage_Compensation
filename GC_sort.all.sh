#!/bin/bash

#PBS -N sort_all
#PBS -q batch
#PBS -l nodes=1:ppn=4:HIGHMEM
#PBS -l mem=200gb
#PBS -l walltime=480:00:00
#PBS -M hmcqueary@uga.edu
#PBS -m ae
#PBS -o /lustre1/hcm14449/SC_RNAseq/sort_all.o
#PBS -e /lustre1/hcm14449/SC_RNAseq/sort_all.e

module load samtools/1.3.1
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_test

while read SampleName
do
time samtools sort -m 100G -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_test/SORT${SampleName}/${SampleName}.sorted.bam \
-T /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_test/${SampleName}/accepted_hits \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_test/${SampleName}/accepted_hits.bam

done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_test/output.txt
