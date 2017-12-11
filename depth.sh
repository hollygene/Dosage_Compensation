#!/bin/bash
#PBS -N depth_test
#PBS -q batch
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=200gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/old/

module load samtools/1.3.1

while read SampleName
do
mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/old/Depth${SampleName}

samtools depth \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/old/${SampleName}/${SampleName}.sorted.bam \
|  awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${SampleName}.txt


done < /home/hcm14449/Github/Dosage_Compensation/bams.txt
