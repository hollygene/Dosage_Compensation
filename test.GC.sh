#!/bin/bash
#PBS -N assembly
#PBS -q batch
#PBS -l nodes=2:ppn=1:AMD
#PBS -l walltime=96:00:00
#PBS -l pmem=2gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

basedir="/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly"
cd $basedir
THREADS=4

# put bowtie 2.2.9 executables in $PATH
module load bowtie2/2.2.9

# build forward and backward indices of reference genome for mapping reads with bowtie
bowtie2-build -f genome.fa genome

#unload bowtie2
module unload bowtie2/2.2.9
#load in tophat
module load tophat/2.1.1

#index transcriptome file
tophat -G genes.gtf --transcriptome-index=transcriptome_data/known genome

#run tophat on all samples
#do GC ones first
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC
mkdir "/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test/"

for file in ./*.fastq

do

FBASE=$(basename $file .fastq)
BASE=${FBASE%.fastq}
tophat -p $THREADS -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test/${BASE}_tophat_out \
-i 10 -I 1000 \
--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
./${BASE}.fastq

done

#unload tophat
module unload tophat/2.1.1

#load in cufflinks
module load cufflinks/2.2.1

#then run cufflinks on all GC samples
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test
printf '%s\n' * > output.txt

mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test

#run cufflinks on all samples
while read SampleName
do
mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinks${SampleName}

cufflinks \
-p $THREADS
-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinks${SampleName} \
-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test/${BASE}_tophat_out

done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test/output.txt

mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/quant_test
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/quant_test
#run cuffquant to estimate expression levels to put into cuffdiff
while read SampleName
do
mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/quant_test/cuffquant${SampleName}

cuffquant \
-p $THREADS
-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/quant_test/cuffquant${SampleName} \
-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test/${BASE}_tophat_out/accepted_hits.bam

done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test/output.txt
