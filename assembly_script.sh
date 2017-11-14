#!/bin/bash
#PBS -N assembly_test
#PBS -q batch
#PBS -l nodes=1:ppn=4:AMD
#PBS -l walltime=480:00:00
#PBS -l mem=5gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

basedir="/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly"
cd $basedir
THREADS=4


#run FASTQC on samples to make sure things look good and I don't need to remove any adapters or anything
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC
#mkdir fastqc
#module load java/jdk1.8.0_20 fastqc
#fastqc *.fastq -o fastqc

#module unload java/jdk1.8.0_20 fastqc

#run trimmomatic to trim sequences a little
module load trimmomatic/0.36
mkdir trimmed

#do this in loop
for file in ./*.fastq

do

FBASE=$(basename $file .fastq)
BASE=${FBASE%.fastq}

java -jar /usr/local/apps/trimmomatic/0.36/trimmomatic-0.36.jar SE -threads $THREADS \
-phred33 ./${BASE}.fastq ./trimmed/${BASE}_trim.fastq
done

module unload trimmomatic/0.36

#load in bwa for mapping
#module load bwa/0.7.15

#bwa index genome.fa

#map reads to ref genome using bwa
#bwa mem Holly_2A_S1_R1_001.fastq

#Holly_2B_S7_R1_001.fastq
#Holly_2C_S4_R1_001.fastq

#load in tophat
#module load tophat/2.1.1

#index transcriptome file
#tophat -G genes.gtf --transcriptome-index=transcriptome_data/known genome

#run tophat on all samples
#do GC ones first
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC
#mkdir "/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test/"

#for file in ./*.fastq

#do

#FBASE=$(basename $file .fastq)
#BASE=${FBASE%.fastq}
#tophat -p $THREADS -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test/${BASE}_tophat_out \
#-i 10 -I 1000 \
#--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
#./${BASE}.fastq

#done

#unload tophat
#module unload tophat/2.1.1

#load in cufflinks
#module load cufflinks/2.2.1

#then run cufflinks on all GC samples
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test
#printf '%s\n' * > output.txt

#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test

#run cufflinks on all samples
#while read SampleName
#do
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinks${SampleName}

#cufflinks \
#-p $THREADS
#-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinks${SampleName} \
#-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test/${BASE}_tophat_out

#done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test/output.txt

#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/quant_test
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/quant_test
#run cuffquant to estimate expression levels to put into cuffdiff

#while read SampleName
#do
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/quant_test/cuffquant${SampleName}

#cuffquant \
#-p $THREADS
#-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/quant_test/cuffquant${SampleName} \
#-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test/${BASE}_tophat_out/accepted_hits.bam

#done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test/output.txt

#module unload cufflinks/2.2.1

#####################################################################################################################
#then MA new
#load in tophat
#module load tophat/2.1.1
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new

#for file in ./*.fastq

#do

#FBASE=$(basename $file .fastq)
#BASE=${FBASE%.fastq}

#tophat -p $THREADS -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/tophat/${BASE}_tophat_out -i 10 -I 1000 \
#--transcriptome-index=transcriptome_data/known \
#genome \
#./${BASE}.fastq

#done
#unload tophat
#module unload tophat/2.1.1

#load in cufflinks
#module load cufflinks/2.2.1
#then run cufflinks on all new MA samples
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/tophat
#printf '%s\n' * > output.txt

#run cufflinks on all samples
#while read SampleName
#do
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/Cufflinks/cufflinks${SampleName}
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/Cufflinks
#cufflinks \
#-p $THREADS
#-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/Cufflinks/cufflinks${SampleName} \
#-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/tophat/${BASE}_tophat_out

#done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/tophat/output.txt

#run cuffquant to estimate expression levels to put into cuffdiff
#while read SampleName
#do
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/quant/cuffquant${SampleName}
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/quant
#cuffquant \
#-p $THREADS
#-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/quant/cuffquant${SampleName} \
#-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/tophat/${BASE}_tophat_out/accepted_hits.bam

#done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/tophat/output.txt

#unload cufflinks
#module unload cufflinks/2.2.1

#####################################################################################################################
#then MA old
#load in tophat
#module load tophat/2.1.1
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old

#for file in ./*.fastq

#do

#FBASE=$(basename $file .fastq)
#BASE=${FBASE%.fastq}

#tophat -p $THREADS -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/tophat/${BASE}_tophat_out -i 10 -I 1000 \
#--transcriptome-index=transcriptome_data/known \
#genome \
#./${BASE}.fastq

#done
#unload tophat
#module unload tophat/2.1.1

#load cufflinks
#module load cufflinks/2.2.1

#then run cufflinks on all new MA samples
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/tophat
#printf '%s\n' * > output.txt

#run cufflinks on all samples
#while read SampleName
#do
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/Cufflinks/cufflinks${SampleName}
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/Cufflinks
#cufflinks \
#-p $THREADS
#-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/Cufflinks/cufflinks${SampleName} \
#-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/tophat/${BASE}_tophat_out

#done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/tophat/output.txt

#run cuffquant to estimate expression levels to put into cuffdiff
#while read SampleName
#do
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/quant/cuffquant${SampleName}
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/quant
#cuffquant \
#-p $THREADS
#-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/quant/cuffquant${SampleName} \
#-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/tophat/${BASE}_tophat_out/accepted_hits.bam

#done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/tophat/output.txt

#run cuffdiff to find differential expression between each line and its ancestor
#cuffdiff -p $THREADS --library-type fr-unstranded -o FNR_cuffdiff --labels wild-type,dFNR MG1655.ref.gtf SRR5344681_cuffquant/abundances.cxb,SRR5344682_cuffquant/abundances.cxb SRR5344683_cuffquant/abundances.cxb,SRR5344684_cuffquant/abundances.cxb

#unload cufflinks
#module unload cufflinks/2.2.1
