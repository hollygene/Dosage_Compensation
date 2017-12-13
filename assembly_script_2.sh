#!/bin/bash
#PBS -N assembly_test_2
#PBS -q batch
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=200gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

basedir="/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly"
cd $basedir
THREADS=4


#run FASTQC on samples to make sure things look good and I don't need to remove any adapters or anything
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC
#trying to rename files so that I can eventually run a loop with cuffdiff
#rename 'GC' A *.fq
#mkdir fastqc
#module load java/jdk1.8.0_20 fastqc
#fastqc *.fastq -o fastqc

#module unload java/jdk1.8.0_20 fastqc

#trim sequences a little
#trying trimgalore since trimmomatic didnt seem to work
#module load trimgalore/0.4.4

#for file in ./*.fastq

#do

#FBASE=$(basename $file .fastq)
#BASE=${FBASE%.fastq}

#trim_galore --phred33 -q 20 -o trimmed ${BASE}.fastq

#done

#module unload trimgalore/0.4.4

#map to reference using bowtie2
#module load bowtie2/2.2.9
#build index for ref genome
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/
#bowtie2-build genome.fa genome
#bowtie2-inspect -s genome


#load in tophat
#module load tophat/2.1.1

#index transcriptome file
#tophat -G genes.gtf --transcriptome-index=transcriptome_data/known genome

#run tophat on all samples
#do GC ones first
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed
#mkdir "/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_test/"

#for file in ./*.fq

#do

#FBASE=$(basename $file .fq)
#BASE=${FBASE%.fq}
#tophat -p $THREADS -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_test/${BASE}_tophat_out \
#-i 10 -I 1000 \
#--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
#./${BASE}.fq

#done

#unload tophat
#module unload tophat/2.1.1

#load in cufflinks
#module load cufflinks/2.2.1

#then run cufflinks on all GC samples
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_test

#ls >> output.txt
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks

#run cufflinks on all samples
#while read SampleName
#do
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks/cufflinks${SampleName}

#cufflinks \
#-p $THREADS \
#-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks/cufflinks${SampleName} \
#-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_test/${SampleName}/accepted_hits.bam

#done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_test/output.txt
#module unload cufflinks/2.2.1

#module load cufflinks/2.2.1

#need to run cuffmerge for later input into edgeR
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks

#cuffmerge \
#-s /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
#-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
#-p $THREADS \
#/home/hcm14449/Github/Dosage_Compensation/assembly_GTF_list.txt


#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/quant_test
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/quant_test
#run cuffquant to estimate expression levels to put into cuffdiff

#while read SampleName
#do
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/quant_test/cuffquant${SampleName}

#cuffquant \
#-p $THREADS \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/quant_test/cuffquant${SampleName} \
#-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/tophat_test/${BASE}_tophat_out/accepted_hits.bam

#done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_test/output.txt

#module unload cufflinks/2.2.1

#####################################################################################################################
#then MA new

#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new
#trying to rename files so that I can eventually run a loop with cuffdiff
#rename 'GC' A *.fq
#these work
#mkdir fastqc
#module load java/jdk1.8.0_20 fastqc
#fastqc *.fastq -o fastqc

#module unload java/jdk1.8.0_20 fastqc

#trim sequences a little
#trying trimgalore since trimmomatic didnt seem to work
#module load trimgalore/0.4.4
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new

#for file in ./*.fastq

#do

#FBASE=$(basename $file .fastq)
#BASE=${FBASE%.fastq}

#trim_galore --phred33 -q 20 -o trimmed ${BASE}.fastq

#done

#module unload trimgalore/0.4.4

#module load tophat/2.1.1
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed
#mkdir "/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/"

#for file in ./*.fq

#do

#FBASE=$(basename $file .fq)
#BASE=${FBASE%.fq}
#tophat -p $THREADS -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/${BASE}_tophat_out \
#-i 10 -I 1000 \
#--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
#./${BASE}.fq

#done

#unload tophat
#module unload tophat/2.1.1

#not sure if this works yet
module load cufflinks/2.2.1

cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat
#ls >> output.txt
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/Cufflinks

#run cufflinks on all samples
while read SampleName
do
mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/Cufflinks/cufflinks${SampleName}

cufflinks \
-p $THREADS \
-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/Cufflinks/cufflinks${SampleName} \
-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/${SampleName}/accepted_hits.bam

done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/output.txt

module unload cufflinks/2.2.1

#run cuffmerge
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/Cufflinks

#cuffmerge \
#-s /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
#-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
#-p $THREADS \
#/home/hcm14449/Github/Dosage_Compensation/assembly_GTF_list.txt


#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/quant_test
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/quant_test


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


#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old

#these codes work

#mkdir fastqc
#module load java/jdk1.8.0_20 fastqc
#fastqc *.fastq -o fastqc

#module unload java/jdk1.8.0_20 fastqc

#trim sequences a little
#trying trimgalore since trimmomatic didnt seem to work
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed

#module load trimgalore/0.4.4

#for file in ./*.fastq

#do

#FBASE=$(basename $file .fastq)
#BASE=${FBASE%.fastq}

#trim_galore --phred33 -q 20 -o trimmed ${BASE}.fastq

#done

#module unload trimgalore/0.4.4


#module load tophat/2.1.1
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed
#mkdir "/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed/tophat/"

#for file in ./*.fq

#do

#FBASE=$(basename $file .fq)
#BASE=${FBASE%.fq}
#tophat -p $THREADS -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed/tophat/${BASE}_tophat_out \
#-i 10 -I 1000 \
#--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
#./${BASE}.fq

#done

#unload tophat
#module unload tophat/2.1.1

#not sure if this works yet

#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed/tophat
#ls >> output.txt
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/Cufflinks

#module load cufflinks/2.2.1

#while read SampleName
#do
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/Cufflinks/cufflinks${SampleName}

#cufflinks \
#-p $THREADS \
#-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/Cufflinks/cufflinks${SampleName} \
#-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed/tophat/${SampleName}/accepted_hits.bam

#done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed/tophat/output.txt

#module unload cufflinks/2.2.1


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
