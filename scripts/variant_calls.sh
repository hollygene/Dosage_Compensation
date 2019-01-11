#!/bin/bash
#PBS -N variant_calling
#PBS -q batch
#PBS -l nodes=2:ppn=4:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=200gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

module load samtools/1.3.1
module load bcftools/1.6
##### variant calling pipeline
## to determine if my replicates for GC/MA lines match what they should be matching or not
##if they don't match, find what they DO match to
#need reference fasta file and alignments (.bam files)
#create a link to those in a certain directory
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants
#cp /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/Sc_genome/genome.fa /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/genome.fa
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants

cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/sam_files
#to call variants: simplest approach is to look for positions where mapped reads consistently have a different base than the reference assembly
# this is called consensus approach
#two steps: samtools mpileup: looks for inconsistencies between ref and aligned reads
# bcftools call: interprets them as variants
#create index of reference genome
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp

#module load samtools/1.3.1
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/sam_files

#while read SampleName
#do

#samtools sort -m 200G -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/sam_files/${SampleName}.sorted.sam \
#-T /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/sam_files/${SampleName} \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/sam_files/${SampleName}

#done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/sam_files/samples.txt

## make sam files into bam files
#for file in ./*.sorted.sam

#do

#FBASE=$(basename $file .sam)
#BASE=${FBASE%.sam}

#samtools view -b ${BASE}.sam  -o ${BASE}.bam
#done

#cp /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/genome.fa /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/genome.fa
#samtools faidx /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/genome.fa

mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/temp/bams

for file in ./*.bam
do
  FBASE=$(basename $file .bam)
  BASE=${FBASE%.bam}

#cp /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/sam_files/${BASE}.bam /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/bams/
#run samtools mpileup

#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/variants3

samtools mpileup -uf /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/genome.fa /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/bams/${BASE}.bam > /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/variants3/${BASE}.raw_calls.bcf
#run bcftools call
#bcftools call -v -m /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/variants3/${BASE}.raw_calls.bcf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/variants3/${BASE}.calls.vcf

done


####################################
#query the vcf files
#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/variants
#module load vcftools/0.1.12b
#have to compress them with bgzip first (in vcftools module on Sapelo)
#for file in ./*.vcf
#do
#  FBASE=$(basename $file .vcf)
#  BASE=${FBASE%.vcf}
#bgzip ${BASE}.vcf

#done

#for file in ./*.vcf.gz
#do
#  FBASE=$(basename $file .vcf.gz)
#  BASE=${FBASE%.vcf.gz}
#bcftools query -o ${BASE}.subset.vcf -f '%CHROM  %POS  %REF  %ALT{0}\n'  ${BASE}.vcf.gz
#done

#bcftools query -f '%CHROM  %POS  %REF  %ALT{0}\n' Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out.calls.vcf > Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out.calls.subset.vcf
