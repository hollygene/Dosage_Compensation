#!/bin/bash
#PBS -N variant_calling
#PBS -q batch
#PBS -l nodes=2:ppn=4:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=200gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

module load samtools/1.6
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

## make sam files into bam files
for file in ./*.sam

do

FBASE=$(basename $file .sam)
BASE=${FBASE%.sam}

samtools view -b ${BASE}.sam  -o ${BASE}.bam
done

cp /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/genome.fa /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/genome.fa
samtools faidx /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/genome.fa

for file in ./*.bam

do
  FBASE=$(basename $file .bam)
  BASE=${FBASE%.bam}
cp /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/sam_files/${BASE}.bam /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/${BASE}.bam
#run samtools mpileup
mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/variants2
samtools mpileup -uf /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/genome.fa /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/${BASE}.bam > /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/variants2/${BASE}.raw_calls.bcf
#run bcftools call
bcftools call -v -m /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/variants2/${BASE}.raw_calls.bcf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/variants/temp/variants2/${BASE}.calls.vcf
done
