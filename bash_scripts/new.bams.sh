#!/bin/bash

#$ -q rcc-30d

# S.c. Gene Conversion
date

AR=(90 91 92)

workdir = /project/dwhlab/Holly/GeneConv/GeneConv_3/Coverage


for i in "${AR[@]}"
do
module load bedtools/2.25.0
  bedtools bamtobed -i /lustre1/hcm14449/SC_RNAseq/RNA_seq/CovFiles/Sample29.sorted.fixed.marked.realigned.bam > /lustre1/hcm14449/SC_RNAseq/RNA_seq/CovFiles/Sample29.bed
  sort -k 1,1 -k2,2n /lustre1/hcm14449/SC_RNAseq/RNA_seq/CovFiles/Sample29.bed > /lustre1/hcm14449/SC_RNAseq/RNA_seq/CovFiles/Sample29.sorted.bed
  #let's try switching a and b FILES
  bedtools coverage -a /lustre1/hcm14449/SC_RNAseq/RNA_seq/CovFiles/Sc10K.bed -b /lustre1/hcm14449/SC_RNAseq/RNA_seq/CovFiles/Sample29.sorted.bed > /lustre1/hcm14449/SC_RNAseq/RNA_seq/CovFiles/Sample29.Coverage10K.txt
cut -f 1,2,4 /lustre1/hcm14449/SC_RNAseq/RNA_seq/CovFiles/Sample29.Coverage10K.txt > /lustre1/hcm14449/SC_RNAseq/RNA_seq/CovFiles/Sample29.Coverage10K.shrt.txt
#cut -f 1,2,3 /lustre1/hcm14449/SC_RNAseq/RNA_seq/CovFiles/CoverageHead.txt > /lustre1/hcm14449/SC_RNAseq/RNA_seq/CovFiles/CoverageHead2.txt
cat /lustre1/hcm14449/SC_RNAseq/RNA_seq/CovFiles/CoverageHead2.txt /lustre1/hcm14449/SC_RNAseq/RNA_seq/CovFiles/Sample29.Coverage10K.shrt.txt > /lustre1/hcm14449/SC_RNAseq/RNA_seq/CovFiles/Sample29.Coverage10K.head.txt


done
