#!/bin/bash
#PBS -N kallisto
#PBS -q batch
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=100gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

basedir="/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly"
cd $basedir
THREADS=4

module load kallisto/0.42.5
#first need to build a kallisto index
/usr/local/apps/kallisto/0.42.5/bin/kallisto index -i index genome.fa

#test it on one line first

#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC

#/usr/local/apps/kallisto/0.42.5/bin/kallisto quant -i genome -o test.kallisto --single -l 200 -s 20 Holly-66B-redo-136232_S14_R1_001.fastq
