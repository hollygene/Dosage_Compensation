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
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_rna.fna.gz

#/usr/local/apps/kallisto/0.42.5/bin/kallisto index -i transcripts.idx GCF_000146045.2_R64_rna.fna.gz

#test it on one line first

cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC

/usr/local/apps/kallisto/0.42.5/bin/kallisto quant -i genome -o test.kallisto --single -l 200 -s 20 Holly-66B-redo-136232_S14_R1_001.fastq
