#!/bin/bash
#PBS -N test_scripts
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=480:00:00
#PBS -l mem=20gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

basedir="/scratch/hcm14449/DC_Feb2019/GC/fastqsGC"
THREADS=4
fastqc_module="FastQC/0.11.8-Java-1.8.0_144"
trimgalore_module="Trim_Galore/0.4.5-foss-2016b"
######### GC LINES ##############
GC_data_dir="/scratch/hcm14449/DC_Feb2019/GC/fastqsGC"
# mkdir /scratch/hcm14449/DC_Feb2019/GC/fastqsGC/fastQC
GC_fastQC_dir= "/scratch/hcm14449/DC_Feb2019/GC/fastqsGC/fastQC"
# ref_genome="/project/dwhlab/Holly/SC_RNAseq/RNA_seq/November_2017_Assembly/SCerRefGenome/genome.fa"
GC_trimmed_dir="/scratch/hcm14449/DC_Feb2019/GC/trim"

#run FASTQC on samples to make sure things look good and I don't need to remove any adapters or anything

# trying to rename files so that I can eventually run a loop with cuffdiff
# rename 'GC' A *.fq

# script is functional
# for file in $data_dir/*.sam
#
# do
#
# FBASE=$(basename $file .sam)
# BASE=${FBASE%.sam}
#
# samtools view -b -t /scratch/hcm14449/TE_MA_Paradoxus/Practice/files/ref_genome/SCerevisiae.RefGenome.fa \
# $data_dir/${BASE}.sam > $data_dir/${BASE}.bam
#
# done

module load ${fastqc_module}

for file in $GC_data_dir/*.fastq

do

FBASE=$(basename $file .fastq)
BASE=${FBASE%.fastq}

fastqc $GC_data_dir/${BASE}.fastq -o $GC_fastQC_dir

done

module unload ${fastqc_module}
