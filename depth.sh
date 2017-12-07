#!/bin/bash
#PBS -N depth_test
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=480:00:00
#PBS -l mem=10gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/old/

module load samtools/1.3.1

samtools depth  /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/old/Holly-66B-redo-136232_S14_R1_001_tophat_out.sorted.bam |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > 66B_depth.txt
