#!/bin/bash

#PBS -N sort_112
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=20gb
#PBS -l walltime=480:00:00
#PBS -M hmcqueary@uga.edu
#PBS -m ae
#PBS -o /lustre1/hcm14449/SC_RNAseq/sort_112.o
#PBS -e /lustre1/hcm14449/SC_RNAseq/sort_112.e

module load samtools/1.3.1
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat

time samtools sort -m 100G -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/SORTHolly_112B_S28_R1_001_tophat_out/Holly_112B_S28_R1_001_tophat_out.sorted.bam \
-T /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/Holly_112B_S28_R1_001_tophat_out/accepted_hits \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/Holly_112B_S28_R1_001_tophat_out/accepted_hits.bam
