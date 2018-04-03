#!/bin/bash
#PBS -N tophat_missing
#PBS -q batch
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l walltime=96:00:00
#PBS -l pmem=200gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/
module load tophat/2.1.1

tophat -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_MA_Anc_S5_R1_001_trimmed \
-i 10 -I 1000 \
--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/Holly_MA_Anc_S5_R1_001_trimmed.fq


tophat -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_MA_Anc_B_S1_R1_001_trimmed \
-i 10 -I 1000 \
--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/Holly_MA_Anc_B_S1_R1_001_trimmed.fq

tophat -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_50C_S18_R1_001_trimmed \
-i 10 -I 1000 \
--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/Holly_50C_S18_R1_001_trimmed.fq

tophat -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/SC006_raw_trimmed \
-i 10 -I 1000 \
--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/SC006_raw_trimmed.fq
