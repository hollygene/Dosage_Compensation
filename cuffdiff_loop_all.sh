#!/bin/bash

#PBS -q batch
#PBS -l nodes=2:ppn=3:HIGHMEM
#PBS -l walltime=96:00:00
#PBS -l pmem=100gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae


#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff

module load Cufflinks/2.2.1-foss-2016b

#time cuffdiff \
#-p 3 \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample11 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_11B_S15_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_11C_S12_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_21B_S17_R1_001_trimmed_tophat_out/accepted_hits.bam \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

#time cuffdiff \
#-p 3 \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample18 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-GC-18-redo_S7_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_18B_redo_S13_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_18C_S13_R1_001_trimmed_tophat_out/accepted_hits.bam \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample21 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-GC-21-redo_S8_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_21C_S14_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample31 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-GC-31-redo_S9_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_31B_S19_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_31C_S16_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample4 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-GC-4-redo_S3_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_4B_S9_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_4C_S6_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample49 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-GC-49-redo_S10_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_49B_S20_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_49C_S17_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample61 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-GC-61-redo_S11_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_61B_S23_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_61C_S20_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample7 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-GC-7-redo_S4_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_7B_S12_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_7C_S9_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample8 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-GC-8-redo_S5_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_8B_S13_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_8C_S10_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample1 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_1B_S6_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_1C_S3_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_1_redo_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample2 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_2A_S1_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_2B_S7_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_2C_S4_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample3 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_3A_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_3B_redo_S12_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_3C_S5_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample59 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_59B_S22_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_59C_S19_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_59_S10_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample5 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_5A_S3_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_5B_S10_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_5C_redo_S16_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample69 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_69B_S25_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_69C_S22_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_69_S1_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample6 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_6A_S4_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_6B_S11_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_6C_S8_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample76 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_76B_S26_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_76C_redo_S18_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_76_S2_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample77 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_77B_S27_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_77C_S24_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_77_S3_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam

time cuffdiff \
-p 3 \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cuffdiff/Sample9 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/noMitoSacCer.gtf \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_9A_S5_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_9B_S14_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_9C_S11_R1_001_trimmed_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/accepted_hits.bam
