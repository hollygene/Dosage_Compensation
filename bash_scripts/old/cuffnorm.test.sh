#!/bin/bash
#PBS -N cuffnormtest
#PBS -q batch
#PBS -l nodes=2:ppn=3:AMD
#PBS -l walltime=96:00:00
#PBS -l pmem=10gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old

module load cufflinks/2.2.1

cuffnorm -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/Cuffnorm_out_test \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/genes.gtf \
-p 3 \
-L SCA,SC119 \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/tophat/SCA_1_raw_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/tophat/SCA_3_raw_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/tophat/SCA_2_raw_tophat_out/accepted_hits.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/tophat/SC119_1_raw_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/tophat/SC119_2_raw_tophat_out/accepted_hits.bam,/lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/tophat/SC119_3_raw_tophat_out/accepted_hits.bam
