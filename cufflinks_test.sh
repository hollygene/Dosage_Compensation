#!/bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=96:00:00
#PBS -l pmem=5gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/Cufflinks

module load cufflinks/2.2.1

mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/Cufflinks/cufflinksSC007_3_raw_tophat_out

time cufflinks \
-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/genes.gtf \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/Cufflinks/cufflinksSC007_3_raw_tophat_out \
-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophatgenome.fa \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/tophat/SC007_3_raw_tophat_out/accepted_hits.bam
