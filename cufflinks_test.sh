#!/bin/bash
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=96:00:00
#PBS -l pmem=5gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/Cufflinks

module load cufflinks/2.2.0

mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/Cufflinks/cufflinksSC001_1_raw.fastq_tophat_out

time cufflinks
-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/genes.gtf \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/Cufflinks/cufflinksSC001_1_raw.fastq_tophat_out \
-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophatgenome.fa \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/MA_old/tophatSC001_1_raw.fastq_tophat_out/accepted_hits.bam
