#!/bin/bash
#PBS -N cuffnorm_all
#PBS -q batch
#PBS -l nodes=2:ppn=3:AMD
#PBS -l walltime=96:00:00
#PBS -l pmem=10gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

cd /scratch/hcm14449/DC_MAoldRedoJune2019/tophat

module load Cufflinks/2.2.1-foss-2016b

cuffnorm -o /scratch/hcm14449/DC_MAoldRedoJune2019/tophat/Cuffnorm_out \
/scratch/hcm14449/DC_MAoldRedoJune2019/genes.gtf \
-p 3 \
-L SCA1,SC001,SC002,SC003,SC004,SC005,SC006,SC007,SC008,SC009,SC011,SC015,SC028,SC088,SC108,SC119 \
/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SCA_1_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SCA_3_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SCA_2_raw_tophat_out/accepted_hits.bam \
/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC001_1_raw_fastq_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC001_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC001_2_raw_tophat_out/accepted_hits.bam \
/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC002_1_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC002_2_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC002_3_raw_tophat_out/accepted_hits.bam \
/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC003_1_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC003_2_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC003_3_raw_tophat_out/accepted_hits.bam \
/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC004_1_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC004_2_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC004_3_raw_tophat_out/accepted_hits.bam \
/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC005_1_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC005_2_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC005_3_raw_tophat_out/accepted_hits.bam \
/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC006_1_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC006_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC006_3_raw_tophat_out/accepted_hits.bam \
/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC007_1_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC007_2_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC007_3_raw_tophat_out/accepted_hits.bam \
/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC008_1_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC008_2_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC008_raw_tophat_out/accepted_hits.bam \
/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC009_1_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC009_2_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC009_3_raw_tophat_out/accepted_hits.bam \
/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC011_1_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC011_2_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC011_3_raw_tophat_out/accepted_hits.bam \
/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC015_1_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC015_2_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC015_3_raw_tophat_out/accepted_hits.bam \
/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC028_1_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC028_2_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC028_3_raw_tophat_out/accepted_hits.bam \
/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC088_1_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC088_2_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC088_3_raw_tophat_out/accepted_hits.bam \
/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC108_1_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC108_2_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC108_3_raw_tophat_out/accepted_hits.bam \
/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC119_1_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC119_2_raw_tophat_out/accepted_hits.bam,/scratch/hcm14449/DC_MAoldRedoJune2019/tophat/SC119_3_raw_tophat_out/accepted_hits.bam
