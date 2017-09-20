##script to run cuffnorm on everything

export PATH=/usr/local/cufflinks/latest/bin/:$PATH 
time /usr/local/cufflinks/latest/bin/cuffnorm -o /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/CufflinksOutput2/cuffnorm_out \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/genes.gtf \
-p 5 \
-L X112,X115,X117,X11,X123,X141,X152,X18,X1,X21,X29,X2,X31,X3,X49,X4,X50,X59,X5,X61,X66,X69,X6,X76,X77,X7,X8,X9,GC_Anc,MA_Anc \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_MA_112_S8_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_112B_S28_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_112C_S25_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_MA_115_S9_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_115B_S29_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_115C_S26_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_MA_117_S10_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_117B_S30_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_117C_S27_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly-GC-11-redo_S6_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_11B_S15_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_11C_S12_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_MA_123_S11_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_123B_S31_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_123C_S28_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_MA_141_S12_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_141B_redo_S15_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_141C_S29_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_MA_152_S13_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_152B_S33_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_152C_S30_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly-GC-18-redo_S7_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_18B_redo_S13_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_18C_S13_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_GC_1_redo_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_1B_S6_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_1C_S3_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly-GC-21-redo_S8_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_21B_S17_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_21C_S14_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_MA_29_S6_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_29B_S18_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_29C_redo_S17_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_2A_S1_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_2B_S7_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_2C_S4_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly-GC-31-redo_S9_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_31B_S19_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_31C_S16_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_3A_S2_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_3B_redo_S12_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_3C_S5_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly-GC-49-redo_S10_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_49B_S20_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_49C_S17_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly-GC-4-redo_S3_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_4B_S9_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_4C_S6_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_MA_50_S7_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_50B_S21_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_50C_S18_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_GC_59_S10_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_59B_S22_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_59C_S19_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_5A_S3_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_5B_S10_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_5C_redo_S16_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly-GC-61-redo_S11_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_61B_S23_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_61C_S20_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_GC_66_S12_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly-66B-redo-136232_S14_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_66C_S21_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_GC_69_S1_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_69B_S25_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_69C_S22_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_6A_S4_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_6B_S11_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_6C_S8_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_GC_76_S2_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_76B_S26_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_76C_redo_S18_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_GC_77_S3_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_77B_S27_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_77C_S24_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly-GC-7-redo_S4_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_7B_S12_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_7C_S9_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly-GC-8-redo_S5_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_8B_S13_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_8C_S10_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_9A_S5_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_9B_S14_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_9C_S11_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_GC_Anc_B_S2_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_GC_Anc_C_S32_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_GC_Anc_S4_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_MA_Anc_B_S1_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly-MA-Anc-C-redo_S19_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/Holly_MA_Anc_S5_R1_001_tophat_out/accepted_hits.bam 