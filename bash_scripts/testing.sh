module load tophat/2.1.1
tophat -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-66B-redo-136232_S14_R1_001_tophat_out \
-i 10 -I 1000 \
--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/Holly-66B-redo-136232_S14_R1_001_trimmed.fq

module load samtools/1.3.1
time samtools sort \
-n \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-66B-redo-136232_S14_R1_001_tophat_out/Holly-66B-redo-136232_S14_R1_001_tophat_out.sorted.bam \
   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-66B-redo-136232_S14_R1_001_tophat_out/accepted_hits.bam

module load python/2.7.8
module load htseq/0.6.1p1

#htseq-count -r name -f bam -s no -i gene_id -t transcript \
#   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out.sorted.bam \
#   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_S4_R1_001_trimmed_tophat_out.txt

#tophat -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out \
#-i 10 -I 1000 \
#--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/Holly_GC_Anc_C_S32_R1_001_trimmed.fq
#time samtools sort \
#-n \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out.sorted.bam \
#   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam

#htseq-count -r name -f bam -r name -s no -i gene_id -t transcript \
#   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out.sorted.bam \
#   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out.txt

#time samtools sort \
#   -n \
#   -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out.sorted.bam \
#      /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_C_S32_R1_001_trimmed_tophat_out/accepted_hits.bam

htseq-count -r name -f bam -r name -s no -i gene_id -t transcript \
      /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-66B-redo-136232_S14_R1_001_tophat_out/Holly-66B-redo-136232_S14_R1_001_tophat_out.sorted.bam \
      /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-66B-redo-136232_S14_R1_001_tophat_out.txt



htseq-count -r name -f bam -s no -i gene_id -t exon \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-66B-redo-136232_S14_R1_001_tophat_out/Holly-66B-redo-136232_S14_R1_001_tophat_out.sorted.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/HTseq_update/Holly-66B-redo-136232_S14_R1_001_trimmed_tophat_out.stranded.txt


#module load samtools/1.3.1
#time samtools sort \
#-n \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-GC-18-redo_S7_R1_001_trimmed_tophat_out/Holly-GC-18A.sorted.bam \
#   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-GC-18-redo_S7_R1_001_trimmed_tophat_out/accepted_hits.bam
#   time samtools sort \
#   -n \
#   -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_18B_redo_S13_R1_001_trimmed_tophat_out/Holly_18B.sorted.bam \
#      /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_18B_redo_S13_R1_001_trimmed_tophat_out/accepted_hits.bam
#      time samtools sort \
#      -n \
#      -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_18C_S13_R1_001_trimmed_tophat_out/Holly_18C.sorted.bam \
#         /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_18C_S13_R1_001_trimmed_tophat_out/accepted_hits.bam
#module load python/2.7.8
#module load htseq/0.6.1p1
#htseq-count -r name -f bam -r name -s no -i gene_id -t transcript \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_18B_redo_S13_R1_001_trimmed_tophat_out/Holly_18B.sorted.bam \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-GC-18B_by.name.txt
#htseq-count -r name -f bam -r name -s no -i gene_id -t transcript \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_18C_S13_R1_001_trimmed_tophat_out/Holly_18C.sorted.bam \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-GC-18C_by.name.txt





##################################################################
#testing cufflinks with stranded and unstranded libraries specified
#18A

module load cufflinks/2.2.1


#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test


#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksHolly-GC-18-redo_S7_R1_001_trimmed_tophat_out

#cufflinks \
#--library-type fr-unstranded -g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksHolly-GC-18-redo_S7_R1_001_trimmed_tophat_out \
#-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-GC-18-redo_S7_R1_001_trimmed_tophat_out/accepted_hits.bam

#18B


mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksHolly_18B_redo_S13_R1_001_trimmed_tophat_out

cufflinks \
--library-type fr-unstranded -g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksHolly_18B_redo_S13_R1_001_trimmed_tophat_out \
-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_18B_redo_S13_R1_001_trimmed_tophat_out/accepted_hits.bam


#18C

mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksHolly_18B_redo_S13_R1_001_trimmed_tophat_out

cufflinks \
--library-type fr-unstranded -g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksHolly_18B_redo_S13_R1_001_trimmed_tophat_out \
-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_18B_redo_S13_R1_001_trimmed_tophat_out/accepted_hits.bam

###################################################
#stranded: first strand
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksFirstStrandHolly-GC-18-redo_S7_R1_001_trimmed_tophat_out


#cufflinks \
#--library-type fr-firststrand -g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksFirstStrandHolly-GC-18-redo_S7_R1_001_trimmed_tophat_out \
#-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-GC-18-redo_S7_R1_001_trimmed_tophat_out/accepted_hits.bam

#18B


mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksFirstStrandHolly_18B_redo_S13_R1_001_trimmed_tophat_out

cufflinks \
--library-type fr-firststrand -g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksFirstStrandHolly_18B_redo_S13_R1_001_trimmed_tophat_out \
-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_18B_redo_S13_R1_001_trimmed_tophat_out/accepted_hits.bam


#18C

mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksFirstStrandHolly_18B_redo_S13_R1_001_trimmed_tophat_out

cufflinks \
--library-type fr-firststrand -g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksFirstStrandHolly_18B_redo_S13_R1_001_trimmed_tophat_out \
-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_18B_redo_S13_R1_001_trimmed_tophat_out/accepted_hits.bam


###################################################
#stranded: second strand
#mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksSecondStrandHolly-GC-18-redo_S7_R1_001_trimmed_tophat_out


#cufflinks \
#--library-type fr-secondstrand -g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksSecondStrandHolly-GC-18-redo_S7_R1_001_trimmed_tophat_out \
#-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-GC-18-redo_S7_R1_001_trimmed_tophat_out/accepted_hits.bam

#18B


mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksSecondStrandHolly_18B_redo_S13_R1_001_trimmed_tophat_out

cufflinks \
--library-type fr-secondstrand -g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksSecondStrandHolly_18B_redo_S13_R1_001_trimmed_tophat_out \
-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_18B_redo_S13_R1_001_trimmed_tophat_out/accepted_hits.bam


#18C

mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksSecondStrandHolly_18B_redo_S13_R1_001_trimmed_tophat_out

cufflinks \
--library-type fr-secondstrand -g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksSecondStrandHolly_18B_redo_S13_R1_001_trimmed_tophat_out \
-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_18B_redo_S13_R1_001_trimmed_tophat_out/accepted_hits.bam


#####################################################################
#cufflinks: no strand specifications

mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksDefHolly-GC-18-redo_S7_R1_001_trimmed_tophat_out


cufflinks \
-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Cufflinks_test/cufflinksDefHolly-GC-18-redo_S7_R1_001_trimmed_tophat_out \
-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly-GC-18-redo_S7_R1_001_trimmed_tophat_out/accepted_hits.bam
