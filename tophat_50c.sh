#!/bin/bash
#PBS -N mssing_lines
#PBS -q batch
#PBS -l nodes=1:ppn=4:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l pmem=200gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/


module load trimgalore/0.4.4

trim_galore --phred33 -q 20 -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/Holly_MA_Anc_S5_R1_001.fastq

trim_galore --phred33 -q 20 -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/Holly_MA_Anc_B_S1_R1_001.fastq

trim_galore --phred33 -q 20 -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/Holly_50C_S18_R1_001.fastq

trim_galore --phred33 -q 20 -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/SC006_raw_trimmed.fastq

trim_galore --phred33 -q 20 -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/Holly_GC_Anc_S4_R1_001.fastq

module unload trimgalore/0.4.4


module load tophat/2.1.1

tophat -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_MA_Anc_S5_R1_001_trimmed \
-i 10 -I 1000 \
--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/Holly_MA_Anc_S5_R1_001_trimmed.fq


tophat -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_MA_Anc_B_S1_R1_001_trimmed \
-i 10 -I 1000 \
--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/Holly_MA_Anc_B_S1_R1_001_trimmed.fq


tophat -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_50C_S18_R1_001_trimmed \
-i 10 -I 1000 \
--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/Holly_50C_S18_R1_001_trimmed.fq


tophat -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed/tophat/SC006_raw_trimmed \
-i 10 -I 1000 \
--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed/SC006_raw_trimmed.fq


tophat -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat/Holly_GC_Anc_S4_R1_001_trimmed \
-i 10 -I 1000 \
--transcriptome-index=/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/transcriptome_data/known \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/Holly_GC_Anc_S4_R1_001_trimmed.fq

module unload tophat/2.1.1


module load samtools/1.3.1

samtools sort -m 200G \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat/Holly_GC_Anc_S4_R1_001_trimmed.sorted.bam \
   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat/Holly_GC_Anc_S4_R1_001_trimmed/accepted_hits.bam


   samtools sort -m 200G \
   -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_MA_Anc_B_S1_R1_001_trimmed.sorted.bam \
      /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_MA_Anc_B_S1_R1_001_trimmed/accepted_hits.bam

samtools sort -m 200G \
      -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_MA_Anc_S5_R1_001_trimmed.sorted.bam \
         /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_MA_Anc_S5_R1_001_trimmed/accepted_hits.bam

         samtools sort -m 200G \
         -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_50C_S18_R1_001_trimmed.sorted.bam \
            /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_50C_S18_R1_001_trimmed/accepted_hits.bam

            samtools sort -m 200G \
            -o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed/tophat/SC006_raw_trimmed.sorted.bam \
               /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed/tophat/SC006_raw_trimmed/accepted_hits.bam
module unload samtools/1.3.1

   module load python/2.7.8
   module load htseq/0.6.1p1

mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/HTseq2

   htseq-count \
     -f bam \
     -r pos \
   -i gene_id -t exon \
   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat/Holly_GC_Anc_S4_R1_001_trimmed.sorted.bam \
   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/HTseq2/HTseq_Holly_GC_Anc_S4.txt

   htseq-count \
     -f bam \
     -r pos \
   -i gene_id -t exon \
   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed/tophat/SC006_raw_trimmed.sorted.bam \
   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/HTseq2/HTseq_SC006.txt

   htseq-count \
     -f bam \
     -r pos \
   -i gene_id -t exon \
   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_50C_S18_R1_001_trimmed.sorted.bam \
   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/HTseq2/HTseq_Holly_50C.txt

   htseq-count \
     -f bam \
     -r pos \
   -i gene_id -t exon \
   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_MA_Anc_B_S1_R1_001_trimmed.sorted.bam \
   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/HTseq2/HTseq_Holly_MA_Anc_B.txt

   htseq-count \
     -f bam \
     -r pos \
   -i gene_id -t exon \
   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_MA_Anc_S5_R1_001_trimmed.sorted.bam \
   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/HTseq2/HTseq_Holly_MA_Anc.txt

   module unload python/2.7.8
   module unload htseq/0.6.1p1
