#!/bin/bash
#PBS -N htseq
#PBS -q batch
#PBS -l nodes=1:ppn=4:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=300gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

basedir="/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly"
cd $basedir
THREADS=4


module load HTSeq/0.9.1-foss-2016b-Python-2.7.14
#module load samtools/1.3.1

#cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update

#ls >> samples_update.txt
#while read SampleName
#do
#time samtools sort -m 200G \
#-n \
#-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/${SampleName}/${SampleName}.sorted.bam \
#   /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/${SampleName}/accepted_hits.bam

#done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/samples_update.txt

#GC
#while read SampleName
#do
#htseq-count -r name -f bam -r name -s no -i gene_id -t exon \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/${SampleName}/${SampleName}.sorted.bam \
#/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/HTseq_update/${SampleName}.exon.txt

#done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/samples_update.txt
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat_update
#MA new
while read SampleName
do
htseq-count -r name -f bam -r name -s no -i gene_id -t exon \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat_update/${SampleName}/${SampleName}.sorted.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/HTseq_update/${SampleName}.exon.unstranded.txt

done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat_update/samples_update.txt

#MA old
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed/tophat_update
while read SampleName
do
htseq-count -r name -f bam -r name -s no -i gene_id -t exon \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed/tophat_update/${SampleName}/${SampleName}.sorted.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/HTseq_update/${SampleName}.exon.txt

done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_old/trimmed/tophat_update/samples_update.txt
