#!/bin/bash
#PBS -N htseq
#PBS -q batch
#PBS -l nodes=2:ppn=4:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=200gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

basedir="/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly"
cd $basedir
THREADS=4


module load python/2.7.8
module load htseq/0.6.1p1



#GC
while read SampleName
do
htseq-count \
  -f bam \
  -r pos \
-i gene_id -t exon \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/${SampleName}/${SampleName}.sorted.bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf > /lustre1/hcm14449/SC_RNAseq/RNA_seq/HTseq_update/HTseq_${SampleName}.txt

done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/samples_update.txt
