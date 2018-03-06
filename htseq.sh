#PBS -S /bin/bash
#PBS -q batch
#PBS -N htseq
#PBS -l nodes=1:ppn=4:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=200gb

basedir="/lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/old/Sort"
cd $basedir

module load python/2.7.8
module load htseq/0.6.1p1


mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/old/Sort/HTseq/

while read SampleName
do
  mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/old/Sort/HTseq/${SampleName} \

htseq-count \
  -f bam \
  /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/old/Sort/SORT${SampleName}/${SampleName}.sorted.bam \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/old/Sort/HTseq/${SampleName}.sam \
-i gene_id -t exon \
--additional-attr gene_name \ 
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf

done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/samples.txt
