#PBS -S /bin/bash
#PBS -q batch
#PBS -N htseq
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l walltime=4:00:00
#PBS -l mem=200gb

cd $PBS_O_WORKDIR

module load python/2.7.8
module load htseq/0.6.1p1

while read SampleName
do
mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/HTseq/

htseq-count \
  -f bam \
  /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/HTseq/${SampleName}.sam \
 /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_test/${BASE}_trimmed_tophat_out/accepted_hits.bam


done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_test/output.txt
