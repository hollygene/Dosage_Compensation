#PBS -S /bin/bash
#PBS -q batch
#PBS -N htseq
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=200gb

basedir="/lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/"
cd $basedir

module load python/2.7.8
module load htseq/0.6.1p1

mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/old/Sort/HTseq/Holly_112B_S28_R1_001_trimmed_tophat_out \

htseq-count \
  -f bam \
  /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/SORTHolly_112B_S28_R1_001_trimmed_tophat_out/Holly_112B_S28_R1_001_trimmed_tophat_out.sorted.bam \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/old/Sort/HTseq/Holly_112B_S28_R1_001_trimmed_tophat_out/Holly_112B_S28_R1_001_trimmed_tophat_out.sorted.sam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf
