#PBS -S /bin/bash
#PBS -q batch
#PBS -N samtools_sort
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=200gb


module load samtools/1.3.1
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat

mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/SORTHolly_112B_S28_R1_001_trimmed_tophat_out

time samtools sort \
-m 200G \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/tophat/SORTHolly_112B_S28_R1_001_trimmed_tophat_out/Holly_112B_S28_R1_001_trimmed_tophat_out.sorted.bam \
-O bam \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/MA_new/trimmed/tophat/Holly_112B_S28_R1_001_trimmed_tophat_out/accepted_hits.bam
