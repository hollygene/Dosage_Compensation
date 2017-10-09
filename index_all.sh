while read SampleName
do
  mkdir /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/INDEX${SampleName}

#!/bin/bash

#PBS -N j_s_samtools
#PBS -q batch
#PBS -l nodes=1:ppn=1:AMD
#PBS -l mem=100gb
#PBS -l walltime=480:00:00

cd $PBS_O_WORKDIR
module load samtools/1.3.1
time  samtools index \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/${SampleName}/accepted_hits.bam

done < /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/samples.txt
 
