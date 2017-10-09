#a	 loop to do all of the files at once (WORKS)
#!/bin/bash

#PBS -q batch
#PBS -l nodes=1:ppn=4:AMD
#PBS -l walltime=96:00:00
#PBS -l pmem=20gb


for file in ./*.fastq

do

FBASE=$(basename $file .fastq)

BASE=${FBASE%.fastq}
THREADS=4
cd $PBS_O_WORKDIR
module load tophat/2.1.1
time  tophat -p $THREADS -o ./${BASE}_tophat_out -i 10 -I 1000 --transcriptome-index=transcriptome_data/known genome /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/${BASE}/${BASE}.fastq

done
