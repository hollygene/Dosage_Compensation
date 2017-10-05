while read SampleName
do
mkdir /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/INDEX${SampleName}

/usr/local/samtools/latest/samtools index \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/${SampleName}/accepted_hits.bam

done < /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/samples.txt
