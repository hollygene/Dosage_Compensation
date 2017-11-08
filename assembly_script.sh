#!/bin/bash
#PBS -N assembly
#PBS -q batch
#PBS -l nodes=1:ppn=4:AMD
#PBS -l walltime=96:00:00
#PBS -l pmem=20gb
#PBS -M hmcqueary@uga.edu
#PBS -m ae

basedir="/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly"
cd $basedir
THREADS=4

# put bowtie 2.2.9 executables in $PATH
module load bowtie2/2.2.9
#load in tophat
module load tophat/2.1.1
#load in cufflinks
module load cufflinks/2.2.1
# put the cufflinks 2.2.1 executables in $PATH
module load cufflinks/2.2.1

# build forward and backward indices of reference genome for mapping reads with bowtie
bowtie2-build -f genome.fa genome

#index transcriptome file
tophat -G genes.gtf --transcriptome-index=transcriptome_data/known genome

#run tophat on all samples
for file in ./*.fastq

do

FBASE=$(basename $file .fastq)
BASE=${FBASE%.fastq}

tophat -p $THREADS -o ./tophat/${BASE}_tophat_out -i 10 -I 1000 \
--transcriptome-index=transcriptome_data/known \
genome \
./${BASE}.fastq

done

cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/tophat
printf '%s\n' * > output.txt

#run cufflinks on all samples
while read SampleName
do
mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/Cufflinks/cufflinks${SampleName}
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/Cufflinks
cufflinks \
-p $THREADS
-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/Cufflinks/cufflinks${SampleName} \
-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/tophat/${BASE}_tophat_out

done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/tophat/output.txt

#run cuffquant to estimate expression levels to put into cuffdiff
while read SampleName
do
mkdir /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/quant/cuffquant${SampleName}
cd /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/quant
cuffquant \
-p $THREADS
-g /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf \
-o /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/quant/cuffquant${SampleName} \
-b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genome.fa \
/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/tophat/${BASE}_tophat_out/accepted_hits.bam

done < /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/tophat/output.txt

#run cuffdiff to find differential expression between each line and its ancestor
cuffdiff -p $THREADS --library-type fr-unstranded -o FNR_cuffdiff --labels wild-type,dFNR MG1655.ref.gtf SRR5344681_cuffquant/abundances.cxb,SRR5344682_cuffquant/abundances.cxb SRR5344683_cuffquant/abundances.cxb,SRR5344684_cuffquant/abundances.cxb
