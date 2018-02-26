#Dosage Compensation Transcriptome Assembly Protocol, April 2017
#to log into the zcluster
#ssh hcm14449@zcluster.rcc.uga.edu
#password
#cd /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/CufflinksOutput

#Need to build an index using Bowtie
#bowtie2-build FILE_A_name.fa FILE_A_name
#FILE_A_name.fa is the reference genome file
#FILE_A_name is the resulting file that has been indexed (this is what you then put into TopHat for mapping)

#reference genome and GTF file obtained from: https://support.illumina.com/sequencing/sequencing_software/igenome.html


#1:
time /usr/local/bowtie2/latest/bin/bowtie2
bowtie2-build -f genome.fa genome
qsub -q bowtie_build.sh

#2 Double check that it worked using bowtie2 inspect
time /usr/local/bowtie2/latest/bin/bowtie2
bowtie2-inspect -s genome
qsub -q rcc-30d bowtie_inspect.sh

#Need to use TopHat to make a transcriptome index from the .gtf file given:
#3: export LD_LIBRARY_PATH=/usr/local/boost/1.54.0/gcc447/lib:/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
time /usr/local/tophat/latest/bin/tophat -G genes.gtf --transcriptome-index=transcriptome_data/known genome

#Then run Tophat with the transcriptome index and parameters needed and reference sequence from bowtie
#4: export LD_LIBRARY_PATH=/usr/local/boost/1.54.0/gcc447/lib:/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
time /usr/local/tophat/latest/bin/tophat -i 10 -I 1000 --transcriptome-index=transcriptome_data/known genome Holly_MA_Anc_B_S1_R1_001.fastq

#a	 loop to do all of the files at once (WORKS)

for file in ./*.fastq

do

FBASE=$(basename $file .fastq)

BASE=${FBASE%.fastq}

tophat -o ./${BASE}_tophat_out -i 10 -I 1000 --transcriptome-index=transcriptome_data/known genome ${BASE}.fastq

done

##A loop to do Cufflinks on all of the tophat outputs at once. Using 5 threads
#!/bin/bash
# Cufflinks script
# reads samples.txt that contains the directory names
# makes a directory for each one (cufflinks${SampleName}
# uses the accepted_hits.bam file for each tophat output as input and runs cufflinks on it
#
#
   while read SampleName
do
	mkdir /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/CufflinksOutput/cufflinks${SampleName}

	/usr/local/cufflinks/latest/bin/cufflinks \
	-g /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/genes.gtf \
	-o /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/CufflinksOutput/cufflinks${SampleName} \
	-b /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/genome.fa \
	-p 5 \
	/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/${SampleName}/accepted_hits.bam

done < /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/samples.txt

#Create a file called assemblies.txt that lists the assembly file for each sample (should look like this:
#cufflinksHolly_112B_S28_R1_001_tophat_out/transcripts.gtf
#cufflinksHolly_112C_S25_R1_001_tophat_out/transcripts.gtf
#cufflinksHolly_115B_S29_R1_001_tophat_out/transcripts.gtf
#cufflinksHolly_115C_S26_R1_001_tophat_out/transcripts.gtf

ls | grep '_tophat_out' > assemblies.txt

#then add /transcripts.gtf to the end of the directory name so that the program knows where to look
awk '{print$"/transcripts.gtf"} > assemblies2.txt

################################################################################################################
################################################################################################################
## script to run cuffmerge

export PATH=/usr/local/cufflinks/latest/bin/:$PATH
time /usr/local/cufflinks/latest/bin/cuffmerge \
-g /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/genes.gtf \
-s /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/genome.fa \
-p 4 \
assemblies2.txt

###on more files
##start with line 1
export PATH=/usr/local/cufflinks/latest/bin/:$PATH
time /usr/local/cufflinks/latest/bin/cuffmerge \
-g /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/genes.gtf \
-s /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/genome.fa \
-p 4 \
X1_assemblies.txt
################################################################################################################
################################################################################################################
#script to run cuffmerge on only one ancestor and one sample file
awk '{print "cufflinksHolly_GC_Anc_S4_R1_001_tophat_out/transcripts.gtf", "cufflinksHolly-GC-8-redo_S5_R1_001_tophat_out/transcripts.gtf"}' > test_assemblies.txt
awk '{print$"/transcripts.gtf"} > test_assemblies2.txt

awk '{print "cufflinksHolly-GC-4-redo_S3_R1_001_tophat_out/transcripts.gtf", "cufflinksHolly-GC-8-redo_S5_R1_001_tophat_out/transcripts.gtf"}' > test_assemblies.txt
awk '{print$"/transcripts.gtf"} > test_assemblies2.txt

################################################################################################################
################################################################################################################

export PATH=/usr/local/cufflinks/latest/bin/:$PATH
time /usr/local/cufflinks/latest/bin/cuffmerge \
-g /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/genes.gtf \
-s /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/genome.fa \
-p 4 \
cufflinksHolly_GC_Anc_S4_R1_001_tophat_out/transcripts.gtf \
cufflinksHolly-GC-8-redo_S5_R1_001_tophat_out/transcripts.gtf

################################################################################################################
################################################################################################################
##script to run cuffdiff on everything

export PATH=/usr/local/cufflinks/latest/bin/:$PATH
time /usr/local/cufflinks/latest/bin/cuffdiff -o diff_out \
-b /panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/genome.fa \
-p 5 \
-L 112,115,117,11,123,141,152,18,1,21,29,2,31,3,49,4,50,59,5,61,66,69,6,76,77,7,8,9,GC_Anc,MA_Anc \
-u merged_asm/merged.gtf \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_MA_112_S8_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_112B_S28_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_112C_S25_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_MA_115_S9_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_115B_S29_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_115C_S26_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_MA_117_S10_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_117B_S30_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_117C_S27_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly-GC-11-redo_S6_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_11B_S15_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_11C_S12_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_MA_123_S11_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_123B_S31_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_123C_S28_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_MA_141_S12_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_141B_S32_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_141C_S29_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_MA_152_S13_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_152B_S33_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_152C_S30_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly-GC-18-redo_S7_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_18B_redo_S13_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_18C_S13_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_GC_1_redo_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_1B_S6_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_1C_S3_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly-GC-21-redo_S8_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_21B_S17_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_21C_S14_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_MA_29_S6_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_29B_S18_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_29C_redo_S17_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_2A_S1_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_2B_S7_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_2C_S4_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly-GC-31-redo_S9_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_31B_S19_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_31C_S16_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_3A_S2_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_3B_redo_S12_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_3C_S5_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly-GC-49-redo_S10_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_49B_S20_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_49C_S17_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly-GC-4-redo_S3_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_4B_S9_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_4C_S6_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_MA_50_S7_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_50B_S21_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_50C_S18_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_GC_59_S10_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_59B_S22_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_59C_S19_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_5A_S3_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_5B_S10_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_5C_redo_S16_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly-GC-61-redo_S11_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_61B_S23_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_61C_S20_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_GC_66_S12_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_66B_S24_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_66C_S21_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_GC_69_S1_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_69B_S25_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_69C_S22_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_6A_S4_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_6B_S11_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_6C_S8_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_GC_76_S2_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_76B_S26_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_76C_redo_S18_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_GC_77_S3_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_77B_S27_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_77C_S24_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly-GC-7-redo_S4_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_7B_S12_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_7C_S9_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly-GC-8-redo_S5_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_8B_S13_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_8C_S10_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_9A_S5_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_9B_S14_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_9C_S11_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_GC_Anc_B_S2_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_GC_Anc_C_S32_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_GC_Anc_S4_R1_001_tophat_out/accepted_hits.bam \
/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_MA_Anc_B_S1_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly-MA-Anc-C-redo_S19_R1_001_tophat_out/accepted_hits.bam,/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/tophat_outputs/Holly_MA_Anc_S5_R1_001_tophat_out/accepted_hits.bam

################################################################################################################
################################################################################################################
#this code replaces the roman numerals in the file with numbers for the chromosome
#need this for downstream analysis in R
awk -F'\t' '{
if ($7~"I" && $7!~"V" && $7!~"X") print $0,"\t","1";
if ($7~"II" && $7!~"III") print $0,"\t","2";
if ($7~"III") print $0,"\t","3";
if ($7~"IV") print $0,"\t","4";
if ($7~"V" && $7!~"I" && $7!~"X") print $0,"\t","5";
if ($7~"VI"&&$7!~"VII"&&$7!~"VIII") print $0,"\t","6";
if ($7~"VII"&&$7!~"VIII") print $0,"\t","7";
if ($7~"VIII") print $0,"\t","8";
if ($7~"IX") print $0,"\t","9";
if ($7~"X" && $7!~"I" && $7!~"V") print $0,"\t","10";
if ($7~"XI"&&$7!~"XII"&&$7!~"XIII") print $0,"\t","11";
if ($7~"XII"&&$7!~"XIII") print $0,"\t","12";
if ($7~"XIII") print $0,"\t","13";
if ($7~"XIV") print $0,"\t","14";
if ($7~"XV"&&$7!~"XVI") print $0,"\t","15";
if ($7~"XVI") print $0,"\t","16";
if ($7~"MT") print $0,"\t","MT";
}' genes.fpkm_tracking.csv > genes.fpkm.csv




##need to join genes.attr_table with genes.fpkm_table
join genes.attr_table.csv genes.fpkm_table.csv > genes.fpkm_table.joined.csv

#doing this again but not .csv format this time 9/25
join genes.attr_table genes.fpkm_table > genes.fpkm_table.joined

#need to go in and replace commas with "/" in original file
tr < genes.fpkm_table.joined.csv "," "/" > genes.fpkm_table.joined.rmcom.csv

#again 9/25
tr < genes.fpkm_table.joined "," "/" > genes.fpkm_table.joined.rmcom

#after files are joined, want to split column that contains the chromosomes
tr < genes.fpkm_table.joined.rmcom.csv ":" "," > genes.fpkm_table.joined.rmcom.split.csv

#9/25
tr < genes.fpkm_table.joined.rmcom ":" "," > genes.fpkm_table.joined.rmcom.split

#add commas where tabs are to make comma separated values file not tab separated values
tr ' ' ',' <genes.fpkm_table.joined.rmcom.split.csv > genes.fpkm_table.joined.rmcom.split.commas.csv

#9/25
tr ' ' ',' < genes.fpkm_table.joined.rmcom.split > genes.fpkm_table.joined.rmcom.split.commas.csv

##remove the mitochondrial sequences
awk -F',' '{if ($7!~"MT") print $0}' genes.fpkm_table.joined.rmcom.split.commas.csv > genes.fpkm_table.joined.rmcom.split.commas.rmvmito.csv

#9/25
awk -F',' '{if ($7!~"MT") print $0}' genes.fpkm_table.joined.rmcom.split.commas.csv > genes.fpkm_table.joined.rmcom.split.commas.nomito.csv


#need this for downstream analysis in R
##awk -F',' '{
#if (NR>1) {
#if ($7~"I" && $7!~"V" && $7!~"X" && $7!~"II" && $7!~"III") print "1", ",", $0;
#if ($7~"II" && $7!~"III" && $7!~"X" && $7!~"V") print "2", ",", $0;
#if ($7~"III" && $7!~"V" && $7!~"X") print "3", ",", $0;
#if ($7~"IV" && $7!~"X") print "4", ",", $0;
#if ($7~"V" && $7!~"I" && $7!~"X") print "5", ",", $0;
#if ($7~"VI" && $ 7!~"VII" && $7!~"VIII" && $7!~"XVI") print "6", ",", $0;
#if ($7~"VII" && $7!~"VIII") print "7", ",", $0;
#if ($7~"VIII") print "8", ",", $0;
#if ($7~"IX") print "9", ",", $0;
#if ($7~"X" && $7!~"I" && $7!~"V") print "10", ",", $0;
#if ($7~"XI" && $7!~"XII" && $7!~"XIII" && $7!~"V") print "11", ",", $0;
#if ($7~"XII" && $7!~"XIII") print "12", ",", $0;
#if ($7~"XIII") print "13", ",", $0;
#if ($7~"XIV") print "14", ",", $0;
#if ($7~"XV" && $7!~"XVI") print "15", ",", $0;
#if ($7~"XVI") print "16", ",", $0  } else print $0
#}' genes.fpkm_table.joined.split.commas.rmvmito.csv > genes.fpkm_table.joined.split.commas.rmvmito.chgchr.csv

##this one works just fine
awk -F',' '{
if (NR>1) {
if ($7=="I") print "1",",",$0;
if ($7=="II") print "2", ",", $0;
if ($7=="III") print "3", ",", $0;
if ($7=="IV") print "4", ",", $0;
if ($7=="V") print "5", ",", $0;
if ($7=="VI") print "6", ",", $0;
if ($7=="VII") print "7", ",", $0;
if ($7=="VIII") print "8", ",", $0;
if ($7=="IX") print "9", ",", $0;
if ($7=="X") print "10", ",", $0;
if ($7=="XI") print "11", ",", $0;
if ($7=="XII") print "12", ",", $0;
if ($7=="XIII") print "13", ",", $0;
if ($7=="XIV") print "14", ",", $0;
if ($7=="XV") print "15", ",", $0;
if ($7=="XVI") print "16", ",", $0  } else print $0
}' genes.fpkm_table.joined.rmcom.split.commas.nomito.csv> genes.fpkm_table.joined.rmcom.split.commas.nomito.chgchr.csv

awk -F',' '{
if (NR>1) {
if ($7=="I") print "1", ",", $0;
if ($7=="II") print "2", ",", $0;
if ($7=="III") print "3", ",", $0;
if ($7=="IV") print "4", ",", $0;
if ($7=="V") print "5", ",", $0;
if ($7=="VI") print "6", ",", $0;
if ($7=="VII") print "7", ",", $0;
if ($7=="VIII") print "8", ",", $0;
if ($7=="IX") print "9", ",", $0;
if ($7=="X") print "10", ",", $0;
if ($7=="XI") print "11", ",", $0;
if ($7=="XII") print "12", ",", $0;
if ($7=="XIII") print "13", ",", $0;
if ($7=="XIV") print "14", ",", $0;
if ($7=="XV") print "15", ",", $0;
if ($7=="XVI") print "16", ",", $0  } else print $0
}' genes.fpkm_table.joined.rmcom.split.commas.nomito.csv > genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.csv



##NOW I HAVE THE RIGHT NUMBER OF LINES BETWEEN MY TWO FILES
#Now I need to figure out how to make my column names match
#How to print one string in one column in one line?
#Or just make a new header?
#this code works but I have to put it into textwrangler to get it all on one line so that command line likes it
#awk 'BEGIN{print #"chr",",","tracking_id",",","class_code",",","nearest_ref_id",",","gene_id",",","gene_short_name",",","tss_id",",","chr_roman",",","length",",","null",",","X112_2",",","X112_0",",","X112_1",",","X115_0",",","X115_1",",","X117_0",",","X117_1",",","X117_2",",","X11_1",",","X11_0",",","X11_2",",","X123_1",",","X123_0,",","X123_2",",","X141_1",",","X141_0",",","X141_2",",","X152_1",",","X152_0",",","X152_2",",","X18_1",",","X18_0",",","X18_2",",","X1_0",",","X1_1",",","X1_2",",","X21_1",",","X21_0"#,",","X21_2",",","X29_1",",","X29_0",",","X29_2",",","X2_0",",","X2_1",",","X2_2",",","X31_1",",","X31_0",",","X31_2",",","X3_0",",","X3_1",",","X3_2",",","X49_0",",","X49_1",",","X49_2",",","X4_1",",","X4_0",",","X4_2",",","X50_1",",","X50_0",",","X50_2",",","X59_0",",","X59_1",",","X59_2",",","X5_1",",","X5_0",",","X5_2",",","X61_1",",","X61_0",",","X61_2",",","X66_0",",","X66_1",",","X66_2",",","X69_1",",","X69_0",",","X69_2",",","X6_0",",","X6_1",",","X6_2",",","X76_0",","
#,"X76_1",",","X76_#2",",","X77_1",",","X77_0",",","X77_2",",","X7_1",",","X7_0",",","X7_2",",","X8_0",",","X_1",",","X8_2",",","X9_0",",","X9_1",",","X9_2",",","GC_Anc_1",",",GC_Anc_0",",","GC_Anc_2",",","MA_Anc_1",",","MA_Anc_0",",","MA_Anc_2"}; {print}' genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.csv > genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.nwhead.csv


#get rid of old header
grep -v "locus"  genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.csv > temp && mv temp genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.csv

#9/25
grep -v "locus"  genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.csv > temp && mv temp genes.fpkm_table.joined.rmcom.split.commas.nomito.chgchr.rmhead.csv

#add new header
echo -e "chr,tracking_id,class_code,nearest_ref_id,gene_id,gene_short_name,tss_id,chr_roman,length,null,X112_2,X112_0,X112_1,X115_0,X115_1,X115_2,X117_0,X117_1,X117_2,X11_1,X11_0,X11_2,X123_1,X123_0,X123_2,X141_1,X141_0,X141_2,X152_1,X152_0,X152_2,X18_1,X18_0,X18_2,X1_0,X1_1,X1_2,X21_1,X21_0,X21_2,X29_1,X29_0,X29_2,X2_0,X2_1,X2_2,X31_1,X31_0,X31_2,X3_0,X3_1,X3_2,X49_0,X49_1,X49_2,X4_1,X4_0,X4_2,X50_1,X50_0,X50_2,X59_0,X59_1,X59_2,X5_1,X5_0,X5_2,X61_1,X61_0,X61_2,X66_0,X66_1,X66_2,X69_1,X69_0,X69_2,X6_0,X6_1,X6_2,X76_0,X76_1,X76_2,X77_1,X77_0,X77_2,X7_1,X7_0,X7_2,X8_0,X8_1,X8_2,X9_0,X9_1,X9_2,GC_Anc_1,GC_Anc_0,GC_Anc_2,MA_Anc_1,MA_Anc_0,MA_Anc_2" | cat - genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.csv > genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.csv

#9/25
echo -e "chr,tracking_id,class_code,nearest_ref_id,gene_id,gene_short_name,tss_id,chr_roman,length,null,SCA1_0,SCA1_2,SCA1_1,SC001_0,SC001_1,SC001_2,SC002_0,SC002_2,SC002_1,SC003_2,SC003_0,SC003_1,SC004_0,SC004_1,SC004_2,SC005_2,SC005_0,SC005_1,SC006_0,SC006_1,SC006_2,SC007_0,SC007_1,SC007_2,SC008_0,SC008_1,SC008_2,SC009_2,SC009_0,SC009_1,SC011_0,SC011_1,SC011_2,SC015_0,SC015_1,SC015_2,SC028_0,SC028_2,SC028_1,SC088_0,SC088_1,SC088_2,SC108_0,SC108_1,SC108_2,SC119_0,SC119_2,SC119_1" | cat - genes.fpkm_table.joined.rmcom.split.commas.nomito.chgchr.rmhead.csv> genes.fpkm_table.joined.rmcom.split.commas.nomito.chgchr.rmhead.nwhd.csv

#need to get rid of columns that I don't actually need
#don't need: 3-10
cut -f1,2,11- -d',' genes.fpkm_table.joined.rmcom.split.commas.nomito.chgchr.rmhead.nwhd.csv > genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.csv

#9/25
#getting rid of anything not MA also
cut -f1,2,5,7,9,11-19,23-31,41-43,59-61,98-100 -d',' genes.fpkm_table.joined.rmcom.split.commas.nomito.chgchr.rmhead.nwhd.csv > genes.fpkm_table.joined.rmcom.split.commas.nomito.chgchr.rmhead.nwhd.MA.csv


#in old data, need to get rid of mitochondrial data and also convert the file to .csv not .tsv
awk -F',' '{if ($1!~"chrmt") print $0}' gene_fpkm_chrm.tsv > genes_fpkm_chrm.nomito.tsv
tr < genes_fpkm_chrm.nomito.tsv "," "/" > genes_fpkm_chrm.nomito.rmcommas.tsv
#need to put back commas into old data to make it actually comma separated values
##this one works!!!
for file in genes_fpkm_chrm.nomito.tsv; do cat $file| tr '[\t]' '[,]' > genes_fpkm_chrm.nomito.commas.tsv; done
##change file to csv
cp genes_fpkm_chrm.nomito.commas.tsv genes_fpkm_chrm.nomito.commas.csv
##
awk 'NR==FNR {a[$2]; next} $2 in a {print}'  genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.csv genes_fpkm_chrm.nomito.commas.csv > genes_fpkm_match_oldtonew.csv
##concatenate them together somehow
#this ALMOST works, but messes up some lines... not sure why yet.will continue tomorrow.
#pr -mts, genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.csv genes_fpkm_match_oldtonew.csv > genes_fpkm_merge.csv
#trying join
cp genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.csv genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.csv_bak
cp genes_fpkm_match_oldtonew.csv genes_fpkm_match_oldtonew.csv_bak
#join -1 2 -2 2 -o 1.1, 1.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 2.10, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.20, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.30, 2.31, 2.32, 2.33, 2.34, 2.35, 2.36, 2.37, 2.38, 2.39, 2.40, 2.41, 2.42, 2.43, 2.44, 2.45, 2.46, 2.47, 2.48, 2.49, 2.50, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 1.10, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.20, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.30, 1.31, 1.32, 1.33, 1.34, #1.35, 1.36, 1.37, 1.38, 1.39, 1.40, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.50, 1.51, 1.52, 1.53, 1.54, 1.55, 1.56, 1.57, 1.58, 1.59, 1.60, 1.61, 1.62, 1.63, 1.64, 1.65, 1.66, 1.67, 1.68, 1.69, 1.70, 1.71, 1.72, 1.73, 1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 1.80, 1.81, 1.82, 1.83, 1.84, 1.85, 1.86, 1.87, 1.88, 1.89, 1.90, 1.91, 1.92 -t , genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.csv genes_fpkm_match_oldtonew.csv > genes_fpkm_merge.csv
#join -1 2 -2 2 -o 1.1, 1.2, 2.3 -t , <genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.csv genes_fpkm_match_oldtonew.csv > genes_fpkm_merge.csv
#maybe my files need to be sorted
#sort -k 2 -t , genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.csv > genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.sorted.csv
#doesnt work
#(head -n 1 genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.csv | sort -k 2 -t ,) > genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.sorted.csv
#trying awk -WORKS
awk 'NR<2{print $0;next}{print$0 | "sort -k 1 -t , -g "}' genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.csv  > genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.sorted.csv
awk 'NR<2{print $0;next}{print$0 | "sort -k 1 -t , -g "}' genes_fpkm_match_oldtonew.csv > genes_fpkm_match_oldtonew_sorted.csv
#sorted on column 2 since I want to join them on column 2
awk 'NR<2{print $0;next}{print$0 | "sort -k 2 -t , -f "}' genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.csv  > genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.sorted2.csv
awk 'NR<2{print $0;next}{print$0 | "sort -k 2 -t , -f "}' genes_fpkm_match_oldtonew.csv > genes_fpkm_match_oldtonew_sorted2.csv

#trying to join them again
#join -1 2 -2 2 -o #2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18,2.19,2.20,2.21,2.22,2.23,2.24,2.25,2.26,2.27,2.28,2.29,2.30,2.31,2.32,2.33,2.34,2.35,2.36,2.37,2.38,2.39,2.40,2.41,2.42,2.43,2.44,2.45,2.46,2.47,2.48,2.49,2.50,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,#1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85,1.86,1.87,1.88,1.89,1.90,1.91,1.92 -t , genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.sorted.csv genes_fpkm_match_oldtonew_sorted.csv > genes_fpkm_sorted_merge.csv
#join -1 2 -2 2 -t , genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.sorted.csv genes_fpkm_match_oldtonew_sorted.csv > genes_fpkm_sorted_merge.csv
#trying awk instead
#awk 'FNR==NR{a[$2]=$2$3;next} $2 in a {print $0, a[$2]}' file2 file1 > output
#using paste-WORKS
paste -d"," genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.sorted.csv genes_fpkm_match_oldtonew_sorted.csv > genes_fpkm_sorted_paste.csv
#making sure there are the same number of lines when I do this with the other sort method (sort by column 2)
paste -d"," genes.fpkm_table.joined.rmcom.split.commas.rmvmito.chgchr.rmhead.nwhd.shrt.sorted2.csv genes_fpkm_match_oldtonew_sorted2.csv > genes_fpkm_sorted_paste2.csv
#yes, so probably best to sort by chromosome
##now I need to get rid of the two columns that I don't need from the old data
cut -f1-92,95- - -d',' genes_fpkm_sorted_paste.csv > genes_fpkm_sorted_paste_colrm.csv
#these don't work
#tr '	' ',' < genes_fpkm_chrm.nomito.csv > genes_fpkm_chrm.nomito.commas.csv
#tr ' ' '\t' < genes_fpkm_chrm.nomito.csv > genes_fpkm_chrm.nomito.commas.csv
##doesn't work
#sed -i -e 's/\t/,/g'  genes_fpkm_chrm.nomito.csv > genes_fpkm_chrm.nomito.commas.csv


##############################################################################################################
#these codes are unnecessary now, were to compare what was in my two files
#print column 2 for comparison
#awk -F"," '{print $7,$8}' genes.fpkm_table.joined.split.commas.rmvmito.csv > comparison.rmvmito.csv
#awk -F"," '{print $8,$9}' genes.fpkm_table.joined.split.commas.rmvmito.chgchr.csv > comparison.rmvmito.chgchr.csv
##find difference between two files
#awk -F"," 'NR==FNR{a[$1]=$2;next}{if (a[$1])print a[$1],$0;else print "Not Found", $0;}' OFS=";" comparison.rmvmito.csv comparison.rmvmito.chgchr.csv > difference.csv
#grep -vFxf comparison.rmvmito.csv comparison.rmvmito.chgchr.csv > difference.csv
#grep -vFxf comparison.rmvmito.chgchr.csv comparison.rmvmito.csv  > difference.csv
##reason for my files being different lengths was that some lines had a gene id in which
##there was a comma in it so it made a new column
#need to go in and replace commas with "/" in original file
#tr < genes.fpkm_table.joined.csv "," "/" > genes.fpkm_table.joined.rmcom.csv
#awk -F'\t' '{
#if ($1=="I" && $1!~"V" && $1!~"X" && $1!~"II" && $1!~"III") print "1", ",", $0;
#}' genes.fpkm_table.csv > test.csv
#gsub(pattern, replacement, x,
#     ignore.case = FALSE, extended = TRUE, perl = FALSE,
#     fixed = FALSE, useBytes = FALSE)
#     awk '{ gsub(/I/&&!/II/, "1"); print }' genes.fpkm_table.csv > genes.fpkm_table_2.csv
#     \t
#     sed -i '' 's/I'\t'/1'\t'/g'  genes.fpkm_table.csv > genes.fpkm_table_2.csv
###excel formula to change roman numerals into numbers
#=MATCH(B4,INDEX(ROMAN(ROW(INDIRECT("1:4000"))),0),0)
#awk -F'\t' '{ if $1!~"#N/A" print $0 }' genes.fpkm_table.txt > genes.fpkm_table_nomito.txt

################################################################################################################
################################################################################################################
###need to make a code to compare files from Cuffdiff to files from edgeR of differentially expressed genes
###want to find those genes that are common and different between the two
##awk -F'\t' 'NR==FNR{c[$3]++;next};c[$2] > 0' gene_exp_112.diff X112.vs.anc.csv > 112_diff_genes.txt
###awk 'FNR==NR {arr[$3];next} $2 in arr' gene_exp_112.diff X112.vs.anc.csv > 112_diff_genes.txt
##this code finds the lines where the p-value is less than or equal to 0.05 and also includes the header in the file
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' gene_exp_112.diff > sig.genes.112.cuffdiff.txt

##do the same for the edgeR file
##have to get the .csv file into a .txt file
##some of the commented out scripts below don't work for that
###awk '{gsub(/\,/,'\t');print;}' X112.vs.anc.csv >X112.vs.anc.txt
###sed 's|,|\t|g' X112.vs.anc.csv >X112.vs.anc.txt
####sed -e 's/,/\t/g' X112.vs.anc.csv >X112.vs.anc.txt
###this one works

cat X112.vs.anc.csv | tr "," "\\t" > X112.vs.anc.txt
awk -F'\t' '{
if (NR>1) {
if ($8<="0.05") print $0 } else print $0}' 112vsanc_noquotes.txt > sig.genes.112.edgeR.txt

###awk 'FNR==NR {arr[$2];next} $3 in arr' sig.genes.112.cuffdiff.txt sig.genes.112.edgeR.txt > 112_diff_genes.txt
###awk -F'\t' 'NR==FNR{c[$3]++;next};c[$2] > 0' sig.genes.112.cuffdiff.txt sig.genes.112.edgeR.txt > 112_diff_genes.txt
###grep -f sig.genes.112.cuffdiff.txt sig.genes.112.edgeR.txt > 112_diff_genes.txt
##awk -F, 'FNR==NR {a[$3]; next}; $2 in a' sig.genes.112.edgeR.txt sig.genes.112.cuffdiff.txt  > 112_diff_genes.txt
###awk -F, 'FNR==NR {a[$2]; next}; $3 in a'  112_diff_genes.txt sig.genes.112.edgeR.txt  > 112_diff_genes_rev.txt
###comm -3 <(sort 112_diff_genes.txt) <(sort sig.genes.112.cuffdiff.txt) > differences.txt
##comm -3 <(sort 112_diff_genes.txt) <(sort sig.genes.112.edgeR.txt) > differences.txt
###awk 'NR==FNR {a[$3]=$2; next} $2 in a {print $0, a[$4,$5,$6,$7,$8,$9,$10]}' OFS='\t' sig.genes.112.edgeR.txt sig.genes.112.cuffdiff.txt  > 112_diff_genes.txt
##awk 'NR==FNR {a[$3]=$2; next} $2 in a {print $0,a[$3],a[$8],a[$9],a[$10],a[$11],a[$12]}' sig.genes.112.edgeR.txt sig.genes.112.cuffdiff.txt  > 112_diff_genes_combo.txt
##awk 'NR==FNR {a[$3]=$2; next} $2 in a {print $0}' sig.genes.112.edgeR.txt sig.genes.112.cuffdiff.txt  > 112_diff_genes_combo.txt

awk '{print $3}' sig.genes.112.cuffdiff.txt > genes.cuffdiff.txt

awk '{print $2}' sig.genes.112.edgeR.txt > genes.edgeR.txt

awk 'NR==FNR {a[$2]; next} $1 in a {print}' sig.genes.112.edgeR.txt genes.cuffdiff.txt > 112_diff_genestest.txt

###THIS ONE WORKS
##NOT SURE WHY YET BUT IT DOES

awk 'NR==FNR {a[$2]; next} $3 in a {print}' sig.genes.112.edgeR.txt sig.genes.112.cuffdiff.txt > 112_diff_genes.txt

###awk 'NR==FNR {a[$3]=$2;next} $2 in a {print $0, a[$3],a[$8],a[$9],a[$10],a[$11],a[$12]}' OFS='\t' sig.genes.112.cuffdiff.txt sig.genes.112.edgeR.txt > 112_diff_genestest4.txt
####{print $0,a[$3],a[$8],a[$9],a[$10],a[$11],a[$12]
###gawk 'NR==FNR {a[$2][$2]++; next} $3 in a {for (x in a[$1]) print $0, x}' OFS="\t" file_2.txt file_1.txt
##trying to remove the quotes
##works

sed 's/\"//g' X112.vs.anc.txt > test_noquotes.txt

################################################################################################################
################################################################################################################
###pipeline, change for each file
##prints only those lines whose p-value is less than 0.05
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' gene_exp_112.diff > sig.genes.112.cuffdiff.txt
##removes commas
cat X112.vs.anc.csv | tr "," "\\t" > X112.vs.anc.txt
###removes quotes
sed 's/\"//g' X112.vs.anc.txt > test_noquotes.txt
##prints only those lines whose p-value is less than 0.05 (edgeR file)
awk -F'\t' '{
if (NR>1) {
if ($8<="0.05") print $0 } else print $0}' 112vsanc_noquotes.txt > sig.genes.112.edgeR.txt
##compares the cuffdiff and edgeR files and prints the matching lines from the second file (cuffdiff file)
awk 'NR==FNR {a[$2]; next} $3 in a {print}' sig.genes.112.edgeR.txt sig.genes.112.cuffdiff.txt > 112_diff_genes.txt
##confirms that the same number of genes are found when reversing the order of the files
##works, same number of genes identified.
awk 'NR==FNR {a[$3]; next} $2 in a {print}'  sig.genes.112.cuffdiff.txt sig.genes.112.edgeR.txt > 112_diff_genes_confirm.txt

##for line 115
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' gene_exp_115.diff > sig.genes.115.cuffdiff.txt
##removes commas
cat X115.vs.anc.csv | tr "," "\\t" > X115.vs.anc.txt
###removes quotes
sed 's/\"//g' X115.vs.anc.txt > 115vsanc_noquotes.txt
##prints only those lines whose p-value is less than 0.05 (edgeR file)
awk -F'\t' '{
if (NR>1) {
if ($8<="0.05") print $0 } else print $0}' 115vsanc_noquotes.txt > sig.genes.115.edgeR.txt
##compares the cuffdiff and edgeR files and prints the matching lines from the second file (cuffdiff file)
awk 'NR==FNR {a[$2]; next} $3 in a {print}' sig.genes.115.edgeR.txt sig.genes.115.cuffdiff.txt > 115_diff_genes.txt
##confirms that the same number of genes are found when reversing the order of the files
##works, same number of genes identified.
awk 'NR==FNR {a[$3]; next} $2 in a {print}'  sig.genes.115.cuffdiff.txt sig.genes.115.edgeR.txt > 115_diff_genes_confirm.txt



##for line 117
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' gene_exp_117.diff > sig.genes.117.cuffdiff.txt

##removes commas
cat X117.vs.anc.csv | tr "," "\\t" > X117.vs.anc.txt

###removes quotes
sed 's/\"//g' X117.vs.anc.txt > 117vsanc_noquotes.txt

##prints only those lines whose p-value is less than 0.05 (edgeR file)
awk -F'\t' '{
if (NR>1) {
if ($8<="0.05") print $0 } else print $0}' 117vsanc_noquotes.txt > sig.genes.117.edgeR.txt

##compares the cuffdiff and edgeR files and prints the matching lines from the second file (cuffdiff file)
awk 'NR==FNR {a[$2]; next} $3 in a {print}' sig.genes.117.edgeR.txt sig.genes.117.cuffdiff.txt > 117_diff_genes.txt

##confirms that the same number of genes are found when reversing the order of the files
##works, same number of genes identified.
awk 'NR==FNR {a[$3]; next} $2 in a {print}'  sig.genes.117.cuffdiff.txt sig.genes.117.edgeR.txt > 117_diff_genes_confirm.txt


##for line 117
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' gene_exp_X21.diff > sig.genes.21.cuffdiff.txt

##removes commas
cat X21.vs.anc.csv | tr "," "\\t" > X21.vs.anc.txt

###removes quotes
sed 's/\"//g' X21.vs.anc.txt > 21vsanc_noquotes.txt

##prints only those lines whose p-value is less than 0.05 (edgeR file)
awk -F'\t' '{
if (NR>1) {
if ($8<="0.05") print $0 } else print $0}' 21vsanc_noquotes.txt > sig.genes.21.edgeR.txt

##compares the cuffdiff and edgeR files and prints the matching lines from the second file (cuffdiff file)
awk 'NR==FNR {a[$2]; next} $3 in a {print}' sig.genes.21.edgeR.txt sig.genes.21.cuffdiff.txt > 21_diff_genes.txt

##confirms that the same number of genes are found when reversing the order of the files
##works, same number of genes identified.
awk 'NR==FNR {a[$3]; next} $2 in a {print}'  sig.genes.21.cuffdiff.txt sig.genes.21.edgeR.txt > 21_diff_genes_confirm.txt

##compare the .diff file to the fpkm file to find genes that are not differentially expressed between the sample and ancestor
awk 'NR==FNR {a[$2]; next} !$4 in a {print}' X4.gene_exp.diff.txt X4_genes.fpkm_tracking.txt  > X4_not_diff_exp.txt

#grep -v -f X4.gene_exp.diff.txt X4_genes.fpkm_tracking.txt > X4_not_diff_exp.txt

#comm -23 X4.gene_exp.diff.txt X4_genes.fpkm_tracking.txt > X4_not_diff_exp.txt
##comm outputs 3 columns: left-only, right-only, both. the -1, -2, -3 switches suppress these columns
##-23 hides the right-only and both columns, showing the lines that appear only in the first (left) file

#diff  X4_genes.fpkm_tracking.txt X4.gene_exp.diff.txt > X4_not_diff_exp.txt | grep \^\<
##compare the two files and output only lines that are in the left one but not the right one
##such lines are flagged by diff with < so it suffices to grep that symbol at the beginning of the line

######################################################################################################################################

##find genes that are not differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12>="0.05") print $0 } else print $0}' X4.gene_exp.diff.txt > X4_not_diff_exp.txt

#find genes that are significantly differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' X4.gene_exp.diff.txt > X4_diff_exp.txt

#find genes that are not significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if ($4~"V" && $4!~"I" && $4!~"X" ) print $0 }' X4_not_diff_exp.txt > X4_diff_exp.chr5.txt

#find genes that are significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if ($4~"V" && $4!~"I" && $4!~"X" ) print $0 }' X4_diff_exp.txt > X4_diff_exp.chr5.txt

######################################################

##find genes that are not differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12>="0.05") print $0 } else print $0}' 1_gene_exp.diff.txt > 1_not_diff_exp.txt

#find genes that are significantly differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' 1_gene_exp.diff.txt > 1_diff_exp.txt

#find genes that are not significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) {  if ($4~"XII" && $4!~"XIII" ) print $0 } else print $0}' 1_not_diff_exp.txt > 1_not_diff_exp.chr12.txt

#find genes that are significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) { if ($4~"XII" && $4!~"XIII" ) print $0 }else print $0}' 1_diff_exp.txt > 1_diff_exp.chr12.txt


######################################################


##add .txt to end of filename and also include line number
cp gene_exp.diff  X2_gene_exp.diff.txt

##find genes that are not differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12>="0.05") print $0 } else print $0}' X2_gene_exp.diff.txt > X2_not_diff_exp.txt

#find genes that are significantly differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' X2_gene_exp.diff.txt > X2_diff_exp.txt

######################################################


##add .txt to end of filename and also include line number
cp gene_exp.diff  X3_gene_exp.diff.txt

##find genes that are not differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12>="0.05") print $0 } else print $0}' X3_gene_exp.diff.txt > X3_not_diff_exp.txt

#find genes that are significantly differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' X3_gene_exp.diff.txt > X3_diff_exp.txt

######################################################


##add .txt to end of filename and also include line number
cp gene_exp.diff  X5_gene_exp.diff.txt

##find genes that are not differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12>="0.05") print $0 } else print $0}' X5_gene_exp.diff.txt > X5_not_diff_exp.txt

#find genes that are significantly differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' X5_gene_exp.diff.txt > X5_diff_exp.txt


######################################################


##add .txt to end of filename and also include line number
cp gene_exp.diff  X6_gene_exp.diff.txt

##find genes that are not differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12>="0.05") print $0 } else print $0}' X6_gene_exp.diff.txt > X6_not_diff_exp.txt

#find genes that are significantly differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' X6_gene_exp.diff.txt > X6_diff_exp.txt

######################################################
##add .txt to end of filename and also include line number
cp gene_exp.diff  X7_gene_exp.diff.txt

##find genes that are not differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12>="0.05") print $0 } else print $0}' X7_gene_exp.diff.txt > X7_not_diff_exp.txt

#find genes that are significantly differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' X7_gene_exp.diff.txt > X7_diff_exp.txt

#find genes that are not significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) {  if ($4~"VII" && $4!~"VIII" ) print $0 } else print $0}' X7_not_diff_exp.txt > X7_not_diff_exp.COI.txt

#find genes that are significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) { if ($4~"VII" && $4!~"VIII" ) print $0 }else print $0}' X7_diff_exp.txt > X7_diff_exp.chrCOI.txt

######################################################
##add .txt to end of filename and also include line number
cp gene_exp.diff  X141_gene_exp.diff.txt

##find genes that are not differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12>="0.05") print $0 } else print $0}' X141_gene_exp.diff.txt > X141_not_diff_exp.txt

#find genes that are significantly differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' X141_gene_exp.diff.txt > X141_diff_exp.txt

#find genes that are not significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) {  if ($4~"XVI" ) print $0 } else print $0}' X141_not_diff_exp.txt > X141_not_diff_exp.COI.txt

#find genes that are significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) { if ($4~"XVI") print $0 }else print $0}' X141_diff_exp.txt > X141_diff_exp.chrCOI.txt


######################################################
##add .txt to end of filename and also include line number
cp gene_exp.diff  X9_gene_exp.diff.txt

##find genes that are not differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12>="0.05") print $0 } else print $0}' X9_gene_exp.diff.txt > X9_not_diff_exp.txt

#find genes that are significantly differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' X9_gene_exp.diff.txt > X9_diff_exp.txt

######################################################
##add .txt to end of filename and also include line number
cp gene_exp.diff  X11_gene_exp.diff.txt

##find genes that are not differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12>="0.05") print $0 } else print $0}' X11_gene_exp.diff.txt > X11_not_diff_exp.txt

#find genes that are significantly differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' X11_gene_exp.diff.txt > X11_diff_exp.txt

#find genes that are not significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) {  if ($4~"I" && $4!~"V" && $4!~"X" && $4!~"II" && $4!~"III") print $0 } else print $0}' X11_not_diff_exp.txt > X11_not_diff_exp.COI1.txt

#find genes that are significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) { if ($4~"I" && $4!~"V" && $4!~"X" && $4!~"II" && $4!~"III") print $0 }else print $0}' X11_diff_exp.txt > X11_diff_exp.chrCOI1.txt

#find genes that are not significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) {  if ($4~"XV" && $4!~"I") print $0 } else print $0}' X11_not_diff_exp.txt > X11_not_diff_exp.COI2.txt

#find genes that are significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) { if ($4~"XV" && $4!~"I") print $0 }else print $0}' X11_diff_exp.txt > X11_diff_exp.chrCOI2.txt

######################################################
##add .txt to end of filename and also include line number
cp gene_exp.diff  X123_gene_exp.diff.txt

##find genes that are not differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12>="0.05") print $0 } else print $0}' X123_gene_exp.diff.txt > X123_not_diff_exp.txt

#find genes that are significantly differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' X123_gene_exp.diff.txt > X123_diff_exp.txt

#find genes that are not significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) {  if ($4~"XII" && $4!~"V" && $4!~"XIII") print $0 } else print $0}' X123_not_diff_exp.txt > X123_not_diff_exp.COI1.txt

#find genes that are significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) { if ($4~"XII" && $4!~"V" && $4!~"XIII") print $0 }else print $0}' X123_diff_exp.txt > X123_diff_exp.chrCOI1.txt


######################################################
##add .txt to end of filename and also include line number
cp gene_exp.diff  X21_gene_exp.diff.txt

##find genes that are not differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12>="0.05") print $0 } else print $0}' X21_gene_exp.diff.txt > X21_not_diff_exp.txt

#find genes that are significantly differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' X21_gene_exp.diff.txt > X21_diff_exp.txt

#find genes that are not significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) {  if ($4~"I" && $4!~"V" && $4!~"X" && $4!~"II" && $4!~"III") print $0 } else print $0}' X21_not_diff_exp.txt > X21_not_diff_exp.COI1.txt

#find genes that are significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) { if ($4~"I" && $4!~"V" && $4!~"X" && $4!~"II" && $4!~"III") print $0 }else print $0}' X21_diff_exp.txt > X21_diff_exp.chrCOI1.txt


######################################################
##add .txt to end of filename and also include line number
cp gene_exp.diff  X76_gene_exp.diff.txt

##find genes that are not differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12>="0.05") print $0 } else print $0}' X76_gene_exp.diff.txt > X76_not_diff_exp.txt

#find genes that are significantly differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' X76_gene_exp.diff.txt > X76_diff_exp.txt

#find genes that are not significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) {  if ($4~"VII" && $4!~"VIII" )print $0 } else print $0}' X115_not_diff_exp.txt > X115_not_diff_exp.COI1.txt

#find genes that are significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) { if ($4~"VII" && $4!~"VIII" )print $0 }else print $0}' X115_diff_exp.txt > X115_diff_exp.chrCOI1.txt


#find genes that are not significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) {  if ($4~"X" && $4!~"V" && $4!~"I")print $0 } else print $0}' X76_not_diff_exp.txt > X76_not_diff_exp.COI2.txt

#find genes that are significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) { if ($4~"X" && $4!~"V" && $4!~"I")print $0 }else print $0}' X76_diff_exp.txt > X76_diff_exp.chrCOI2.txt

#find genes that are not significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) {  if ($4~"XIV")print $0 } else print $0}' X76_not_diff_exp.txt > X76_not_diff_exp.COI3.txt

#find genes that are significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) { if ($4~"XIV")print $0 }else print $0}' X76_diff_exp.txt > X76_diff_exp.chrCOI3.txt


######################################################
##add .txt to end of filename and also include line number
cp gene_exp.diff  X152_gene_exp.diff.txt

##find genes that are not differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12>="0.05") print $0 } else print $0}' X152_gene_exp.diff.txt > X152_not_diff_exp.txt

#find genes that are significantly differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' X152_gene_exp.diff.txt > X152_diff_exp.txt

#find genes that are not significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) {  if ($4~"VIII" )print $0 } else print $0}' X152_not_diff_exp.txt > X152_not_diff_exp.COI1.txt

#find genes that are significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) { if ($4~"VIII")print $0 }else print $0}' X152_diff_exp.txt > X152_diff_exp.chrCOI1.txt


#find genes that are not significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) {  if ($4~"I" && $4!~"V" && $4!~"X" && $4!~"II" && $4!~"III") print $0 } else print $0}' X152_not_diff_exp.txt > X152_not_diff_exp.COI2.txt

#find genes that are significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) { if ($4~"I" && $4!~"V" && $4!~"X" && $4!~"II" && $4!~"III") print $0 }else print $0}' X152_diff_exp.txt > X152_diff_exp.chrCOI2.txt

######################################################
##add .txt to end of filename and also include line number
cp gene_exp.diff  X117_gene_exp.diff.txt

##find genes that are not differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12>="0.05") print $0 } else print $0}' X117_gene_exp.diff.txt > X117_not_diff_exp.txt

#find genes that are significantly differentially expressed
awk -F'\t' '{
if (NR>1) {
if ($12<="0.05") print $0 } else print $0}' X117_gene_exp.diff.txt > X117_diff_exp.txt

#find genes that are not significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) {  if ($7~"V" && $7!~"I" && $7!~"X") print $0 } else print $0}' X117_not_diff_exp.txt > X117_not_diff_exp.COI1.txt

#find genes that are significantly differentially expressed on the chromosome of interest
awk -F'\t' '{ if (NR>1) { if ($7~"V" && $7!~"I" && $7!~"X") print $0 }else print $0}' X117_diff_exp.txt > X117_diff_exp.chrCOI1.txt



($7~"V" && $7!~"I" && $7!~"X")


/panfs/pstor.storage/scratch/dwhlab/Holly/RNA_seq/fastq_files/cuffnorm_out_2/genes.count_table

cp genes.fpkm_table genes.fpkm_table.txt

cp genes.attr_table genes.attr_table.txt



##need to index bam files to view in IGV

/usr/local/samtools/latest/samtools index Sample${i}.sorted.bam


#################################################################################
#NEXT STEPS:
#want to look for genes that are known to be upregulated or downregulated during
#the ESR
#also want to look for histone genes
#so need to come up with pipeline to do that.

####################################################################
#sort the files by yes or no according to cuffdiff's idea of significant
awk -F'\t' '{
if (NR>1) {
if ($14=="yes") print $0 } else print $0}' gene*.txt > gene*_sig.txt
#above code doesn't work to do it to all files individually, so need to put it in a for loop

for file in ./*.txt
do
  FBASE=$(basename $file .txt)
  BASE=${FBASE%.txt}
  awk -F'\t' '{
  if (NR>1) {
  if ($14=="yes") print $0 } else print $0}' ${BASE}.txt > ${BASE}.txt
done


###merge HTseq data into one file
FILES=$(ls -t -v *.sam | tr '\n' ' ');awk 'NF > 1{ a[$1] = a[$1]"\t"$2} END {for( i in a ) print i a[i]}' $FILES > merged.sam
