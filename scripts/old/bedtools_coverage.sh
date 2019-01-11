#!/bin/bash
#PBS -t $'\t'N assembly_test_2
#PBS -t $'\t'q batch
#PBS -t $'\t'l nodes=1:ppn=4:HIGHMEM
#PBS -t $'\t'l walltime=480:00:00
#PBS -t $'\t'l mem=300gb
#PBS -t $'\t'M hmcqueary@uga.edu
#PBS -t $'\t'm ae

basedir="/lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly"
cd $basedir

module load bedtools/2.25.0

bedtools coverage -t $'\t'a /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/genes.gtf -t $'\t'b /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/accepted_hits.bam > /lustre1/hcm14449/SC_RNAseq/RNA_seq/November_2017_Assembly/GC/trimmed/tophat_update/Holly_GC_Anc_B_S2_R1_001_trimmed_tophat_out/GC_anc_B_counts.gtf

awk -t $'\t'F";" '$1=$1' OFS="\t" 7A_counts.txt > 7A_counts_sep.txt

awk -t $'\t'F"\t" '{for(i=1;i<=NF;i++){if ($i ~/gene_id/){print i}}}' 7A_counts_sep.txt > 7A_counts_shrt.txt

awk -t $'\t'F"\t" '{print $9, $10}' 7A_counts.txt > 7A_counts.red.txt

sed '//' 7A_counts.txt > 7A_counts.short.txt

#this one works
join -t $'\t' Holly_1_redo_R1_001_trimmed_tophat_out.stranded.txt Holly_1B_S6_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_1C_S3_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_2A_S1_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_2B_S7_R1_001_trimmed_tophat_out.stranded.txt |
join -t $'\t' - Holly_2C_S4_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_3A_S2_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_3B_redo_S12_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_3C_S5_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_4-redo_S3_R1_001_trimmed_tophat_out.stranded.txt |
join -t $'\t' - Holly_4B_S9_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_4C_S6_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_5A_S3_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_5B_S10_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_5C_redo_S16_R1_001_trimmed_tophat_out.stranded.txt |
join -t $'\t' - Holly_6A_S4_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_6B_S11_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_6C_S8_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_7-redo_S4_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_7B_S12_R1_001_trimmed_tophat_out.stranded.txt |
join -t $'\t' - Holly_7C_S9_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_8-redo_S5_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_8B_S13_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_8C_S10_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_9A_S5_R1_001_trimmed_tophat_out.stranded.txt |
join -t $'\t' - Holly_9B_S14_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_9C_S11_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_11_S17_R1_001_trimmed_tophat_out.stranded.txt |
join -t $'\t' - Holly_11B_S15_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_11C_S12_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_18-redo_S7_R1_001_trimmed_tophat_out.stranded.txt |
join -t $'\t' - Holly_18B_redo_S13_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_18C_S13_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_21-redo_S8_R1_001_trimmed_tophat_out.stranded.txt |
join -t $'\t' - Holly_21C_S14_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_31-redo_S9_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_31B_S19_R1_001_trimmed_tophat_out.stranded.txt |
join -t $'\t' - Holly_31C_S16_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_49-redo_S10_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_49B_S20_R1_001_trimmed_tophat_out.stranded.txt |
join -t $'\t' - Holly_49C_S17_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_59_S12_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_59B_S22_R1_001_trimmed_tophat_out.stranded.txt |
join -t $'\t' - Holly_59C_S19_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_61-redo_S11_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_61B_S23_R1_001_trimmed_tophat_out.stranded.txt |
join -t $'\t' - Holly_61C_S20_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_69_S1_R1_001_tophat_out.stranded.txt | join -t $'\t' - Holly_69B_S25_R1_001_trimmed_tophat_out.stranded.txt |
join -t $'\t' - Holly_69C_S22_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_76_S2_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_76B_S26_R1_001_trimmed_tophat_out.stranded.txt |
join -t $'\t' - Holly_76C_redo_S18_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_77_S3_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_77B_S27_R1_001_trimmed_tophat_out.stranded.txt |
join -t $'\t' - Holly_77C_S24_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_Anc_A_R1_001_trimmed_tophat_out.stranded.txt | join -t $'\t' - Holly_Anc_B_S2_R1_001_trimmed_tophat_out.stranded.txt |
join -t $'\t' - Holly_Anc_C_S32_R1_001_trimmed_tophat_out.stranded.txt > GC_counts.txt


function multijoin() {
    out=$1
    shift 1
    cat $1 | awk '{print $1}' > $out
    for f in $*; do join $out $f > tmp; mv tmp $out; done
}

multijoin GC_counts.txt Holly*
