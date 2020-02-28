# Dosage_Compensation


**Quality control:** FastQC version 1.8.0_20 using default parameters (available at https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 

**Low-quality bases trimmed:** Trimgalore v. 0.4.4 using -phred 33, -q 20 (available at https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) 

RNA samples were **aligned:* *Tophat v. 2.1.1 with -i 10 -I 10000 (TRAPNELL et al. 2012)

**annotated transcripts** (https://support.illumina.com/sequencing/sequencing_software/igenome.html)

**assemble the transcriptomes:** Cufflinks v. 2.2.1 default parameters (TRAPNELL et al. 2012). 

**Normalization:** Cuffnorm v. 2.2.1 with default parameters

**differential expression:** was found using Cuffdiff v. 2.2.1 with default parameters 

**read counts:** HTseq v. 0.6.1pl (Python v. 2.7.8) (ANDERS et al. 2015) 

**convert** .sam files into .bam files and **sort** the resulting .bam files:  Samtools (version 1.3.1) (LI et al. 2009) 
Scripts can be found at https://github.com/hollygene/Dosage_Compensation/blob/master/assembly_script.sh 
edgeR (ROBINSON et al. 2010) was used to calculate counts per million, then filtered based on high variance 

**FPKM:** Cuffnorm (TRAPNELL et al. 2012) 
A homebrew bash script was used to join the FPKM values for each line with the gene attributes file, turn the file into a .csv, remove mitochondrial sequences, and change the chromosome names from Roman numerals to numbers (script can be found at https://github.com/hollygene/Dosage_Compensation/blob/master/DC_workflow_April2017.sh)

R scripts are located at https://github.com/hollygene/Dosage_Compensation/blob/master/DC_workflow_old_MA.Rmd and https://github.com/hollygene/Dosage_Compensation/blob/master/DC_workflow.Rmd
  
_**Individual gene analysis**_

**raw count data:** HTseq 

**differential expression analysis:** DESeq2 default parameters, FDR of 0.1
https://github.com/hollygene/Dosage_Compensation/blob/master/R/scripts/R_scripts/DESeq2.Rmd

**Annotations:** TxDb S cerevisiae database (CARLSON M 2015)

**Histograms:** ggplot2 in the R package (WICKHAM 2016)
