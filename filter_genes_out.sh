# need to match the first column of the counts file to the first column of the genes file to filter out genes that I don't want
# best to do this in a loop since I need to do it for all of my samples
# then once I have files with only the genes I want, can put them into DESeq2/edgeR and do analyses


#loction of count files
/Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Indiv_Genes/HTseq/HTseq_update_May2018/GC/exons_unstranded

#location of gene files
/Users/hollymcqueary/Documents/GitHub/Dosage_Compensation/R/scripts/

# to remove quotation marks
sed 's/"//g' #-- removes all quotation marks on each line

#change from comma separated to tab-delimited
 < genesGC.fil.2.csv  tr "," "\t" > genesGC.fil.2.txt

# compare the two files and print the lines of the second that match the first column of the first files
awk -F'|' 'NR==FNR{c[$1]++;next};c[$1] > 0' file1 file2 > output.txt


#remove quotation marks from the gene lists

sed 's/"//g' genesGC.fil.csv > genesGC.fil.2.csv

sed 's/"//g'  genesMANew.fil.csv > genesMANew.fil.2.csv

sed 's/"//g'  genesMAOld.fil.csv > genesMAOld.fil.2.csv

#change gene files from comma separated to tab delimited
 < genesGC.fil.2.csv  tr "," "\t" > genesGC.fil.2.txt
 < genesMANew.fil.2.csv  tr "," "\t" > genesMANew.fil.2.txt
  < genesMAOld.fil.2.csv  tr "," "\t" > genesMAOld.fil.2.txt

  # compare the two files and print the lines of the second that match the first column of the first files
  for file in /Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Indiv_Genes/HTseq/HTseq_update_May2018/GC/exons_unstranded/*.txt

  do

  FBASE=$(basename $file .txt)

  BASE=${FBASE%.txt}

  awk -F'\t' 'NR==FNR{c[$2]++;next};c[$1] > 0' /Users/hollymcqueary/Documents/GitHub/Dosage_Compensation/R/scripts/genesFil/genesGC.fil.2.txt /Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Indiv_Genes/HTseq/HTseq_update_May2018/GC/exons_unstranded/${BASE}.txt > /Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Indiv_Genes/HTseq/HTseq_update_May2018/GC/exons_unstranded/filtered_genes/${BASE}.txt

  done

#########################################################
  for file in /Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Indiv_Genes/HTseq/HTseq_update_May2018/MA_new/exons_unstranded/*.txt

  do

  FBASE=$(basename $file .txt)

  BASE=${FBASE%.txt}

  awk -F'\t' 'NR==FNR{c[$2]++;next};c[$1] > 0' /Users/hollymcqueary/Documents/GitHub/Dosage_Compensation/R/scripts/genesFil/genesMANew.fil.2.txt /Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Indiv_Genes/HTseq/HTseq_update_May2018/MA_new/exons_unstranded/${BASE}.txt > /Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Indiv_Genes/HTseq/HTseq_update_May2018/MA_new/exons_unstranded/filtered_genes/${BASE}.txt

  done

  #########################################################
    for file in /Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Indiv_Genes/HTseq/HTseq_update_May2018/MA_old/*.txt

    do

    FBASE=$(basename $file .txt)

    BASE=${FBASE%.txt}

    awk -F'\t' 'NR==FNR{c[$2]++;next};c[$1] > 0' /Users/hollymcqueary/Documents/GitHub/Dosage_Compensation/R/scripts/genesFil/genesMAOld.fil.2.txt /Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Indiv_Genes/HTseq/HTseq_update_May2018/MA_old/${BASE}.txt > /Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Indiv_Genes/HTseq/HTseq_update_May2018/MA_old/filtered_genes/${BASE}.txt

    done
