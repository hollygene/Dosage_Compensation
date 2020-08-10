#!/bin/bash
FILES=/Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Indiv_Genes/Cuffdiff_GC_MA/aneuploid/GC/*.txt
for f in $FILES;
do
  awk -F'\t' '{
  if (NR>1) {
  if ($14=="yes") print $0 } else print $0}' "$f" > "$f"_significant.txt;
done

#finding same genes DE based on "gene" column
awk 'NR==FNR {a[$3]; next} $3 in a {print}' gene_exp_S11_GCcopycopy2.txt_significant.txt gene_exp_S18_GCcopycopy2.txt_significant.txt  > 18_11_sig.txt

awk 'NR==FNR {a[$3]; next} $3 in a {print}' gene_exp_S21_GCcopycopy2.txt_significant.txt 18_11_sig.txt  > 18_11_21_sig.txt

awk 'NR==FNR {a[$3]; next} $3 in a {print}' gene_exp_S49_GCcopycopy2.txt_significant.txt 18_11_21_sig.txt > 18_11_21_49_sig.txt

awk 'NR==FNR {a[$3]; next} $3 in a {print}' gene_exp_S4_GCcopycopy2.txt_significant.txt 18_11_21_49_sig.txt > 18_11_21_49_4_sig.txt

awk 'NR==FNR {a[$3]; next} $3 in a {print}' gene_exp_S59_GCcopycopy2.txt_significant.txt 18_11_21_49_4_sig.txt > 18_11_21_49_4_59_sig.txt

awk 'NR==FNR {a[$3]; next} $3 in a {print}' gene_exp_S61_GCcopycopy2.txt_significant.txt 18_11_21_49_4_59_sig.txt > 18_11_21_49_4_59_61_sig.txt

awk 'NR==FNR {a[$3]; next} $3 in a {print}' gene_exp_S66_GCcopycopy2.txt_significant.txt 18_11_21_49_4_59_61_sig.txt > 18_11_21_49_4_59_61_66_sig.txt

awk 'NR==FNR {a[$3]; next} $3 in a {print}' gene_exp_S76_GCcopycopy2.txt_significant.txt 18_11_21_49_4_59_61_66_sig.txt > 18_11_21_49_4_59_61_66_76_sig.txt
