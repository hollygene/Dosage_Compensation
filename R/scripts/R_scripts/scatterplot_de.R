## Code to make a scatterplot with two samples for one chromosome
# need to modify later to: 
  # use more samples
  # maybe for more than one chromosome
  

# put all sample expression data in the same data frame 
allSamplesExpr <- rbind(res1Dat, res2Dat, res3Dat, res4Dat, res5Dat, res7Dat, res8Dat, res9Dat, res11Dat, res18Dat, res49Dat, res59Dat, res61Dat, res76Dat, res77Dat)
# add chromosome numbers to this file
match <- intersect(allSamplesExpr$GENEID, genesNoAneu$GENEID)
# tell R which genes to keep 
keep <- genesNoAneu$GENEID %in% match
# subset list of genes in analysis to include only the ones matching those found in          sample file 
genes <- genesNoAneu[keep,]
# add the correct chromosomes in 
test <- cbind(genes$GENEID, genes$TXCHROM, allSamplesExpr)
chr1.all <- subset(allSamplesExpr, chr==1)
sample1 <- subset(allSamplesExpr, Sample==1)
View(sample1)
sample7 <- subset(allSamplesExpr, Sample==7)

{ plot(sample1$log2FoldChange,sample7$log2FoldChange)
  # want to draw two lines to indicate differential expression cutoffs
  # one for monosomy 
  abline(-1, 1, col='purple', lty="dashed")
  #one for trisomy 
  abline(0.583, 1, col="purple", lty="dashed") 
  # # one for tetrasomy 
  # abline(1, 1, col="purple", lty="dashed")
  # tell R to color dots above or below these lines certain colors 
  above <- sample7$log2FoldChange-sample1$log2FoldChange > 1
  points(sample7$log2FoldChange[above], sample1$log2FoldChange[above], col="blue", cex=0.5)
  below <- sample1$log2FoldChange-sample7$log2FoldChange < -1
  points(sample1$log2FoldChange[below], sample7$log2FoldChange[below], col="red", cex=0.5) }
