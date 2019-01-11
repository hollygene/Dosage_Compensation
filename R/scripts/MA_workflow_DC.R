#R codes 
#New for MA data only 
#Starting it Sept 2017

#first I need to load in the required packages
library(moments)
library("Hmisc")
library(ggplot2)

#then I need to load in my data
sc_counts <- read.csv("/Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Redo/September_2017/Sept25/genes.fpkm_table.joined.rmcom.split.commas.nomito.chgchr.rmhead.nwhd.MA.csv",header=T)
#check it out 
class(sc_counts)
str(sc_counts)
#looks good

#to prevent problems, add .5 to entire dataset so where it is =0 it doesn't cause issues later 
sc_counts.z <- sc_counts
sc_counts.z[,6:32] <- sc_counts.z[,6:32]+.5 #add .5 to dataset so the zeroes dont mess up

gene.means <- rowMeans(sc_counts.z[,6:29]) #average across the sample
anc.means <- rowMeans(sc_counts.z[,30:32]) #average in the ancestor

#remove lowly-expressed genes
#if the average fpkm is less than 5 in the ancestor or less than 5 in the overall, remove the gene
sc_counts.final <- sc_counts.z[-which(gene.means<5 | anc.means<5),]

#final matrix for ancestor and descendans
sc_anc <- sc_counts.final[,c(1:5,30,31,32)]
sc_des <- sc_counts.final[,c(1:5,6:29)] #average for each GC line

#the purpose of these is to get the rows of each chromosome

index1 <- which(sc_anc[,1]==1)
index2 <- which(sc_anc[,1]==2)
index3 <- which(sc_anc[,1]==3)
index4 <- which(sc_anc[,1]==4)
index5 <- which(sc_anc[,1]==5)
index6 <- which(sc_anc[,1]==6)
index7 <- which(sc_anc[,1]==7)
index8 <- which(sc_anc[,1]==8)
index9 <- which(sc_anc[,1]==9)
index10 <- which(sc_anc[,1]==10)
index11 <- which(sc_anc[,1]==11)
index12 <- which(sc_anc[,1]==12)
index13 <- which(sc_anc[,1]==13)
index14 <- which(sc_anc[,1]==14)
index15 <- which(sc_anc[,1]==15)
index16 <- which(sc_anc[,1]==16)

# FOLD CHANGES: average of line over average of ancestor
###ancestor average
ancAvg <- rowMeans(sc_anc[,6:8],na.rm=T)
ancAvg <- cbind(sc_anc[,1:2],ancAvg)

#first write some functions that help automate the process
#this gets the average FPKM for a given line from the sc_des matrix
getAvgFPKM.MA <- function(line) { #for line as "###" string
  r1 <- paste(c("X",line,"_0"),collapse="") #triplicates are named X[LINE]_[rep]
  r2 <- paste(c("X",line,"_1"),collapse="")
  r3 <- paste(c("X",line,"_2"),collapse="")
  reps <- sc_des[,c("chr","tracking_id",r1,r2,r3)] #get those 4 columns out of the descendant dataset
  reps[reps==0] <- NA #replace 0 values with NA as before
  desAvg.MA <- rowMeans(reps[,c(-1,-2)],na.rm=T)
  desAvg.MA <- cbind(sc_des[,c(1,2)],desAvg.MA)
  #colnames(desAvg) <- c("chr","gene",paste(c(line,"_Avg"),collapse=""))
  return(desAvg.MA) #returned is a matrix equivalent to ancAvg: chrm, gene, and line average FPKM 
}


#this gets the fold change for the whole line
getFPKMRatio.MA <- function(line) { #line must be a 3 digit character
  desAvg.MA <- getAvgFPKM.MA(line) #first call avg fpkm script, which returns a 3 column matrix
  ratioByGene.MA <- desAvg.MA[,3]/ancAvg[,3]
  rat.MA <- cbind(sc_des[,1:2],ratioByGene.MA) #add gene and chrm values
  colnames(rat.MA) <- c("chr","tracking_id","ratio")
  return(rat.MA)
}

#this gets the fold change for just the specified chromosome
getChrmRatio.MA <- function(line, chr){ #input line as "XXX", chr as a number
  rat.MA <- getFPKMRatio.MA(line) #get list of ratios which includes chrm and gene
  rat.c.MA <- rat.MA[(rat.MA[,1]==chr),] #get just the rows for the appropriate chrom
  rat.cs.MA <- rat.c.MA[order(rat.c.MA[,3]),]#order by column 3
  return(rat.cs.MA)
}

#now we can look at the fold changes for each chromosome and each line
#make a matrix of median ratio:
y.MA=NULL
rep.MA <- c("112","115","117","123","141","152","29","50")
#View(rep.MA)
for(line in rep.MA){
  for(chr in 1:16){
    y.MA <- c(y.MA,median(getChrmRatio.MA(line,chr)[,3],na.rm=T))
  }
}

ratMatMed.MA <- matrix(y.MA,nrow=16,ncol=8)
colnames(ratMatMed.MA) <- rep.MA
ratMatMed.MA <- cbind(c(1:16),ratMatMed.MA)
colnames(ratMatMed.MA)[1] <- "chr"

#make a matrix that is the LOG2 RATIO
yl.MA=NULL
for(line in rep.MA){
  for(chrm in 1:16){
    yl.MA <- c(yl.MA,median(log2(getChrmRatio.MA(line,chrm)[,3])))
  }
}

warnings()
logMat.MA <- matrix(yl.MA,nrow=16,ncol=8)
colnames(logMat.MA) <- rep.MA
logMat.MA <- data.frame(logMat.MA)
chrms <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
logMat2.MA <-cbind(chrms,logMat.MA)

##make a matrix of the standard deviations for each line and each chromosome
ysd.MA=NULL
for(line in rep.MA){
  for(chrm in 1:16){
    ysd.MA <- c(ysd.MA,sd(log2(getChrmRatio.MA(line,chrm)[,3])))
  }
}

logMatsd.MA <- matrix(ysd.MA,nrow=16,ncol=8)
colnames(logMatsd.MA) <- rep.MA
logMatsd.MA <- data.frame(logMatsd.MA)
chrms <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
logMatsd2.MA <-cbind(chrms,logMatsd.MA)
write.csv(logMat2.MA,file="logMat.MA")
shapiro.test(logMat2.MA)

######################################################################################
#how to make a graph of the chromosomes and their respective error bars 
#changes the margins to the standard margin width
par(mar=c(5, 4, 4, 2) + 0.1)
par(mfrow=c(1,1))



#line 112
#get the standard deviations for each chromosome for line 112 only
X112 <- subset(logMat.MA, select=c("X112"))
X112.sd <- subset(logMatsd2.MA, select=c("X112"))
#merge the stdev with the mean and the chromosome number
X112.all <- cbind(chrms,X112,X112.sd)
#change the column names appropriately
colnames(X112.all)[1] <- "chr"; colnames(X112.all)[2] <- "mean";colnames(X112.all)[3] <- "sd"

#plots the chromosomes with error bars for each one 
with (
  data = X112.all
  , expr = errbar(main="MA Line 112",chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a title to the graph
title(main="MA Line 112")
abline(h=0,lty=3)
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 115
par(mar=c(5, 4, 4, 2) + 0.1)
#get the standard deviations for each chromosome for line 115 only
X115 <- subset(logMat.MA, select=c("X115"))
X115.sd <- subset(logMatsd2.MA, select=c("X115"))
#merge the stdev with the mean and the chromosome number
X115.all <- cbind(chrms,X115,X115.sd)
#change the column names appropriately
colnames(X115.all)[1] <- "chr"; colnames(X115.all)[2] <- "mean";colnames(X115.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X115.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="MA Line 115")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 117
#get the standard deviations for each chromosome for line 117 only
X117 <- subset(logMat.MA, select=c("X117"))
X117.sd <- subset(logMatsd2.MA, select=c("X117"))
#merge the stdev with the mean and the chromosome number
X117.all <- cbind(chrms,X117,X117.sd)
#change the column names appropriately
colnames(X117.all)[1] <- "chr"; colnames(X117.all)[2] <- "mean";colnames(X117.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X117.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="MA Line 117")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 123
X123 <- subset(logMat.MA, select=c("X123"))
X123.sd <- subset(logMatsd2.MA, select=c("X123"))
#merge the stdev with the mean and the chromosome number
X123.all <- cbind(chrms,X123,X123.sd)
#change the column names appropriately
colnames(X123.all)[1] <- "chr"; colnames(X123.all)[2] <- "mean";colnames(X123.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X123.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="MA Line 123")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 141
X141 <- subset(logMat.MA, select=c("X141"))
X141.sd <- subset(logMatsd2.MA, select=c("X141"))
#merge the stdev with the mean and the chromosome number
X141.all <- cbind(chrms,X141,X141.sd)
#change the column names appropriately
colnames(X141.all)[1] <- "chr"; colnames(X141.all)[2] <- "mean";colnames(X141.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X141.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="MA Line 141")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 152
X152 <- subset(logMat.MA, select=c("X152"))
X152.sd <- subset(logMatsd2.MA, select=c("X152"))
#merge the stdev with the mean and the chromosome number
X152.all <- cbind(chrms,X152,X152.sd)
#change the column names appropriately
colnames(X152.all)[1] <- "chr"; colnames(X152.all)[2] <- "mean";colnames(X152.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X152.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="MA Line 152")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 29
X29 <- subset(logMat.MA, select=c("X29"))
X29.sd <- subset(logMatsd2.MA, select=c("X29"))
#merge the stdev with the mean and the chromosome number
X29.all <- cbind(chrms,X29,X29.sd)
#change the column names appropriately
colnames(X29.all)[1] <- "chr"; colnames(X29.all)[2] <- "mean";colnames(X29.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X29.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="MA Line 29")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 50
X50 <- subset(logMat.MA, select=c("X50"))
X50.sd <- subset(logMatsd2.MA, select=c("X50"))
#merge the stdev with the mean and the chromosome number
X50.all <- cbind(chrms,X50,X50.sd)
#change the column names appropriately
colnames(X50.all)[1] <- "chr"; colnames(X50.all)[2] <- "mean";colnames(X50.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X50.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="MA Line 50")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

##############################################################################################
l112.1 <- getChrmRatio.MA("112",1)
l115.1 <- getChrmRatio.MA("115",1)
l117.1 <- getChrmRatio.MA("117",1)
l123.1 <- getChrmRatio.MA("123",1)
l141.1 <- getChrmRatio.MA("141",1)
l152.1 <- getChrmRatio.MA("152",1)
l29.1 <- getChrmRatio.MA("29",1)
l50.1 <- getChrmRatio.MA("50",1)

l112.5 <- getChrmRatio.MA("112",5)
l115.5 <- getChrmRatio.MA("115",5)
l117.5 <- getChrmRatio.MA("117",5)
l123.5 <- getChrmRatio.MA("123",5)
l141.5 <- getChrmRatio.MA("141",5)
l152.5 <- getChrmRatio.MA("152",5)
l29.5 <- getChrmRatio.MA("29",5)
l50.5 <- getChrmRatio.MA("50",5)

l112.7 <- getChrmRatio.MA("112",7)
l115.7 <- getChrmRatio.MA("115",7)
l117.7 <- getChrmRatio.MA("117",7)
l123.7 <- getChrmRatio.MA("123",7)
l141.7 <- getChrmRatio.MA("141",7)
l152.7 <- getChrmRatio.MA("152",7)
l29.7 <- getChrmRatio.MA("29",7)
l50.7 <- getChrmRatio.MA("50",7)

l112.8 <- getChrmRatio.MA("112",8)
l115.8 <- getChrmRatio.MA("115",8)
l117.8 <- getChrmRatio.MA("117",8)
l123.8 <- getChrmRatio.MA("123",8)
l141.8 <- getChrmRatio.MA("141",8)
l152.8 <- getChrmRatio.MA("152",8)
l29.8 <- getChrmRatio.MA("29",8)
l50.8 <- getChrmRatio.MA("50",8)

l112.9 <- getChrmRatio.MA("112",9)
l115.9 <- getChrmRatio.MA("115",9)
l117.9 <- getChrmRatio.MA("117",9)
l123.9 <- getChrmRatio.MA("123",9)
l141.9 <- getChrmRatio.MA("141",9)
l152.9 <- getChrmRatio.MA("152",9)
l29.9 <- getChrmRatio.MA("29",9)
l50.9 <- getChrmRatio.MA("50",9)

l112.10 <- getChrmRatio.MA("112",10)
l115.10 <- getChrmRatio.MA("115",10)
l117.10 <- getChrmRatio.MA("117",10)
l123.10 <- getChrmRatio.MA("123",10)
l141.10 <- getChrmRatio.MA("141",10)
l152.10 <- getChrmRatio.MA("152",10)
l29.10 <- getChrmRatio.MA("29",10)
l50.10 <- getChrmRatio.MA("50",10)

l112.12 <- getChrmRatio.MA("112",12)
l115.12 <- getChrmRatio.MA("115",12)
l117.12 <- getChrmRatio.MA("117",12)
l123.12 <- getChrmRatio.MA("123",12)
l141.12 <- getChrmRatio.MA("141",12)
l152.12 <- getChrmRatio.MA("152",12)
l29.12 <- getChrmRatio.MA("29",12)
l50.12 <- getChrmRatio.MA("50",12)

l112.14 <- getChrmRatio.MA("112",14)
l115.14 <- getChrmRatio.MA("115",14)
l117.14 <- getChrmRatio.MA("117",14)
l123.14 <- getChrmRatio.MA("123",14)
l141.14 <- getChrmRatio.MA("141",14)
l152.14 <- getChrmRatio.MA("152",14)
l29.14 <- getChrmRatio.MA("29",14)
l50.14 <- getChrmRatio.MA("50",14)

l112.15 <- getChrmRatio.MA("112",15)
l115.15 <- getChrmRatio.MA("115",15)
l117.15 <- getChrmRatio.MA("117",15)
l123.15 <- getChrmRatio.MA("123",15)
l141.15 <- getChrmRatio.MA("141",15)
l152.15 <- getChrmRatio.MA("152",15)
l29.15 <- getChrmRatio.MA("29",15)
l50.15 <- getChrmRatio.MA("50",15)

l112.16 <- getChrmRatio.MA("112",16)
l115.16 <- getChrmRatio.MA("115",16)
l117.16 <- getChrmRatio.MA("117",16)
l123.16 <- getChrmRatio.MA("123",16)
l141.16 <- getChrmRatio.MA("141",16)
l152.16 <- getChrmRatio.MA("152",16)
l29.16 <- getChrmRatio.MA("29",16)
l50.16 <- getChrmRatio.MA("50",16)

#BOXPLOTS PER LINE AND CHROMOSOME 
#CHROMOSOME 1 

#because for this experiment I don't have any euploid replicates (without using the old data), I have 
#to use lines that are euploid at that specific chromosome
#line 152 is 3n for c1
par(mfrow=c(1,1))
par(mar = c(5, 4, 4, 2) + 0.1)

l152.1 <- getChrmRatio.MA("152",1)
euRatio.c1 <- c(l112.1[,3],l115.1[,3],l117.1[,3],l123.1[,3],l141.1[,3],l29.1[,3],l50.1[,3])

boxplot((log2(l152.1[,3])),(log2(euRatio.c1)),names=c("Trisomic","Disomic"),ylab="log2(fold change)", col=c("purple","turquoise"), main="Chr 1 Line 152 vs Avg Euploid Lines")
abline(h=0,lty=3)
#run a t-test between the two 
wilcox.test((log2(l152.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
p <- "< 2.2e-16"         
mylabel = bquote(italic(p) == .(format(p, digits = 9)))
text(x=2.3,y = 1.6, labels = mylabel)

#CHROMOSOME 5
#line 50 and line 117 are 3n for chromosome 5 
euRatio.c5 <- c(l112.5[,3],l115.5[,3],l152.5[,3],l123.5[,3],l141.5[,3],l29.5[,3])
#boxplot
boxplot((log2(l117.5[,3])),(log2(l50.5[,3])),(log2(euRatio.c5)),names=c("Line 117","Line 50","Disomic"),ylab="log2(fold change)", col=c("purple","turquoise","lawngreen"), main="Chr 5 Line 117 and 50 vs Euploid Lines")
abline(h=0,lty=3)
#run a t-test between the two 
wilcox.test((log2(l117.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
p <- "< 2.2e-16"         
mylabel = bquote(italic(p.117) == .(format(p, digits = 9)))
text(x=3.3,y = 4, labels = mylabel)
#between 50 and the euploids 
wilcox.test((log2(l50.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
p <- "0.01778"         
mylabel = bquote(italic(p.50) == .(format(p, digits = 9)))
text(x=3.3,y = 3.6, labels = mylabel)
##even though line 50 doesn't look trisomic, the p-value is still significant 

#CHROMOSOME 7 
#line 115 is trisomic for chromosome 7
euRatio.c7 <- c(l112.7[,3],l152.7[,3],l117.7[,3],l123.7[,3],l141.7[,3],l29.7[,3],l50.7[,3])

boxplot((log2(l115.7[,3])),(log2(euRatio.c7)),names=c("Trisomic","Disomic"),ylab="log2(fold change)", col=c("purple","turquoise"), main="Chr 7 Line 115 vs Avg Euploid Lines")
abline(h=0,lty=3)
#run a t-test between the two 
wilcox.test((log2(l115.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
p <- "0.606"         
mylabel = bquote(italic(p) == .(format(p, digits = 9)))
text(x=2.3,y = 3.7, labels = mylabel)

#CHROMOSOME 8
#line 152 is trisomic for chromosome 8
euRatio.c8 <- c(l112.8[,3],l115.8[,3],l117.8[,3],l123.8[,3],l141.8[,3],l29.8[,3],l50.8[,3])

boxplot((log2(l152.8[,3])),(log2(euRatio.c8)),names=c("Trisomic","Disomic"),ylab="log2(fold change)", col=c("purple","turquoise"), main="Chr 8 Line 152 vs Avg Euploid Lines")
abline(h=0,lty=3)
#run a t-test between the two 
wilcox.test((log2(l152.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
p <- "< 2.2e-16"         
mylabel = bquote(italic(p) == .(format(p, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel)

#CHROMOSOME 9
#line 29 is monosomic for chromosome 9
euRatio.c9 <- c(l112.9[,3],l115.9[,3],l117.9[,3],l123.9[,3],l141.9[,3],l152.9[,3],l50.9[,3])

boxplot((log2(l29.9[,3])),(log2(euRatio.c9)),names=c("Monosomic","Disomic"),ylab="log2(fold change)", col=c("purple","turquoise"), main="Chr 9 Line 29 vs Avg Euploid Lines")
abline(h=0,lty=3)
#run a t-test between the two 
wilcox.test((log2(l29.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
p <- "< 2.2e-16"         
mylabel = bquote(italic(p) == .(format(p, digits = 9)))
text(x=2.3,y = 6, labels = mylabel)

#CHROMOSOME 12
#line 123 is trisomic for chromosome 12
euRatio.c12 <- c(l112.12[,3],l115.12[,3],l117.12[,3],l29.12[,3],l141.12[,3],l152.12[,3],l50.12[,3])

boxplot((log2(l123.12[,3])),(log2(euRatio.c12)),names=c("Trisomic","Disomic"),ylab="log2(fold change)", col=c("purple","turquoise"), main="Chr 12 Line 123 vs Avg Euploid Lines")
abline(h=0,lty=3)
#run a t-test between the two 
wilcox.test((log2(l123.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
p <- "0.0158"         
mylabel = bquote(italic(p) == .(format(p, digits = 9)))
text(x=2.3,y = 4, labels = mylabel)

#CHROMOSOME 16
#line 112 and line 141 are trisomic for chromosome 16
euRatio.c16 <- c(l123.16[,3],l115.16[,3],l117.16[,3],l29.16[,3],l152.16[,3],l50.16[,3])

boxplot((log2(l141.16[,3])),(log2(l112.16[,3])),(log2(euRatio.c16)),names=c("141","112","Disomic"),ylab="log2(fold change)", col=c("purple","turquoise","lawngreen"), main="Chr 16 Trisomic vs Avg Euploid Lines")
abline(h=0,lty=3)
#run a t-test between the two 
wilcox.test((log2(l112.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
p <- "< 2.2e-16"         
mylabel = bquote(italic(p.112) == .(format(p, digits = 9)))
text(x=3.3,y = 4, labels = mylabel)
#and for 141
wilcox.test((log2(l141.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
p <- "< 2.2e-16"         
mylabel = bquote(italic(p.141) == .(format(p, digits = 9)))
text(x=3.3,y = 4.5, labels = mylabel)

##############################################################################################

#make a histogram of log ratio of all genes across 
#euploid lines
##for chromosome one 
##can't do this for euploid ratios
#has to be individual euploid lines 

#going to do this for each line at each chromosome that is aneuploid
#CHROMOSOME 1
par(mfrow=c(4,2))
par(mar = rep(2, 4))
#line 112
hist(log2(l112.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 1 Line 112")
lines(density(log2(l112.1[,3])),col="green")
mc11 <- mean(log2(l112.1[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l112.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l112.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#line 115
hist(log2(l115.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 1 Line 115")
lines(density(log2(l115.1[,3])),col="green")
mc11 <- mean(log2(l115.1[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l115.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l115.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#line 115
hist(log2(l117.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 1 Line 117")
lines(density(log2(l117.1[,3])),col="green")
mc11 <- mean(log2(l117.1[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l117.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l117.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#line 123
hist(log2(l123.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 1 Line 123")
lines(density(log2(l123.1[,3])),col="green")
mc11 <- mean(log2(l123.1[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l123.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l123.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#line 141
hist(log2(l141.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 1 Line 141")
lines(density(log2(l141.1[,3])),col="green")
mc11 <- mean(log2(l141.1[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l141.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l141.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#line 29
hist(log2(l29.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 1 Line 29")
lines(density(log2(l29.1[,3])),col="green")
mc11 <- mean(log2(l29.1[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l29.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l29.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#line 50
hist(log2(l50.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 1 Line 50")
lines(density(log2(l50.1[,3])),col="green")
mc11 <- mean(log2(l50.1[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l50.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l50.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#line 152
hist(log2(l152.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Trisomic Chr 1 Line 152")
lines(density(log2(l152.1[,3])),col="green")
mc11 <- mean(log2(l152.1[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l152.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l152.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

##CHROMOSOME 5 
#line 112
hist(log2(l112.5[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 5 Line 112")
lines(density(log2(l112.5[,3])),col="green")
mc11 <- mean(log2(l112.5[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l112.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l112.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#line 115
hist(log2(l115.5[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 5 Line 115")
lines(density(log2(l115.5[,3])),col="green")
mc11 <- mean(log2(l115.5[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l115.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l115.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .6, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .65, labels = mylabel3,cex=.8)

#line 117
hist(log2(l117.5[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Trisomic Chr 5 Line 117")
lines(density(log2(l117.5[,3])),col="green")
mc11 <- mean(log2(l117.5[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l117.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .3, labels = mylabel,cex=.8)
sd2 <- sd(log2(l117.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .35, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .4, labels = mylabel3,cex=.8)

#line 123
hist(log2(l123.5[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 5 Line 123")
lines(density(log2(l123.5[,3])),col="green")
mc11 <- mean(log2(l123.5[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l123.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l123.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#line 141
hist(log2(l141.5[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 5 Line 141")
lines(density(log2(l141.5[,3])),col="green")
mc11 <- mean(log2(l141.5[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l141.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l141.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#line 29
hist(log2(l29.5[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 5 Line 29")
lines(density(log2(l29.5[,3])),col="green")
mc11 <- mean(log2(l29.5[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l29.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .3, labels = mylabel,cex=.8)
sd2 <- sd(log2(l29.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .5, labels = mylabel3,cex=.8)

#line 50
hist(log2(l50.5[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 5 Line 50")
lines(density(log2(l50.5[,3])),col="green")
mc11 <- mean(log2(l50.5[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l50.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .65, labels = mylabel,cex=.8)
sd2 <- sd(log2(l50.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .75, labels = mylabel3,cex=.8)

#line 152
hist(log2(l152.5[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 5 Line 152")
lines(density(log2(l152.5[,3])),col="green")
mc11 <- mean(log2(l152.5[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l152.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l152.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#CHROMOSOME 7
#line 112
hist(log2(l112.7[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 7 Line 112")
lines(density(log2(l112.7[,3])),col="green")
mc11 <- mean(log2(l112.7[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l112.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .65, labels = mylabel,cex=.8)
sd2 <- sd(log2(l112.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .75, labels = mylabel3,cex=.8)

#line 115
hist(log2(l115.7[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Trisomic Chr 7 Line 115")
lines(density(log2(l115.7[,3])),col="green")
mc11 <- mean(log2(l115.7[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l115.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .65, labels = mylabel,cex=.8)
sd2 <- sd(log2(l115.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .75, labels = mylabel3,cex=.8)

#line 117
hist(log2(l117.7[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 7 Line 117")
lines(density(log2(l117.7[,3])),col="green")
mc11 <- mean(log2(l117.7[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l117.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .3, labels = mylabel,cex=.8)
sd2 <- sd(log2(l117.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .5, labels = mylabel3,cex=.8)

#line 123
hist(log2(l123.7[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 7 Line 123")
lines(density(log2(l123.7[,3])),col="green")
mc11 <- mean(log2(l123.7[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l123.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .5, labels = mylabel,cex=.8)
sd2 <- sd(log2(l123.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .6, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .7, labels = mylabel3,cex=.8)

#line 141
hist(log2(l141.7[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 7 Line 141")
lines(density(log2(l141.7[,3])),col="green")
mc11 <- mean(log2(l141.7[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l141.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l141.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#line 29
hist(log2(l29.7[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 7 Line 29")
lines(density(log2(l29.7[,3])),col="green")
mc11 <- mean(log2(l29.7[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l29.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .5, labels = mylabel,cex=.8)
sd2 <- sd(log2(l29.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .6, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .7, labels = mylabel3,cex=.8)

#line 50
hist(log2(l50.7[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 7 Line 50")
lines(density(log2(l50.7[,3])),col="green")
mc11 <- mean(log2(l50.7[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l50.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .3, labels = mylabel,cex=.8)
sd2 <- sd(log2(l50.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .5, labels = mylabel3,cex=.8)

#line 152
hist(log2(l152.7[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 7 Line 152")
lines(density(log2(l152.7[,3])),col="green")
mc11 <- mean(log2(l152.7[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l152.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l152.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#CHROMOSOME 8 
#line 112
hist(log2(l112.8[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 8 Line 112")
lines(density(log2(l112.8[,3])),col="green")
mc11 <- mean(log2(l112.8[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l112.8[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .5, labels = mylabel,cex=.8)
sd2 <- sd(log2(l112.8[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .6, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .7, labels = mylabel3,cex=.8)

#line 115
hist(log2(l115.8[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 8 Line 115")
lines(density(log2(l115.8[,3])),col="green")
mc11 <- mean(log2(l115.8[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l115.8[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .5, labels = mylabel,cex=.8)
sd2 <- sd(log2(l115.8[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .6, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .7, labels = mylabel3,cex=.8)

#line 115
hist(log2(l117.8[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 8 Line 117")
lines(density(log2(l117.8[,3])),col="green")
mc11 <- mean(log2(l117.8[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l117.8[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .3, labels = mylabel,cex=.8)
sd2 <- sd(log2(l117.8[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .5, labels = mylabel3,cex=.8)

#line 123
hist(log2(l123.8[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 8 Line 123")
lines(density(log2(l123.8[,3])),col="green")
mc11 <- mean(log2(l123.8[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l123.8[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .5, labels = mylabel,cex=.8)
sd2 <- sd(log2(l123.8[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .6, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .7, labels = mylabel3,cex=.8)

#line 141
hist(log2(l141.8[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 8 Line 141")
lines(density(log2(l141.8[,3])),col="green")
mc11 <- mean(log2(l141.8[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l141.8[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l141.8[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#line 29
hist(log2(l29.8[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 8 Line 29")
lines(density(log2(l29.8[,3])),col="green")
mc11 <- mean(log2(l29.8[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l29.8[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .5, labels = mylabel,cex=.8)
sd2 <- sd(log2(l29.8[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .6, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .7, labels = mylabel3,cex=.8)

#line 50
hist(log2(l50.8[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 8 Line 50")
lines(density(log2(l50.8[,3])),col="green")
mc11 <- mean(log2(l50.8[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l50.8[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .3, labels = mylabel,cex=.8)
sd2 <- sd(log2(l50.8[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .5, labels = mylabel3,cex=.8)

#line 152
hist(log2(l152.8[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Trisomic Chr 8 Line 152")
lines(density(log2(l152.8[,3])),col="green")
mc11 <- mean(log2(l152.8[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l152.8[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l152.8[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#CHROMOSOME 9
#line 112
hist(log2(l112.9[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 9 Line 112")
lines(density(log2(l112.9[,3])),col="green")
mc11 <- mean(log2(l112.9[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l112.9[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l112.9[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#line 115
hist(log2(l115.9[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 9 Line 115")
lines(density(log2(l115.9[,3])),col="green")
mc11 <- mean(log2(l115.9[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l115.9[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .5, labels = mylabel,cex=.8)
sd2 <- sd(log2(l115.9[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .6, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .7, labels = mylabel3,cex=.8)

#line 115
hist(log2(l117.9[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 9 Line 117")
lines(density(log2(l117.9[,3])),col="green")
mc11 <- mean(log2(l117.9[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l117.9[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .5, labels = mylabel,cex=.8)
sd2 <- sd(log2(l117.9[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .6, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .7, labels = mylabel3,cex=.8)

#line 123
hist(log2(l123.9[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 9 Line 123")
lines(density(log2(l123.9[,3])),col="green")
mc11 <- mean(log2(l123.9[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l123.9[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .3, labels = mylabel,cex=.8)
sd2 <- sd(log2(l123.9[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .5, labels = mylabel3,cex=.8)

#line 141
hist(log2(l141.9[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 9 Line 141")
lines(density(log2(l141.9[,3])),col="green")
mc11 <- mean(log2(l141.9[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l141.9[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l141.9[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#line 29
hist(log2(l29.9[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Monosomic Chr 9 Line 29")
lines(density(log2(l29.9[,3])),col="green")
mc11 <- mean(log2(l29.9[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l29.9[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .3, labels = mylabel,cex=.8)
sd2 <- sd(log2(l29.9[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .5, labels = mylabel3,cex=.8)

#line 50
hist(log2(l50.9[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 9 Line 50")
lines(density(log2(l50.9[,3])),col="green")
mc11 <- mean(log2(l50.9[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l50.9[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel,cex=.8)
sd2 <- sd(log2(l50.9[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .5, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .6, labels = mylabel3,cex=.8)

#line 152
hist(log2(l152.9[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 9 Line 152")
lines(density(log2(l152.9[,3])),col="green")
mc11 <- mean(log2(l152.9[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l152.9[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l152.9[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#CHROMOSOME 12
#line 112
hist(log2(l112.12[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 12 Line 112")
lines(density(log2(l112.12[,3])),col="green")
mc11 <- mean(log2(l112.12[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l112.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel,cex=.8)
sd2 <- sd(log2(l112.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .5, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .6, labels = mylabel3,cex=.8)

#line 115
hist(log2(l115.12[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 12 Line 115")
lines(density(log2(l115.12[,3])),col="green")
mc11 <- mean(log2(l115.12[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l115.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .3, labels = mylabel,cex=.8)
sd2 <- sd(log2(l115.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .5, labels = mylabel3,cex=.8)

#line 115
hist(log2(l117.12[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 12 Line 117")
lines(density(log2(l117.12[,3])),col="green")
mc11 <- mean(log2(l117.12[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l117.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .3, labels = mylabel,cex=.8)
sd2 <- sd(log2(l117.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .5, labels = mylabel3,cex=.8)

#line 123
hist(log2(l123.12[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Trisomic Chr 12 Line 123")
lines(density(log2(l123.12[,3])),col="green")
mc11 <- mean(log2(l123.12[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l123.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .3, labels = mylabel,cex=.8)
sd2 <- sd(log2(l123.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .5, labels = mylabel3,cex=.8)

#line 141
hist(log2(l141.12[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 12 Line 141")
lines(density(log2(l141.12[,3])),col="green")
mc11 <- mean(log2(l141.12[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l141.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel,cex=.8)
sd2 <- sd(log2(l141.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .5, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .6, labels = mylabel3,cex=.8)

#line 29
hist(log2(l29.12[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 12 Line 29")
lines(density(log2(l29.12[,3])),col="green")
mc11 <- mean(log2(l29.12[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l29.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel,cex=.8)
sd2 <- sd(log2(l29.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .5, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .6, labels = mylabel3,cex=.8)

#line 50
hist(log2(l50.12[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 12 Line 50")
lines(density(log2(l50.12[,3])),col="green")
mc11 <- mean(log2(l50.12[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l50.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .3, labels = mylabel,cex=.8)
sd2 <- sd(log2(l50.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .5, labels = mylabel3,cex=.8)

#line 152
hist(log2(l152.12[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 12 Line 152")
lines(density(log2(l152.12[,3])),col="green")
mc11 <- mean(log2(l152.12[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l152.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l152.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)


#CHROMOSOME 16
#line 112
hist(log2(l112.16[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Trisomic Chr 16 Line 112")
lines(density(log2(l112.16[,3])),col="green")
mc11 <- mean(log2(l112.16[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l112.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l112.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#line 115
hist(log2(l115.16[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 16 Line 115")
lines(density(log2(l115.16[,3])),col="green")
mc11 <- mean(log2(l115.16[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l115.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .3, labels = mylabel,cex=.8)
sd2 <- sd(log2(l115.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .5, labels = mylabel3,cex=.8)

#line 117
hist(log2(l117.16[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 16 Line 117")
lines(density(log2(l117.16[,3])),col="green")
mc11 <- mean(log2(l117.16[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l117.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .25, labels = mylabel,cex=.8)
sd2 <- sd(log2(l117.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .35, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .45, labels = mylabel3,cex=.8)

#line 123
hist(log2(l123.16[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 16 Line 123")
lines(density(log2(l123.16[,3])),col="green")
mc11 <- mean(log2(l123.16[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l123.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .3, labels = mylabel,cex=.8)
sd2 <- sd(log2(l123.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .5, labels = mylabel3,cex=.8)

#line 141
hist(log2(l141.16[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Trisomic Chr 16 Line 141")
lines(density(log2(l141.16[,3])),col="green")
mc11 <- mean(log2(l141.16[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l141.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l141.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)

#line 29
hist(log2(l29.16[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 16 Line 29")
lines(density(log2(l29.16[,3])),col="green")
mc11 <- mean(log2(l29.16[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l29.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .4, labels = mylabel,cex=.8)
sd2 <- sd(log2(l29.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .5, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .6, labels = mylabel3,cex=.8)

#line 50
hist(log2(l50.16[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 16 Line 50")
lines(density(log2(l50.16[,3])),col="green")
mc11 <- mean(log2(l50.16[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l50.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .25, labels = mylabel,cex=.8)
sd2 <- sd(log2(l50.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .35, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .45, labels = mylabel3,cex=.8)

#line 152
hist(log2(l152.16[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 16 Line 152")
lines(density(log2(l152.16[,3])),col="green")
mc11 <- mean(log2(l152.16[,3]))
mc11
abline(v = mc11, col = "blue")
s2 <- skewness(log2(l152.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l152.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .85, labels = mylabel3,cex=.8)


##############################################################################################
# let's look at the bottom X% of genes in the samples
#############
go1 <- NULL;go2 <- NULL;go3 <- NULL;go4 <- NULL;go5 <- NULL;go6 <- NULL;go7 <- NULL;go8 <- NULL;
go9 <- NULL;go10 <- NULL;go11 <- NULL;go12 <- NULL;go13 <- NULL;go14 <- NULL;go15 <- NULL;go16 <- NULL


########################
for(i in 1:16){ #for each go1 go2 go3 etc
  c <- NULL
  for (line in rep.MA){ #for each line
    c <- cbind(c,as.character(getChrmRatio.MA(line,i)[,3]))
    #assign(paste("go",i,sep=""),value=as.character(getChrmRatio(line,i)[,2]))
  }
  assign(paste("go",i,sep=""),value=c)
} #this gets the gene order on each line for each chromsome

#########################
colnames(go1) <- rep.MA
colnames(go2) <- rep.MA
colnames(go3) <- rep.MA
colnames(go4) <- rep.MA
colnames(go5) <- rep.MA
colnames(go6) <- rep.MA
colnames(go7) <- rep.MA
colnames(go8) <- rep.MA
colnames(go9) <- rep.MA
colnames(go10) <- rep.MA
colnames(go11) <- rep.MA
colnames(go12) <- rep.MA
colnames(go13) <- rep.MA
colnames(go14) <- rep.MA
colnames(go15) <- rep.MA
colnames(go16) <- rep.MA
##########################
#makes colnames X001, X002, etc
go1 <- data.frame(go1)
length(intersect(intersect(intersect(intersect(intersect(intersect(intersect(go1$X112[1:50],go1$X115[1:50]),go1$X117[1:50]),go1$X123[1:50]),go1$X141[1:50]),go1$X152[1:50]),go1$X29[1:50]),go1$X50[1:50]))

#no shared genes between the anueploid lines
#look at euploid lines next 

length(intersect(intersect(intersect(intersect(go1$X2[1:50],go1$X3[1:50]),go1$X5[1:50]),go1$X6[1:50]),go1$X9[1:50]))
#only one gene shared between 5 euploid lines 

#if we compare sets of 3 of the euploids, what do we get?

#make the list of reps we should choose from
xrep <- c(1:8,10,12)

#now make a loop to randomly sample these
overlap <- NULL
for(i in 1:1000){
  triple <- sample(xrep,3,replace=F)
  bounds <- 1:75
  overlap[i] <- length(intersect(intersect(go1[bounds,triple[1]],go1[bounds,triple[2]]),go1[bounds,triple[3]]))
}
mean(overlap)
max(overlap)
min(overlap)
