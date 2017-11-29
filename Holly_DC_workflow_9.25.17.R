#Holly's DC scripts, modified from DCfinal.R
#DC with the library-size normalized counts

###Rate of aneuploidy for GC and MA data 
#are they statistically significantly different?
#GC_rate <- 1.9969e-4 
#MA_rate <- 9.7e-5
#rate <- c(GC_rate,MA_rate)
#chisq.test(rate)

library(moments)
##########read the data
sc_counts <- read.csv("/Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Redo/September_2017/genes_fpkm_sorted_paste_colrm.csv",header=T)
###colnames(sc_counts) <- as.character(unlist(sc_counts[1,]))
##need to find rows in which second column is MT 
###sc_counts = sc_counts[-1,]
class(sc_counts)
str(sc_counts)
#colnames(sc_counts)[1] <- "gene"; colnames(sc_counts)[2] <- "class_code";colnames(sc_counts)[3] <- "nearest_ref_id"
#colnames(sc_counts)[4] <- "ref_id"; colnames(sc_counts)[5] <- "gene_short_name"; colnames(sc_counts)[6] <- "tss_id";
#colnames(sc_counts)[7] <- "locus";colnames(sc_counts)[8] <- "length";colnames(sc_counts)[9] <- "null";colnames(sc_counts)[10] <- "X112_0";
#colnames(sc_counts)[11] <- "X112_1";colnames(sc_counts)[12] <- "X112_2";colnames(sc_counts)[13] <- "X115_0";
#colnames(sc_counts)[14] <- "X115_1";colnames(sc_counts)[15] <- "X115_2";colnames(sc_counts)[16] <- "X117_0";
#colnames(sc_counts)[17] <- "X117_1";colnames(sc_counts)[18] <- "X117_2";colnames(sc_counts)[19] <- "X11_1";
#colnames(sc_counts)[20] <- "X11_0";colnames(sc_counts)[21] <- "X11_2";colnames(sc_counts)[22] <- "X123_0";
#colnames(sc_counts)[23] <- "X123_1";colnames(sc_counts)[24] <- "X123_2";colnames(sc_counts)[25] <- "X141_0";
#colnames(sc_counts)[26] <- "X141_1";colnames(sc_counts)[27] <- "X141_2";colnames(sc_counts)[28] <- "X152_1";
#colnames(sc_counts)[29] <- "X152_0";colnames(sc_counts)[30] <- "X152_2";colnames(sc_counts)[31] <- "X18_1";
#colnames(sc_counts)[32] <- "X18_0";colnames(sc_counts)[33] <- "X18_2";colnames(sc_counts)[34] <- "X1_0";
#colnames(sc_counts)[35] <- "X1_1";colnames(sc_counts)[36] <- "X1_2";colnames(sc_counts)[37] <- "X21_1";
#colnames(sc_counts)[38] <- "X21_0";colnames(sc_counts)[39] <- "X21_2";colnames(sc_counts)[40] <- "X29_1";
#colnames(sc_counts)[41] <- "X29_0";colnames(sc_counts)[42] <- "X29_2";colnames(sc_counts)[43] <- "X2_0";
#colnames(sc_counts)[44] <- "X2_1";colnames(sc_counts)[45] <- "X2_2";colnames(sc_counts)[46] <- "X31_1";
#colnames(sc_counts)[47] <- "X31_0";colnames(sc_counts)[48] <- "X31_2";colnames(sc_counts)[49] <- "X3_0";
#colnames(sc_counts)[50] <- "X3_1";colnames(sc_counts)[51] <- "X3_2";colnames(sc_counts)[52] <- "X49_0";
#colnames(sc_counts)[53] <- "X49_1";colnames(sc_counts)[54] <- "X49_2";colnames(sc_counts)[55] <- "X4_1";
#colnames(sc_counts)[56] <- "X4_0";colnames(sc_counts)[57] <- "X4_2";colnames(sc_counts)[58] <- "X50_1";
#colnames(sc_counts)[59] <- "X50_0";colnames(sc_counts)[60] <- "X50_2";colnames(sc_counts)[61] <- "X59_0";
#colnames(sc_counts)[62] <- "X59_1";colnames(sc_counts)[63] <- "X59_2";colnames(sc_counts)[64] <- "X5_1";
#colnames(sc_counts)[65] <- "X5_0";colnames(sc_counts)[66] <- "X5_2";colnames(sc_counts)[67] <- "X61_1";
#colnames(sc_counts)[68] <- "X61_0";colnames(sc_counts)[69] <- "X61_2";colnames(sc_counts)[70] <- "X66_0";
#colnames(sc_counts)[71] <- "X66_1";colnames(sc_counts)[72] <- "X66_2";colnames(sc_counts)[73] <- "X69_1";
#colnames(sc_counts)[74] <- "X69_0";colnames(sc_counts)[75] <- "X69_2";colnames(sc_counts)[76] <- "X6_0";
#colnames(sc_counts)[77] <- "X6_1";colnames(sc_counts)[78] <- "X6_2";colnames(sc_counts)[79] <- "X76_0";
#colnames(sc_counts)[80] <- "X76_1";colnames(sc_counts)[81] <- "X76_2";colnames(sc_counts)[82] <- "X77_1";
#colnames(sc_counts)[83] <- "X77_0";colnames(sc_counts)[84] <- "X77_2";colnames(sc_counts)[85] <- "X7_1";
#colnames(sc_counts)[86] <- "X7_0";colnames(sc_counts)[87] <- "X7_2";colnames(sc_counts)[88] <- "X8_0";
#colnames(sc_counts)[89] <- "X8_1";colnames(sc_counts)[90] <- "X8_2";colnames(sc_counts)[91] <- "X9_0";
#colnames(sc_counts)[92] <- "X9_1";colnames(sc_counts)[93] <- "X9_2";colnames(sc_counts)[94] <- "GC_Anc_1";
#colnames(sc_counts)[95] <- "GC_Anc_0";colnames(sc_counts)[96] <- "GC_Anc_2";colnames(sc_counts)[97] <- "MA_Anc_1";
#colnames(sc_counts)[98] <- "MA_Anc_0";colnames(sc_counts)[99] <- "MA_Anc_2";

##remove the top row which is the old header
#sc_counts <- sc_counts[-c(1), ]
#remove problem rows-ones that have nothing in them
#remove the columns we don't need
#sc_counts.a <- sc_counts[,-c(1,2,3,4,6,8,9) ]
#class(sc_counts.a)
#str(sc_counts.a)
colnames(sc_counts)[1] <- "chr"; colnames(sc_counts)[2] <- "tracking_id"
#decided to keep length
#sc_counts.b <- sc_counts.a[,-c(3) ]



#sc_counts.a <- data.matrix(sc_counts.a)
#sc_counts.a <- data.frame(sc_counts.a)
#sc_counts.a$X112_0 <-as.numeric(as.character(sc_counts.a[,3]));sc_counts.a$X112_1 <-as.numeric(as.character(sc_counts.a[,4]));
#sc_counts.a$X112_2 <-as.numeric(as.character(sc_counts.a[,2]));
#test to see if you can add .5 now to columns you changed to numeric 
#sc_counts.z[,3:4] <- sc_counts.z[,3:4]+.5 

sc_counts.z <- sc_counts
sc_counts.z[,3:140] <- sc_counts.z[,3:140]+.5 #add .5 to dataset so the zeroes dont mess up


warnings()
#we want to exclude genes that are very lowly expressed
gene.means <- rowMeans(sc_counts.z[,3:140]) #average across the sample
anc.means.MA <- rowMeans(sc_counts.z[,90:92]) #average in the ancestor
anc.means.GC <- rowMeans(sc_counts.z[87:89]) #average in the GC ancestor
anc.means.old <- rowMeans(sc_counts.z[,138:140]) #average in the old ancestor

#if the average fpkm is less than 5 in the ancestor or less than 5 in the overall, remove the gene
sc_counts.final <- sc_counts.z[-which(gene.means<5 | anc.means.MA<5 | anc.means.GC<5 |anc.means.old<5),]

#final matrix for ancestor and descendans
sc_anc <- sc_counts.final[,c(1,2,87,88,89,90,91,92,138,139,140)]
sc_des.GC <- sc_counts.final[,c(1:2,12:14,24:32,36:50,54:86)] #average for each GC line
sc_des.MA <- sc_counts.final[,c(1:11,15:23,33:35,51:53)] #average for each MA line
sc_des.old <- sc_counts.final[,c(1:2,93:137)] #average for each old (MA) line 
###########
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
#PROBLEM HERE! I was using the old code and not taking into account that I had two ancestors! now I have 3...
ancAvg.GC <- rowMeans(sc_anc[,3:5],na.rm=T)
ancAvg.GC <- cbind(sc_anc[,1:2],ancAvg.GC)
ancAvg.MA <- rowMeans(sc_anc[,6:8],na.rm=T)
ancAvg.MA <- cbind(sc_anc[,1:2],ancAvg.MA)
ancAvg.old <- rowMeans(sc_anc[,9:11],na.rm=T)
ancAvg.old <- cbind(sc_anc[,1:2],ancAvg.old)
#ancMin <- min(sc_anc[,3:6],na.rm=T) 
#ancMin <- cbind(sc_anc[,1:2],ancMin)
#ancMax <- max(sc_anc[,3:6],na.rm=T)
#ancMax <- cbind(sc_anc[,1:2],ancMax)

#first write some functions that help automate the process
##########
#this gets the average FPKM for a given line from the sc_des matrix
getAvgFPKM.MA <- function(line) { #for line as "###" string
  r1 <- paste(c("X",line,"_0"),collapse="") #triplicates are named X[LINE]_[rep]
  r2 <- paste(c("X",line,"_1"),collapse="")
  r3 <- paste(c("X",line,"_2"),collapse="")
  reps <- sc_des.MA[,c("chr","tracking_id",r1,r2,r3)] #get those 4 columns out of the descendant dataset
  reps[reps==0] <- NA #replace 0 values with NA as before
  desAvg.MA <- rowMeans(reps[,c(-1,-2)],na.rm=T)
  desAvg.MA <- cbind(sc_des.MA[,c(1,2)],desAvg.MA)
  #colnames(desAvg) <- c("chr","gene",paste(c(line,"_Avg"),collapse=""))
  return(desAvg.MA) #returned is a matrix equivalent to ancAvg: chrm, gene, and line average FPKM 
}


getAvgFPKM.GC <- function(line) { #for line as "###" string
  r1 <- paste(c("X",line,"_0"),collapse="") #triplicates are named X[LINE]_[rep]
  r2 <- paste(c("X",line,"_1"),collapse="")
  r3 <- paste(c("X",line,"_2"),collapse="")
  reps <- sc_des.GC[,c("chr","tracking_id",r1,r2,r3)] #get those 4 columns out of the descendant dataset
  reps[reps==0] <- NA #replace 0 values with NA as before
  desAvg.GC <- rowMeans(reps[,c(-1,-2)],na.rm=T)
  desAvg.GC <- cbind(sc_des.GC[,c(1,2)],desAvg.GC)
  #colnames(desAvg) <- c("chr","gene",paste(c(line,"_Avg"),collapse=""))
  return(desAvg.GC) #returned is a matrix equivalent to ancAvg: chrm, gene, and line average FPKM 
}

getAvgFPKM.old <- function(line) { #for line as "###" string
  r1 <- paste(c("X",line,"_0"),collapse="") #triplicates are named X[LINE]_[rep]
  r2 <- paste(c("X",line,"_1"),collapse="")
  r3 <- paste(c("X",line,"_2"),collapse="")
  reps <- sc_des.old[,c("chr","tracking_id",r1,r2,r3)] #get those 4 columns out of the descendant dataset
  reps[reps==0] <- NA #replace 0 values with NA as before
  desAvg.old <- rowMeans(reps[,c(-1,-2)],na.rm=T)
  desAvg.old <- cbind(sc_des.old[,c(1,2)],desAvg.old)
  #colnames(desAvg) <- c("chr","gene",paste(c(line,"_Avg"),collapse=""))
  return(desAvg.old) #returned is a matrix equivalent to ancAvg: chrm, gene, and line average FPKM 
}

##to get minimum FPKM for a given line from the sc_des matrix
#getminFPKM <- function(line) {
# r1 <- paste(c("X",line,"_0"),collapse="") #triplicates are named X[LINE]_[rep]
#r2 <- paste(c("X",line,"_1"),collapse="")
#  r3 <- paste(c("X",line,"_2"),collapse="")
#  reps <- sc_des[,c("tracking_id","chr",r1,r2,r3)] #get those 4 columns out of the descendant dataset
#  reps[reps==0] <- NA #replace 0 values with NA as before
#  desMin <- pmin(reps[,c(-1,-2)],na.rm=T)
#  desMin <- cbind(sc_des[,c(1,2)],desMin)
#  return(desMin) #returned is a matrix equivalent to ancAvg: chrm, gene, and line average FPKM 
#}
##to get maximum FPKM for a given line from the sc_des matrix
#getmaxFPKM <- function(line) {
#  r1 <- paste(c("X",line,"_0"),collapse="") #triplicates are named X[LINE]_[rep]
#  r2 <- paste(c("X",line,"_1"),collapse="")
#  r3 <- paste(c("X",line,"_2"),collapse="")
#  reps <- sc_des[,c("tracking_id","chr",r1,r2,r3)] #get those 4 columns out of the descendant dataset
#  reps[reps==0] <- NA #replace 0 values with NA as before
#  desMax <- pmax(reps[,c(-1,-2)],na.rm=T)
#  desMax <- cbind(sc_des[,c(1,2)],desMax)
#  return(desMax) #returned is a matrix equivalent to ancAvg: chrm, gene, and line average FPKM 
#}

##########
#this gets the fold change for the whole line
getFPKMRatio.MA <- function(line) { #line must be a 3 digit character
  desAvg.MA <- getAvgFPKM.MA(line) #first call avg fpkm script, which returns a 3 column matrix
  ratioByGene.MA <- desAvg.MA[,3]/ancAvg.MA[,3]
  rat.MA <- cbind(sc_des.MA[,1:2],ratioByGene.MA) #add gene and chrm values
  colnames(rat.MA) <- c("chr","tracking_id","ratio")
  return(rat.MA)
}

getFPKMRatio.GC <- function(line) { #line must be a 3 digit character
  desAvg.GC <- getAvgFPKM.GC(line) #first call avg fpkm script, which returns a 3 column matrix
  ratioByGene.GC <- desAvg.GC[,3]/ancAvg.GC[,3]
  rat.GC <- cbind(sc_des.GC[,1:2],ratioByGene.GC) #add gene and chrm values
  colnames(rat.GC) <- c("chr","tracking_id","ratio")
  return(rat.GC)
}

getFPKMRatio.old <- function(line) { #line must be a 3 digit character
  desAvg.old <- getAvgFPKM.old(line) #first call avg fpkm script, which returns a 3 column matrix
  ratioByGene.old <- desAvg.old[,3]/ancAvg.old[,3]
  rat.old <- cbind(sc_des.old[,1:2],ratioByGene.old) #add gene and chrm values
  colnames(rat.old) <- c("chr","tracking_id","ratio")
  return(rat.old)
}

##########
#this gets the fold change for the minimum
#getMinRatio <- function(line) { #line must be a 3 digit character
#  desMin <- getminFPKM(line) #first call avg fpkm script, which returns a 3 column matrix
#  ratioByGene <- desMin[,3]/ancMin[,3]
#  rat <- cbind(sc_des[,1:2],ratioByGene) #add gene and chrm values
#  colnames(rat) <- c("chrm","gene","ratio")
#  return(rat)
#}

##this gets the fold change for the maximum 
#getMaxRatio <- function(line) { #line must be a 3 digit character
#  desMax <- getmaxFPKM(line) #first call avg fpkm script, which returns a 3 column matrix
#  ratioByGene <- desMax[,3]/ancMax[,3]
#  rat <- cbind(sc_des[,1:2],ratioByGene) #add gene and chrm values
#  colnames(rat) <- c("chrm","gene","ratio")
#  return(rat)
#}

##########
#this gets the fold change for just the specified chromosome
##because my columns of gene and chr were switched, I had to specify in this function that chr was now column2 
#had to switch them back #thanksforcommentingoldme
getChrmRatio.MA <- function(line, chr){ #input line as "XXX", chr as a number
  rat.MA <- getFPKMRatio.MA(line) #get list of ratios which includes chrm and gene
  rat.c.MA <- rat.MA[(rat.MA[,1]==chr),] #get just the rows for the appropriate chrom
  rat.cs.MA <- rat.c.MA[order(rat.c.MA[,3]),]#order by column 3
  return(rat.cs.MA)
}

getChrmRatio.GC <- function(line, chr){ #input line as "XXX", chr as a number
  rat.GC <- getFPKMRatio.GC(line) #get list of ratios which includes chrm and gene
  rat.c.GC <- rat.GC[(rat.GC[,1]==chr),] #get just the rows for the appropriate chrom
  rat.cs.GC <- rat.c.GC[order(rat.c.GC[,3]),]#order by column 3
  return(rat.cs.GC)
}

getChrmRatio.old <- function(line, chr){ #input line as "XXX", chr as a number
  rat.old <- getFPKMRatio.old(line) #get list of ratios which includes chrm and gene
  rat.c.old <- rat.old[(rat.old[,1]==chr),] #get just the rows for the appropriate chrom
  rat.cs.old <- rat.c.old[order(rat.c.old[,3]),]#order by column 3
  return(rat.cs.old)
}

rat.old <- getFPKMRatio.old(line) #get list of ratios which includes chrm and gene
rat.c.old <- rat.old[(rat.old[,1]==chr),] #get just the rows for the appropriate chrom
rat.cs.old <- rat.c.old[order(rat.c.old[,3]),]#order by column 3
##########

##########
#this gets the minimum fold change for just the specified chromosome
#getMinChrmRatio <- function(line, chr){ #input line as "XXX", chr as a number
#  rat <- getMinRatio(line) #get list of ratios which includes chrm and gene
#  rat.c <- rat[(rat[,1]==chr),] #get just the rows for the appropriate chrom
#  rat.cs <- rat.c[order(rat.c[,3]),]#order by column 3
#  return(rat.cs)
#}

#getMaxChrmRatio <- function(line, chr){ #input line as "XXX", chr as a number
#  rat <- getMaxRatio(line) #get list of ratios which includes chrm and gene
#  rat.c <- rat[(rat[,1]==chr),] #get just the rows for the appropriate chrom
#  rat.cs <- rat.c[order(rat.c[,3]),]#order by column 3
#  return(rat.cs)
#}
#now we can look at the fold changes for each chromosome and each line

#make a matrix of median ratio:

#make a matrix of min ratio:
#ymin=NULL
#rep <- c("112","115","117","11","123","141","152","18","1","21","29","2","31","3","49","4","50","59","5","61","66","69","6","76","77","7","8","9")
#for(line in rep){
#  for(chrm in 1:16){
#    ymin <- c(ymin,getMinChrmRatio(line,chrm)[,3],na.rm=T)
#  }
#}

#make a matrix of max ratio:
#ymax=NULL
#rep <- c("112","115","117","11","123","141","152","18","1","21","29","2","31","3","49","4","50","59","5","61","66","69","6","76","77","7","8","9")
#for(line in rep){
#  for(chrm in 1:16){
#    ymax <- c(ymax,getMaxChrmRatio(line,chrm)[,3],na.rm=T)
#  }
#}

#make a matrix of median ratio:
#y=NULL
#rep <- c("X112","X115","X117","X11","X123","X141","X152","X18","X1","X21","X29","X2","X31","X3","X49","X4","X50","X59","X5","X61","X66","X69","X6","X76","X77","X7","X8","X9")
#for(line in rep){
# for(chrm in 1:16){
#  y <- c(y,median(getChrmRatio(line,chrm)[,3],na.rm=T))
#}
#}

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

#remove the ratios that are about 1
ratMatMedRed.MA <- ratMatMed.MA
ratMatMedRed.MA[ratMatMedRed.MA<1.35] <- NA
#ratMatMedRed.MA[9,15] <- ratMatMedRed.MA[9,15]

y.GC=NULL
rep.GC <- c("11","18","1","21","2","31","3","49","4","59","5","61","66","69","6","76","77","7","8","9")
for(line in rep.GC){
  for(chrm in 1:16){
    y.GC <- c(y.GC,median(getChrmRatio.GC(line,chrm)[,3],na.rm=T))
  }
}

ratMatMed.GC <- matrix(y.GC,nrow=16,ncol=20)
colnames(ratMatMed.GC) <- rep.GC
ratMatMed.GC <- cbind(c(1:16),ratMatMed.GC)
colnames(ratMatMed.GC)[1] <- "chr"

#remove the ratios that are about 1
#ratMatMedRed.GC <- ratMatMed.GC
#ratMatMedRed.GC[ratMatMedRed.GC<1.35] <- NA
#ratMatMedRed.MA[9,15] <- ratMatMedRed.MA[9,15]

y.old=NULL
rep.old <- c("001","002","003","004","005","006","007","008","009","011","015","028","088","108","119")
for(line in rep.old){
  for(chrm in 1:16){
    y.old <- c(y.old,median(getChrmRatio.old(line,chrm)[,3],na.rm=T))
  }
}

#make a matrix out of it
ratMatMed.old <- matrix(y.old,nrow=16,ncol=15)
colnames(ratMatMed.old) <- rep.old
ratMatMed.old <- cbind(c(1:16),ratMatMed.old)
colnames(ratMatMed.old)[1] <- "chr"

#remove the ratios that are about 1
#ratMatMedRed.old <- ratMatMed.old
#ratMatMedRed.old[ratMatMed.old<1.35] <- NA


#make a matrix out of it- min 
#ratMatMedMin <- matrix(ymin,nrow=16,ncol=28)
#colnames(ratMatMedMin) <- rep
#ratMatMedMin <- cbind(c(1:16),ratMatMedMin)


#make a matrix out of it- min 
#ratMatMedMax <- matrix(ymax,nrow=16,ncol=28)
#colnames(ratMatMedMax) <- rep
#ratMatMedMax <- cbind(c(1:16),ratMatMedMax)



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

yl.GC=NULL
for(line in rep.GC){
  for(chrm in 1:16){
    yl.GC <- c(yl.GC,median(log2(getChrmRatio.GC(line,chrm)[,3])))
  }
}

warnings()
logMat.GC <- matrix(yl.GC,nrow=16,ncol=20)
colnames(logMat.GC) <- rep.GC
logMat.GC <- data.frame(logMat.GC)
chrms <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
logMat2.GC <-cbind(chrms,logMat.GC)

#this makes a matrix with the log2 ratio between ancestor and sample of each gene 

getChrmRatio.GC.test <- function(line, chr){ #input line as "XXX", chr as a number
  rat.GC <- getFPKMRatio.GC(line) #get list of ratios which includes chrm and gene
  rat.cs.GC <- rat.GC[order(rat.GC[,3]),]#order by column 3
  return(rat.cs.GC)
}


yl.GC=NULL
for(line in rep.GC){
    y.GC <- c(yl.GC,log2(getChrmRatio.GC.test(line,chrm)[,3]))
}

warnings()
logMat.GC.test<- matrix(y.GC,nrow=5431,ncol=20)
colnames(logMat.GC.test) <- rep.GC
logMat.GC.test <- data.frame(logMat.GC.test)
genes <- ancAvg.GC[,2]
logMat2.GC.test <-cbind(genes,logMat.GC.test)

yl.old=NULL
for(line in rep.old){
  for(chrm in 1:16){
    yl.old <- c(yl.old,median(log2(getChrmRatio.old(line,chrm)[,3])))
  }
}
View(yl.old)

warnings()
logMat.old <- matrix(yl.old,nrow=16,ncol=15)
colnames(logMat.old) <- rep.old
logMat.old <- data.frame(logMat.old)
chrms <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
logMat2.old <-cbind(chrms,logMat.old)

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

#for GC lines
ysd.GC=NULL
for(line in rep.GC){
  for(chrm in 1:16){
    ysd.GC <- c(ysd.GC,sd(log2(getChrmRatio.GC(line,chrm)[,3])))
  }
}


logMatsd.GC <- matrix(ysd.GC,nrow=16,ncol=20)
colnames(logMatsd.GC) <- rep.GC
logMatsd.GC <- data.frame(logMatsd.GC)
chrms <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
logMatsd2.GC <-cbind(chrms,logMatsd.GC)

#for old lines
ysd.old=NULL
for(line in rep.old){
  for(chrm in 1:16){
    ysd.old <- c(ysd.old,sd(log2(getChrmRatio.old(line,chrm)[,3])))
  }
}

logMatsd.old <- matrix(ysd.old,nrow=16,ncol=15)
colnames(logMatsd.old) <- rep.old
logMatsd.old <- data.frame(logMatsd.old)
chrms <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
logMatsd2.old <-cbind(chrms,logMatsd.old)

warnings()
##slightly confusing, but had to change the minimum and maximim logMatrices because they would've been opposite and not useful for error bars

#make a matrix that is the LOG2 RATIO for the minimums 
#ymin=NULL
#for(line in rep){
#  for(chrm in 1:16){
#    ymin <- c(ymin,(log2(getMinRatio(line)[,3])))
#  }
#}

#logMinMat <- matrix(ymin,nrow=16,ncol=28)
#colnames(logMinMat) <- c("112","115","117","11","123","141","152","18","1","21","29","2","31","3","49","4","50","59","5","61","66","69","6","76","77","7","8","9")
#logMinMat <- data.frame(logMinMat)
#chrms <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
#logMaxMat2 <-cbind(chrms,logMinMat)


#make a matrix that is the LOG2 RATIO for the maximums 
#ymax=NULL
#for(line in rep){
#  for(chrm in 1:16){
#    ymax <- c(ymax,(log2(getMaxRatio(line)[,3])))
#  }
#}

#logMaxMat <- matrix(ymax,nrow=16,ncol=28)
#colnames(logMaxMat) <- c("112","115","117","11","123","141","152","18","1","21","29","2","31","3","49","4","50","59","5","61","66","69","6","76","77","7","8","9")
#logMaxMat <- data.frame(logMaxMat)
#chrms <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
#logMinMat2 <-cbind(chrms,logMaxMat)


#make a matrix that is the maximum fpkm 
#fmax=NULL
#for(line in rep){
 # for(chrm in 1:16){
  #  fmax <- c(fmax,(log2(getmaxFPKM(line)[,3])))
  #}
#}

#MaxMat <- matrix(fmax,nrow=16,ncol=28)
#colnames(logMaxMat) <- c("112","115","117","11","123","141","152","18","1","21","29","2","31","3","49","4","50","59","5","61","66","69","6","76","77","7","8","9")
#logMaxMat <- data.frame(logMaxMat)
#chrms <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
#logMaxMat2 <-cbind(chrms,logMaxMat)

##get the column of log2 values of line 112 from the logMat (matrix) of all of the log2 values
#X112 <- subset(logMat, select=c("X112"))
#make a table with the chromosome numer in there too 
#X112c <- cbind(chrms,X112)
##get the maximum log2 fold changes out 
#X112max <- subset(logMaxMat2, select=c("X112"))
#X112maxc <- cbind(chrms,X112max)
#get the minimum log2 fold changes out 
#X112min <- subset(logMinMat2, select=c("X112"))
#X112minc <- cbind(chrms,X112min)
#####these work, but don't do exactly what I want them to do 
##library(lattice)
##xyplot(X112~chrms, data=X112,ylim = (-1.5:1.5))
###y.bar <- mean()
###y.SE <- c(sd(yL)/sqrt(length(yL)), sd(yH)/sqrt(length(yH)))
##plot the log2 values by the chromosome # 
#plot(X112~chrms, data=X112,ylim=c(-1.5,1.5),xlab="Chromosome", ylab="log2(fold change)",type="p",main="MA Line 112",pch=16,col="darkturquoise")
#this is just an example plot script
####plot(1:2, y.bar, ylim=c(0, 25), xlim=c(0.5, 2.5), xaxt="n", pch=16,
     ###xlab="", ylab="Tree density", cex=1.5)
##this makes the x axis have the chromosome # labels
#axis(1, 1:16, c("1", "2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"))
##next I would like to figure out how to put error bars, but it is not completely necessary right now 
##add error bars
#d = data.frame(
#  x  = c(1:16)
#  , y  = c(1.1, 1.5, 2.9, 3.8, 5.2)
#  , sd = c(0.2, 0.3, 0.2, 0.0, 0.4)
#)

##install.packages("Hmisc", dependencies=T)

X112 <- subset(logMat.MA, select=c("X112"))
#X112.1 <-X112[1,]
#X112c <- cbind(chrms,X112)
#X112.sd <- subset(logMatsd2, select=c("X112"))
#X112.1.sd <- X112.sd[1,]
#arrows(1, X112.1-X112.1.sd, 1, X112.1+X112.1.sd, code=3,angle=90,length=0.05)

#View(X112.1)
#X112c <- cbind(chrms,X112)

#X112.sd <- subset(logMatsd2, select=c("X112"))
#X112.all <- cbind(chrms,X112,X112.sd)
#colnames(X112.all)[1] <- "chr"; colnames(X112.all)[2] <- "mean";colnames(X112.all)[3] <- "sd"
#X112.sd1 <- X112.sd[1,]
#arrows(1:2, X112-X112.sd1, 1:2, X112+X112.sd1, code=3,angle=90,length=0.05)
#arrows(x, avg-sdev, x, avg+sdev, length=0.05, angle=90, code=3)

###arrows(1:2, X112-X112min, 1:2, X112+X112max, code=3, angle=90, length=0.05)
#x <- 1:16
#plot(X112~chrms, data=X112,
 #    ylim=range(c(X112-X112.sd, X112+X112.sd)),
  #   pch=16, xlab="Chromosome", ylab="log2(fold change",
   #  main="MA Line 112", col="darkturquoise"
#)

###arrows(x, X112-X112.sd, x, X112+X112.sd, length=0.05, angle=90, code=3)

######################################################################################
#how to make a graph of the chromosomes and their respective error bars 
library("Hmisc")
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

#lets try making a loop just for fun 
for (line in c(112,115,117,123,141,152,29,50)) {
  
line <- subset(logMat.MA, select=line)
line.sd <- subset(logMatsd2.MA, select=line)
#merge the stdev with the mean and the chromosome number
line.all <- cbind(chrms,line,line.sd)
#change the column names appropriately
colnames(line.all)[1] <- "chr"; colnames(line.all)[2] <- "mean";colnames(line.all)[3] <- "sd"

#plots the chromosomes with error bars for each one 
with (
  data = line.all
  , expr = errbar(main="MA",chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)


#adds a title to the graph
abline(h=0,lty=3)
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

}
#well, that doesn't work...


#line 115
par(mar=c(5, 4, 4, 2) + 0.1)
#get the standard deviations for each chromosome for line 112 only
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
#get the standard deviations for each chromosome for line 112 only
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

#line 11 
X11 <- subset(logMat.GC, select=c("X11"))
X11.sd <- subset(logMatsd2.GC, select=c("X11"))
#merge the stdev with the mean and the chromosome number
X11.all <- cbind(chrms,X11,X11.sd)
#change the column names appropriately
colnames(X11.all)[1] <- "chr"; colnames(X11.all)[2] <- "mean";colnames(X11.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X11.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 11")
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

#line 18
X18 <- subset(logMat.GC, select=c("X18"))
X18.sd <- subset(logMatsd2.GC, select=c("X18"))
#merge the stdev with the mean and the chromosome number
X18.all <- cbind(chrms,X18,X18.sd)
#change the column names appropriately
colnames(X18.all)[1] <- "chr"; colnames(X18.all)[2] <- "mean";colnames(X18.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X18.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 18")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 1
X1 <- subset(logMat.GC, select=c("X1"))
X1.sd <- subset(logMatsd2.GC, select=c("X1"))
#merge the stdev with the mean and the chromosome number
X1.all <- cbind(chrms,X1,X1.sd)
#change the column names appropriately
colnames(X1.all)[1] <- "chr"; colnames(X1.all)[2] <- "mean";colnames(X1.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X1.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 1")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 21
X21 <- subset(logMat.GC, select=c("X21"))
X21.sd <- subset(logMatsd2.GC, select=c("X21"))
#merge the stdev with the mean and the chromosome number
X21.all <- cbind(chrms,X21,X21.sd)
#change the column names appropriately
colnames(X21.all)[1] <- "chr"; colnames(X21.all)[2] <- "mean";colnames(X21.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X21.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 21")
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

#line 2
X2 <- subset(logMat.GC, select=c("X2"))
X2.sd <- subset(logMatsd2.GC, select=c("X2"))
#merge the stdev with the mean and the chromosome number
X2.all <- cbind(chrms,X2,X2.sd)
#change the column names appropriately
colnames(X2.all)[1] <- "chr"; colnames(X2.all)[2] <- "mean";colnames(X2.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X2.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 2")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 31
X31 <- subset(logMat.GC, select=c("X31"))
X31.sd <- subset(logMatsd2.GC, select=c("X31"))
#merge the stdev with the mean and the chromosome number
X31.all <- cbind(chrms,X31,X31.sd)
#change the column names appropriately
colnames(X31.all)[1] <- "chr"; colnames(X31.all)[2] <- "mean";colnames(X31.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X31.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 31")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 3
X3 <- subset(logMat.GC, select=c("X3"))
X3.sd <- subset(logMatsd2.GC, select=c("X3"))
#merge the stdev with the mean and the chromosome number
X3.all <- cbind(chrms,X3,X3.sd)
#change the column names appropriately
colnames(X3.all)[1] <- "chr"; colnames(X3.all)[2] <- "mean";colnames(X3.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X3.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 3")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 49
X49 <- subset(logMat.GC, select=c("X49"))
X49.sd <- subset(logMatsd2.GC, select=c("X49"))
#merge the stdev with the mean and the chromosome number
X49.all <- cbind(chrms,X49,X49.sd)
#change the column names appropriately
colnames(X49.all)[1] <- "chr"; colnames(X49.all)[2] <- "mean";colnames(X49.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X49.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 49")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 4
X4 <- subset(logMat.GC, select=c("X4"))
X4.sd <- subset(logMatsd2.GC, select=c("X4"))
#merge the stdev with the mean and the chromosome number
X4.all <- cbind(chrms,X4,X4.sd)
#change the column names appropriately
colnames(X4.all)[1] <- "chr"; colnames(X4.all)[2] <- "mean";colnames(X4.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X4.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 4")
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

#line 59
X59 <- subset(logMat.GC, select=c("X59"))
X59.sd <- subset(logMatsd2.GC, select=c("X59"))
#merge the stdev with the mean and the chromosome number
X59.all <- cbind(chrms,X59,X59.sd)
#change the column names appropriately
colnames(X59.all)[1] <- "chr"; colnames(X59.all)[2] <- "mean";colnames(X59.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X59.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 59")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 5
X5 <- subset(logMat.GC, select=c("X5"))
X5.sd <- subset(logMatsd2.GC, select=c("X5"))
#merge the stdev with the mean and the chromosome number
X5.all <- cbind(chrms,X5,X5.sd)
#change the column names appropriately
colnames(X5.all)[1] <- "chr"; colnames(X5.all)[2] <- "mean";colnames(X5.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X5.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 5")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 61
X61 <- subset(logMat.GC, select=c("X61"))
X61.sd <- subset(logMatsd2.GC, select=c("X61"))
#merge the stdev with the mean and the chromosome number
X61.all <- cbind(chrms,X61,X61.sd)
#change the column names appropriately
colnames(X61.all)[1] <- "chr"; colnames(X61.all)[2] <- "mean";colnames(X61.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X61.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 61")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 66
X66 <- subset(logMat.GC, select=c("X66"))
X66.sd <- subset(logMatsd2.GC, select=c("X66"))
#merge the stdev with the mean and the chromosome number
X66.all <- cbind(chrms,X66,X66.sd)
#change the column names appropriately
colnames(X66.all)[1] <- "chr"; colnames(X66.all)[2] <- "mean";colnames(X66.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X66.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 66")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 69
X69 <- subset(logMat.GC, select=c("X69"))
X69.sd <- subset(logMatsd2.GC, select=c("X69"))
#merge the stdev with the mean and the chromosome number
X69.all <- cbind(chrms,X69,X69.sd)
#change the column names appropriately
colnames(X69.all)[1] <- "chr"; colnames(X69.all)[2] <- "mean";colnames(X69.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X69.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 69")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 6
X6 <- subset(logMat.GC, select=c("X6"))
X6.sd <- subset(logMatsd2.GC, select=c("X6"))
#merge the stdev with the mean and the chromosome number
X6.all <- cbind(chrms,X6,X6.sd)
#change the column names appropriately
colnames(X6.all)[1] <- "chr"; colnames(X6.all)[2] <- "mean";colnames(X6.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X6.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 6")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 76
X76 <- subset(logMat.GC, select=c("X76"))
X76.sd <- subset(logMatsd2.GC, select=c("X76"))
#merge the stdev with the mean and the chromosome number
X76.all <- cbind(chrms,X76,X76.sd)
#change the column names appropriately
colnames(X76.all)[1] <- "chr"; colnames(X76.all)[2] <- "mean";colnames(X76.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X76.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 76")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 77
X77 <- subset(logMat.GC, select=c("X77"))
X77.sd <- subset(logMatsd2.GC, select=c("X77"))
#merge the stdev with the mean and the chromosome number
X77.all <- cbind(chrms,X77,X77.sd)
#change the column names appropriately
colnames(X77.all)[1] <- "chr"; colnames(X77.all)[2] <- "mean";colnames(X77.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X77.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 77")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 7
X7 <- subset(logMat.GC, select=c("X7"))
X7.sd <- subset(logMatsd2.GC, select=c("X7"))
#merge the stdev with the mean and the chromosome number
X7.all <- cbind(chrms,X7,X7.sd)
#change the column names appropriately
colnames(X7.all)[1] <- "chr"; colnames(X7.all)[2] <- "mean";colnames(X7.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X7.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 7")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 8
X8 <- subset(logMat.GC, select=c("X8"))
X8.sd <- subset(logMatsd2.GC, select=c("X8"))
#merge the stdev with the mean and the chromosome number
X8.all <- cbind(chrms,X8,X8.sd)
#change the column names appropriately
colnames(X8.all)[1] <- "chr"; colnames(X8.all)[2] <- "mean";colnames(X8.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X8.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 8")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)

#line 9
X9 <- subset(logMat.GC, select=c("X9"))
X9.sd <- subset(logMatsd2.GC, select=c("X9"))
#merge the stdev with the mean and the chromosome number
X9.all <- cbind(chrms,X9,X9.sd)
#change the column names appropriately
colnames(X9.all)[1] <- "chr"; colnames(X9.all)[2] <- "mean";colnames(X9.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X9.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="GC Line 9")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)


#line 1
X001 <- subset(logMat.old, select=c("X001"))
X001.sd <- subset(logMatsd2.old, select=c("X001"))
#merge the stdev with the mean and the chromosome number
X001.all <- cbind(chrms,X001,X001.sd)
#change the column names appropriately
colnames(X001.all)[1] <- "chr"; colnames(X001.all)[2] <- "mean";colnames(X001.all)[3] <- "sd"
#plots the chromosomes with error bars for each one 
with (
  data = X001.all
  , expr = errbar(chr, mean, mean+sd, mean-sd, add=F, pch=16, cap=.015, errbar.col = "mediumvioletred",col="mediumvioletred",xlab="Chromosome",ylab="log2(fold change)",cex.axis=1,las=3,xaxt="n")
)
#adds a line at y=0 
abline(h=0,lty=3)
#adds a title to the graph
title(main="MA Line 1")
##makes it so that the x-axis is all of the chromosomes and none are missing 
labs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1:16)
axis(side=1, at=1:16, labels=labs[0:16],cex.axis=0.6)
############################################################################################################

########### compare to non amplified chromosome ratios
##starting with chromosome I, going in order
##going to need these for all of the chromosomes that have a CNV, but can come back to this later
###probably an easier way to do this, but will come back to this later
#these are the euploid lines 
l002.1 <- getChrmRatio.GC("2",1)
colnames(l002.1) <- c("chr","tracking_id","Line 2")
l003.1 <- getChrmRatio.GC("3",1)
colnames(l003.1) <- c("chr","tracking_id","Line 3")
l005.1 <- getChrmRatio.GC("5",1)
colnames(l005.1) <- c("chr","tracking_id","Line 5")
l006.1 <- getChrmRatio.GC("6",1)
colnames(l006.1) <- c("chr","tracking_id","Line 6")
l008.1 <- getChrmRatio.GC("8",1)
colnames(l008.1) <- c("chr","tracking_id","Line 8")
l009.1 <- getChrmRatio.GC("9",1)
colnames(l009.1) <- c("chr","tracking_id","Line 9")

chrm.1.Ratios <- merge(l001.1,l002.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l003.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l004.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l005.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l006.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l007.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l008.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l009.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l011.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l018.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l021.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l031.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l049.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l059.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l061.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l066.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l069.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l076.1,by=2)
chrm.1.Ratios <- merge(chrm.1.Ratios,l077.1,by=2)
#chrm.1List <- lapply(c(2,3,5,6,8,9)),function(i)) { df <- read.table()}

files <- c(l002.1,l003.1,l005.1,l006.1,l008.1,l009.1)
big.list.of.data.frames <- lapply(files, read.table, header = TRUE)

big.data.frame <- do.call(rbind,big.list.of.data.frames)

str(l002.1)
chrm.1 <- matrix(l002.1,l003.1)
#aneuploid lines 
l011.1 <- getChrmRatio.GC("11",1)
colnames(l011.1) <- c("chr","tracking_id","Line 11")
l018.1 <- getChrmRatio.GC("18",1)
colnames(l018.1) <- c("chr","tracking_id","Line 18")
l001.1 <- getChrmRatio.GC("1",1)
colnames(l001.1) <- c("chr","tracking_id","Line 1")
l021.1 <- getChrmRatio.GC("21",1)
colnames(l021.1) <- c("chr","tracking_id","Line 21")
l031.1 <- getChrmRatio.GC("31",1)
colnames(l031.1) <- c("chr","tracking_id","Line 31")
l049.1 <- getChrmRatio.GC("49",1)
colnames(l049.1) <- c("chr","tracking_id","Line 49")
l004.1 <- getChrmRatio.GC("4",1)
colnames(l004.1) <- c("chr","tracking_id","Line 4")
l059.1 <- getChrmRatio.GC("59",1)
colnames(l059.1) <- c("chr","tracking_id","Line 59")
l061.1 <- getChrmRatio.GC("61",1)
colnames(l061.1) <- c("chr","tracking_id","Line 61")
l066.1 <- getChrmRatio.GC("66",1)
colnames(l066.1) <- c("chr","tracking_id","Line 66")
l069.1 <- getChrmRatio.GC("69",1)
colnames(l069.1) <- c("chr","tracking_id","Line 69")
l076.1 <- getChrmRatio.GC("76",1)
colnames(l076.1) <- c("chr","tracking_id","Line 76")
l077.1 <- getChrmRatio.GC("77",1)
colnames(l077.1) <- c("chr","tracking_id","Line 77")
l007.1 <- getChrmRatio.GC("7",1)
colnames(l007.1) <- c("chr","tracking_id","Line 7")

#loop by line
#for (line in c(2,3,5,6,9)) { getChrmRatio.GC(line,1) }

#although loop by chromosome might be better
#eu.chrm <- vector("list",16)
#for (chrm in 1:16) 
#{ eu.chrm[[chrm]] <- getChrmRatio.GC(2,chrm)}

#chrm1 <- lapply(eu.chrm, "[",1:3)

l002.2 <- getChrmRatio.GC("2",2)
l003.2 <- getChrmRatio.GC("3",2)
l005.2 <- getChrmRatio.GC("5",2)
l006.2 <- getChrmRatio.GC("6",2)
l008.2 <- getChrmRatio.GC("8",2)
l009.2 <- getChrmRatio.GC("9",2)

l011.2 <- getChrmRatio.GC("11",2)
l018.2 <- getChrmRatio.GC("18",2)
l001.2 <- getChrmRatio.GC("1",2)
l021.2 <- getChrmRatio.GC("21",2)
l031.2 <- getChrmRatio.GC("31",2)
l049.2 <- getChrmRatio.GC("49",2)
l004.2 <- getChrmRatio.GC("4",2)
l059.2 <- getChrmRatio.GC("59",2)
l061.2 <- getChrmRatio.GC("61",2)
l066.2 <- getChrmRatio.GC("66",2)
l069.2 <- getChrmRatio.GC("69",2)
l076.2 <- getChrmRatio.GC("76",2)
l077.2 <- getChrmRatio.GC("77",2)
l007.2 <- getChrmRatio.GC("7",2)

l002.3 <- getChrmRatio.GC("2",3)
l003.3 <- getChrmRatio.GC("3",3)
l005.3 <- getChrmRatio.GC("5",3)
l008.3 <- getChrmRatio.GC("8",3)
l006.3 <- getChrmRatio.GC("6",3)
l009.3 <- getChrmRatio.GC("9",3)

l011.3 <- getChrmRatio.GC("11",3)
l018.3 <- getChrmRatio.GC("18",3)
l001.3 <- getChrmRatio.GC("1",3)
l021.3 <- getChrmRatio.GC("21",3)
l031.3 <- getChrmRatio.GC("31",3)
l049.3 <- getChrmRatio.GC("49",3)
l004.3 <- getChrmRatio.GC("4",3)
l059.3 <- getChrmRatio.GC("59",3)
l061.3 <- getChrmRatio.GC("61",3)
l066.3 <- getChrmRatio.GC("66",3)
l069.3 <- getChrmRatio.GC("69",3)
l076.3 <- getChrmRatio.GC("76",3)
l077.3 <- getChrmRatio.GC("77",3)
l007.3 <- getChrmRatio.GC("7",3)

l002.4 <- getChrmRatio.GC("2",4)
l003.4 <- getChrmRatio.GC("3",4)
l005.4 <- getChrmRatio.GC("5",4)
l006.4 <- getChrmRatio.GC("6",4)
l008.4 <- getChrmRatio.GC("8",4)
l009.4 <- getChrmRatio.GC("9",4)

l011.4 <- getChrmRatio.GC("11",4)
l018.4 <- getChrmRatio.GC("18",4)
l001.4 <- getChrmRatio.GC("1",4)
l021.4 <- getChrmRatio.GC("21",4)
l031.4 <- getChrmRatio.GC("31",4)
l049.4 <- getChrmRatio.GC("49",4)
l004.4 <- getChrmRatio.GC("4",4)
l059.4 <- getChrmRatio.GC("59",4)
l061.4 <- getChrmRatio.GC("61",4)
l066.4 <- getChrmRatio.GC("66",4)
l069.4 <- getChrmRatio.GC("69",4)
l076.4 <- getChrmRatio.GC("76",4)
l077.4 <- getChrmRatio.GC("77",4)
l007.4 <- getChrmRatio.GC("7",4)

l002.5 <- getChrmRatio.GC("2",5)
l003.5 <- getChrmRatio.GC("3",5)
l005.5 <- getChrmRatio.GC("5",5)
l006.5 <- getChrmRatio.GC("6",5)
l008.5 <- getChrmRatio.GC("8",5)
l009.5 <- getChrmRatio.GC("9",5)

l011.5 <- getChrmRatio.GC("11",5)
l018.5 <- getChrmRatio.GC("18",5)
l001.5 <- getChrmRatio.GC("1",5)
l021.5 <- getChrmRatio.GC("21",5)
l031.5 <- getChrmRatio.GC("31",5)
l049.5 <- getChrmRatio.GC("49",5)
l004.5 <- getChrmRatio.GC("4",5)
l059.5 <- getChrmRatio.GC("59",5)
l061.5 <- getChrmRatio.GC("61",5)
l066.5 <- getChrmRatio.GC("66",5)
l069.5 <- getChrmRatio.GC("69",5)
l076.5 <- getChrmRatio.GC("76",5)
l077.5 <- getChrmRatio.GC("77",5)
l007.5 <- getChrmRatio.GC("7",5)

l002.6 <- getChrmRatio.GC("2",6)
l003.6 <- getChrmRatio.GC("3",6)
l005.6 <- getChrmRatio.GC("5",6)
l006.6 <- getChrmRatio.GC("6",6)
l008.6 <- getChrmRatio.GC("8",6)
l009.6 <- getChrmRatio.GC("9",6)

l011.6 <- getChrmRatio.GC("11",6)
l018.6 <- getChrmRatio.GC("18",6)
l001.6 <- getChrmRatio.GC("1",6)
l021.6 <- getChrmRatio.GC("21",6)
l031.6 <- getChrmRatio.GC("31",6)
l049.6 <- getChrmRatio.GC("49",6)
l004.6 <- getChrmRatio.GC("4",6)
l059.6 <- getChrmRatio.GC("59",6)
l061.6 <- getChrmRatio.GC("61",6)
l066.6 <- getChrmRatio.GC("66",6)
l069.6 <- getChrmRatio.GC("69",6)
l076.6 <- getChrmRatio.GC("76",6)
l077.6 <- getChrmRatio.GC("77",6)
l007.6 <- getChrmRatio.GC("7",6)

l002.7 <- getChrmRatio.GC("2",7)
l003.7 <- getChrmRatio.GC("3",7)
l005.7 <- getChrmRatio.GC("5",7)
l006.7 <- getChrmRatio.GC("6",7)
l008.7 <- getChrmRatio.GC("8",7)
l009.7 <- getChrmRatio.GC("9",7)

l011.7 <- getChrmRatio.GC("11",7)
l018.7 <- getChrmRatio.GC("18",7)
l001.7 <- getChrmRatio.GC("1",7)
l021.7 <- getChrmRatio.GC("21",7)
l031.7 <- getChrmRatio.GC("31",7)
l049.7 <- getChrmRatio.GC("49",7)
l004.7 <- getChrmRatio.GC("4",7)
l059.7 <- getChrmRatio.GC("59",7)
l061.7 <- getChrmRatio.GC("61",7)
l066.7 <- getChrmRatio.GC("66",7)
l069.7 <- getChrmRatio.GC("69",7)
l076.7 <- getChrmRatio.GC("76",7)
l077.7 <- getChrmRatio.GC("77",7)
l007.7 <- getChrmRatio.GC("7",7)

l002.8 <- getChrmRatio.GC("2",8)
l003.8 <- getChrmRatio.GC("3",8)
l005.8 <- getChrmRatio.GC("5",8)
l006.8 <- getChrmRatio.GC("6",8)
l008.8 <- getChrmRatio.GC("8",8)
l009.8 <- getChrmRatio.GC("9",8)

l011.8 <- getChrmRatio.GC("11",8)
l018.8 <- getChrmRatio.GC("18",8)
l001.8 <- getChrmRatio.GC("1",8)
l021.8 <- getChrmRatio.GC("21",8)
l031.8 <- getChrmRatio.GC("31",8)
l049.8 <- getChrmRatio.GC("49",8)
l004.8 <- getChrmRatio.GC("4",8)
l059.8 <- getChrmRatio.GC("59",8)
l061.8 <- getChrmRatio.GC("61",8)
l066.8 <- getChrmRatio.GC("66",8)
l069.8 <- getChrmRatio.GC("69",8)
l076.8 <- getChrmRatio.GC("76",8)
l077.8 <- getChrmRatio.GC("77",8)
l007.8 <- getChrmRatio.GC("7",8)

l002.9 <- getChrmRatio.GC("2",9)
l003.9 <- getChrmRatio.GC("3",9)
l005.9 <- getChrmRatio.GC("5",9)
l006.9 <- getChrmRatio.GC("6",9)
l008.9 <- getChrmRatio.GC("8",9)
l009.9 <- getChrmRatio.GC("9",9)

l011.9 <- getChrmRatio.GC("11",9)
l018.9 <- getChrmRatio.GC("18",9)
l001.9 <- getChrmRatio.GC("1",9)
l021.9 <- getChrmRatio.GC("21",9)
l031.9 <- getChrmRatio.GC("31",9)
l049.9 <- getChrmRatio.GC("49",9)
l004.9 <- getChrmRatio.GC("4",9)
l059.9 <- getChrmRatio.GC("59",9)
l061.9 <- getChrmRatio.GC("61",9)
l066.9 <- getChrmRatio.GC("66",9)
l069.9 <- getChrmRatio.GC("69",9)
l076.9 <- getChrmRatio.GC("76",9)
l077.9 <- getChrmRatio.GC("77",9)
l007.9 <- getChrmRatio.GC("7",9)

l002.10 <- getChrmRatio.GC("2",10)
l003.10 <- getChrmRatio.GC("3",10)
l005.10 <- getChrmRatio.GC("5",10)
l006.10 <- getChrmRatio.GC("6",10)
l008.10 <- getChrmRatio.GC("8",10)
l009.10 <- getChrmRatio.GC("9",10)

l011.10 <- getChrmRatio.GC("11",10)
l018.10 <- getChrmRatio.GC("18",10)
l001.10 <- getChrmRatio.GC("1",10)
l021.10 <- getChrmRatio.GC("21",10)
l031.10 <- getChrmRatio.GC("31",10)
l049.10 <- getChrmRatio.GC("49",10)
l004.10 <- getChrmRatio.GC("4",10)
l059.10 <- getChrmRatio.GC("59",10)
l061.10 <- getChrmRatio.GC("61",10)
l066.10 <- getChrmRatio.GC("66",10)
l069.10 <- getChrmRatio.GC("69",10)
l076.10 <- getChrmRatio.GC("76",10)
l077.10 <- getChrmRatio.GC("77",10)
l007.10 <- getChrmRatio.GC("7",10)

l002.11 <- getChrmRatio.GC("2",11)
l003.11 <- getChrmRatio.GC("3",11)
l005.11 <- getChrmRatio.GC("5",11)
l006.11 <- getChrmRatio.GC("6",11)
l008.11 <- getChrmRatio.GC("8",11)
l009.11 <- getChrmRatio.GC("9",11)

l011.11 <- getChrmRatio.GC("11",11)
l018.11 <- getChrmRatio.GC("18",11)
l001.11 <- getChrmRatio.GC("1",11)
l021.11 <- getChrmRatio.GC("21",11)
l031.11 <- getChrmRatio.GC("31",11)
l049.11 <- getChrmRatio.GC("49",11)
l004.11 <- getChrmRatio.GC("4",11)
l059.11 <- getChrmRatio.GC("59",11)
l061.11 <- getChrmRatio.GC("61",11)
l066.11 <- getChrmRatio.GC("66",11)
l069.11 <- getChrmRatio.GC("69",11)
l076.11 <- getChrmRatio.GC("76",11)
l077.11 <- getChrmRatio.GC("77",11)
l007.11 <- getChrmRatio.GC("7",11)

l002.12 <- getChrmRatio.GC("2",12)
l003.12 <- getChrmRatio.GC("3",12)
l005.12 <- getChrmRatio.GC("5",12)
l006.12 <- getChrmRatio.GC("6",12)
l008.12 <- getChrmRatio.GC("8",12)
l009.12 <- getChrmRatio.GC("9",12)

l011.12 <- getChrmRatio.GC("11",12)
l018.12 <- getChrmRatio.GC("18",12)
l001.12 <- getChrmRatio.GC("1",12)
l021.12 <- getChrmRatio.GC("21",12)
l031.12 <- getChrmRatio.GC("31",12)
l049.12 <- getChrmRatio.GC("49",12)
l004.12 <- getChrmRatio.GC("4",12)
l059.12 <- getChrmRatio.GC("59",12)
l061.12 <- getChrmRatio.GC("61",12)
l066.12 <- getChrmRatio.GC("66",12)
l069.12 <- getChrmRatio.GC("69",12)
l076.12 <- getChrmRatio.GC("76",12)
l077.12 <- getChrmRatio.GC("77",12)
l007.12 <- getChrmRatio.GC("7",12)

l002.13 <- getChrmRatio.GC("2",13)
l003.13 <- getChrmRatio.GC("3",13)
l005.13 <- getChrmRatio.GC("5",13)
l006.13 <- getChrmRatio.GC("6",13)
l008.13 <- getChrmRatio.GC("8",13)
l009.13 <- getChrmRatio.GC("9",13)

l011.13 <- getChrmRatio.GC("11",13)
l018.13 <- getChrmRatio.GC("18",13)
l001.13 <- getChrmRatio.GC("1",13)
l021.13 <- getChrmRatio.GC("21",13)
l031.13 <- getChrmRatio.GC("31",13)
l049.13 <- getChrmRatio.GC("49",13)
l004.13 <- getChrmRatio.GC("4",13)
l059.13 <- getChrmRatio.GC("59",13)
l061.13 <- getChrmRatio.GC("61",13)
l066.13 <- getChrmRatio.GC("66",13)
l069.13 <- getChrmRatio.GC("69",13)
l076.13 <- getChrmRatio.GC("76",13)
l077.13 <- getChrmRatio.GC("77",13)
l007.13 <- getChrmRatio.GC("7",13)

l002.14 <- getChrmRatio.GC("2",14)
l003.14 <- getChrmRatio.GC("3",14)
l005.14 <- getChrmRatio.GC("5",14)
l006.14 <- getChrmRatio.GC("6",14)
l008.14 <- getChrmRatio.GC("8",14)
l009.14 <- getChrmRatio.GC("9",14)

l011.14 <- getChrmRatio.GC("11",14)
l018.14 <- getChrmRatio.GC("18",14)
l001.14 <- getChrmRatio.GC("1",14)
l021.14 <- getChrmRatio.GC("21",14)
l031.14 <- getChrmRatio.GC("31",14)
l049.14 <- getChrmRatio.GC("49",14)
l004.14 <- getChrmRatio.GC("4",14)
l059.14 <- getChrmRatio.GC("59",14)
l061.14 <- getChrmRatio.GC("61",14)
l066.14 <- getChrmRatio.GC("66",14)
l069.14 <- getChrmRatio.GC("69",14)
l076.14 <- getChrmRatio.GC("76",14)
l077.14 <- getChrmRatio.GC("77",14)
l007.14 <- getChrmRatio.GC("7",14)

l002.15 <- getChrmRatio.GC("2",15)
l003.15 <- getChrmRatio.GC("3",15)
l005.15 <- getChrmRatio.GC("5",15)
l006.15 <- getChrmRatio.GC("6",15)
l008.15 <- getChrmRatio.GC("8",15)
l009.15 <- getChrmRatio.GC("9",15)

l011.15 <- getChrmRatio.GC("11",15)
l018.15 <- getChrmRatio.GC("18",15)
l001.15 <- getChrmRatio.GC("1",15)
l021.15 <- getChrmRatio.GC("21",15)
l031.15 <- getChrmRatio.GC("31",15)
l049.15 <- getChrmRatio.GC("49",15)
l004.15 <- getChrmRatio.GC("4",15)
l059.15 <- getChrmRatio.GC("59",15)
l061.15 <- getChrmRatio.GC("61",15)
l066.15 <- getChrmRatio.GC("66",15)
l069.15 <- getChrmRatio.GC("69",15)
l076.15 <- getChrmRatio.GC("76",15)
l077.15 <- getChrmRatio.GC("77",15)
l007.15 <- getChrmRatio.GC("7",15)

l002.16 <- getChrmRatio.GC("2",16)
l003.16 <- getChrmRatio.GC("3",16)
l005.16 <- getChrmRatio.GC("5",16)
l006.16 <- getChrmRatio.GC("6",16)
l008.16 <- getChrmRatio.GC("8",16)
l009.16 <- getChrmRatio.GC("9",16)

l011.16 <- getChrmRatio.GC("11",16)
l018.16 <- getChrmRatio.GC("18",16)
l001.16 <- getChrmRatio.GC("1",16)
l021.16 <- getChrmRatio.GC("21",16)
l031.16 <- getChrmRatio.GC("31",16)
l049.16 <- getChrmRatio.GC("49",16)
l004.16 <- getChrmRatio.GC("4",16)
l059.16 <- getChrmRatio.GC("59",16)
l061.16 <- getChrmRatio.GC("61",16)
l066.16 <- getChrmRatio.GC("66",16)
l069.16 <- getChrmRatio.GC("69",16)
l076.16 <- getChrmRatio.GC("76",16)
l077.16 <- getChrmRatio.GC("77",16)
l007.16 <- getChrmRatio.GC("7",16)

#this makes a matrix with the log2 ratio between ancestor and sample of each gene 

getChrmRatio.GC.test <- function(line, chr){ #input line as "XXX", chr as a number
  rat.GC <- getFPKMRatio.GC(line) #get list of ratios which includes chrm and gene
  rat.cs.GC <- rat.GC[order(rat.GC[,3]),]#order by column 3
  return(rat.cs.GC)
}

yl.GC=NULL
for(line in rep.GC){
  y.GC <- c(yl.GC,log2(getChrmRatio.GC.test(line,chrm)[,3]))
}

warnings()
logMat.GC.test<- matrix(y.GC,nrow=5431,ncol=20)
colnames(logMat.GC.test) <- rep.GC
logMat.GC.test <- data.frame(logMat.GC.test)
genes <- ancAvg.GC[,2]
chrms <- ancAvg.GC[,1]
logMat2.GC.test <-cbind(genes,logMat.GC.test)
logMat2.GC.test <-cbind(chrms,logMat2.GC.test)
######################################################################
#testing to see if I specify the chromosome it'll work and then I can just make like 16 of these matrices and specify the column instead of the file from above

getChrmRatio.GC.test <- function(line, chr){ #input line as "XXX", chr as a number
  rat.GC <- getFPKMRatio.GC(line) #get list of ratios which includes chrm and gene
  rat.cs.GC <- rat.GC[order(rat.GC[,3]),]#order by column 3
  return(rat.cs.GC)
}



yl.GC.1=NULL
for(line in 1){
  for(chr in 1:16) {
  y.GC.1 <- c(yl.GC.1,log2(getChrmRatio.GC.test(line,chr)[,3]))
} }


warnings()
logMat.GC.test.1<- matrix(y.GC.1,nrow=5431,ncol=20)
colnames(logMat.GC.test.1) <- rep.GC
logMat.GC.test <- data.frame(logMat.GC.test)
genes <- ancAvg.GC[,2]
chrms <- ancAvg.GC[,1]
logMat2.GC.test <-cbind(genes,logMat.GC.test)
logMat2.GC.test <-cbind(chrms,logMat2.GC.test)

#testing out tapply
#test.1.11 <- tapply(logMat2.GC.test$chrms, logMat2.GC.test$X11, getChrmRatio.GC.test)
#nope not what I want it to do
#trying to select the columns and rows I want with the bracket method
#test.1.11.2 <- logMat2.GC.test["1","X11"]
#test.1.11 <- subset(logMat2.GC.test,chrms="1")

##THIS ONE WORKS
test.1.11 <- logMat2.GC.test[logMat2.GC.test$chrms == 1,]

#logMAt2.GC.test = merge(logMat2.GC.test,ancAvg.GC,by.x="genes",by.y="tracking_id")

#this works, just gives an interesting curve 
#plot((log2(l002.1[,3])),ylab="log2(fold change)", xlab="Chromosome")
library(ggplot2)
#autoplot((log2(l002.1[,3])),ylab="log2(fold change)", xlab="Chromosome")

#ancestor ratios
ancAvg <- rowMeans(sc_anc[,3:6],na.rm=T)
ancAvg <- cbind(sc_anc[,1:2],ancAvg)

##############################################################################################

#CHROMOSOME 1 

#RATIOS FOR CHROMOSOME 1 
anRatio.c1 <- c(l152.1[,3],l018.1[,3],l021.1[,3],l007.1[,3])
euRatio.c1 <- c(l002.1[,3],l003.1[,3],l005.1[,3],l006.1[,3],l009.1[,3])

euRatio.c1.test<- c(logMat2.GC.test[logMat2.GC.test$chrms == 1,logMat2.GC.test$X2],logMat2.GC.test[logMat2.GC.test$chrms == 1,logMat2.GC.test$X3],logMat2.GC.test[logMat2.GC.test$chrms == 1,logMat2.GC.test$X5],logMat2.GC.test[logMat2.GC.test$chrms == 1,logMat2.GC.test$X6],logMat2.GC.test[logMat2.GC.test$chrms == 1,logMat2.GC.test$X8],logMat2.GC.test[logMat2.GC.test$chrms == 1,logMat2.GC.test$X9])

euploids <- c(2,3,5,6,8,9)
for (line in euploids){
  for (chrms in 1:16){
    line.chrms <- logMat2.GC.test[logMat2.GC.test$chrms == chrms ,logMat2.GC.test$Xline]
  }
}


aneuc1 <- cbind(anRatio.c1,euRatio.c1)
aneuc2 <- data.frame(aneuc1)
anovac1 <- aov(anRatio.c1~euRatio.c1, data=aneuc2)
summary(anovac1)

###test if data is normally distributed
shapiro.test(aneuc2$anRatio.c1)
resids <- resid(anovac1)
hist(resids, col="pink", breaks=10, xlab="residuals")
shapiro.test(resids)
##reject H0, failed the test
##boxplot and t-test for chromosome 1 trisomic
boxplot((log2(anRatio.c1)), (log2(euRatio.c1)),names=c("Trisomic", "Euploid"),ylab="log2(fold change)", col=c("lightblue", "purple"), main="Chr 1 Euploid vs Aneuploid")
t.test((log2(anRatio.c1)), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
wilcox.test((log2(anRatio.c1)), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
##p value wilcox test: 3.369e-5 
##adds a line at the 0 point on the y axis to delimit where there would be no fold change 
##dashed line (lty=3), can find more types of lines and information under graphical parameters in the help
abline(h=0,lty=3)
##want to add text to graph that states the p-value of the difference between trisomic and euploid lines 
p <- 0.0002234
mylabel = bquote(italic(p) == .(format(p, digits = 9)))
text(x=2.3,y = 3.5, labels = mylabel)
##
###############################################################################################
###############################################################################################
##BOXPLOTS PER LINE AND CHROMOSOME
##boxplot and t-test for chromosome 1 monosomic, tri, disomic 
boxplot((log2(l011.1[,3])), (log2(euRatio.c1)),(log2(anRatio.c1)),names=c("Monosomic", "Disomic","Trisomic"),ylab="log2(fold change)", col=c("turquoise", "hotpink","purple"), main="Chr 1 Relative Gene Expression Levels")
abline(h=0,lty=3)
t.test((log2(l011.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
##p value for monosomic: 2.841x10^-9
wilcox.test((log2(l011.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
##p value wilcox tets: <2.2e-16
p <- 0.0002234
mylabel = bquote(italic(p) == .(format(p, digits = 9)))
##text(x=3.2,y = 3.7, labels = mylabel,cex=.8)
text(x=1, y=0.5, "*", pos=3, cex=1.5, col="red")
#text(x=3, y=4, "*", pos=3, cex=1.2)
mtext("*", side=3, line=0, at=3, cex=1.5, col="red")
##going to need these for all of the chromosomes that have a CNV, but can come back to this later

##############################################################################################

#CHROMOSOME 1 
par(mar=c(5, 4, 4, 2) + 0.1)
par(mfrow=c(1,1))
#going to do this for all lines for chromosome 1 
#try to plot them all on the same graph? 

#euRatio.test <- c(logMat2.GC.test[,7],logMat2.GC.test[,7])

boxplot((log2(l001.1[,3])),(log2(l002.1[,3])),(log2(l003.1[,3])),(log2(l004.1[,3])),(log2(l005.1[,3])),
        (log2(l006.1[,3])),(log2(l007.1[,3])),(log2(l008.1[,3])),(log2(l009.1[,3])),(log2(l011.1[,3])),
        (log2(l018.1[,3])),(log2(l021.1[,3])),(log2(l031.1[,3])),(log2(l049.1[,3])),(log2(l059.1[,3])),
        (log2(l061.1[,3])),(log2(l066.1[,3])),(log2(l069.1[,3])),(log2(l076.1[,3])),(log2(l077.1[,3]))
        ,names=c( "1","2","3","4","5","6","7","8","9","11","18","21","31","49","59","61","66","69","76","77"),
        ylab="log2(fold change)", col=c("green", "green","green","cyan","green","green","deeppink1","cyan",
                                        "green",
                        "deeppink1","deeppink1","deeppink1","green","cyan","cyan","cyan","cyan",
                        "green","cyan","cyan","cyan"),main="GC Lines Chromsome 1",las=3)
euRatio.c1 <- c(l002.1[,3],l003.1[,3],l005.1[,3],l006.1[,3],l009.1[,3],l069.1[,3])
#use non-parametric because these are not normally distributed
#run same test on all lines, even those that are not "supposed" to be aneuploid
t.test((log2(l001.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
t.test((log2(l002.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
t.test((log2(l003.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
t.test((log2(l004.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
t.test((log2(l005.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
t.test((log2(l006.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
wilcox.test((log2(l007.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
t.test((log2(l008.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
t.test((log2(l009.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
wilcox.test((log2(l011.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
wilcox.test((log2(l018.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
wilcox.test((log2(l021.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
t.test((log2(l031.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
t.test((log2(l049.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
#0.0468 weird. maybe because it's non parametric?
t.test((log2(l059.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
t.test((log2(l061.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
#I think so, because now the p-value is 0.4
t.test((log2(l066.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
t.test((log2(l069.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
t.test((log2(l076.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
t.test((log2(l077.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)

abline(h=0,lty=3)
#p3 <- 0.0006221
#mylabel = bquote(italic(p) == .(format(p3, digits = 9)))
#text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=7, y=1, "*", pos=3, cex=1.5, col="red")
text(x=10, y=.5, "*", pos=3, cex=1.5, col="red")
text(x=11, y=1.2, "*", pos=3, cex=1.5, col="red")
text(x=12, y=1.2, "*", pos=3, cex=1.5, col="red")

#############################################################################################



##############################################################################################
#CHROMOSOME 2

#RATIOS FOR CHROMOSOME 2

euRatio.c2 <- c(l002.2[,3],l003.2[,3],l005.2[,3],l006.2[,3],l008.2[,3],l009.2[,3],l069.2[,3])


boxplot((log2(l001.2[,3])),(log2(l002.2[,3])),(log2(l003.2[,3])),(log2(l004.2[,3])),
        (log2(l005.2[,3])),(log2(l006.2[,3])),(log2(l007.2[,3])),(log2(l008.2[,3])),
        (log2(l009.2[,3])),(log2(l011.2[,3])),(log2(l018.2[,3])),(log2(l021.2[,3])),
        (log2(l031.2[,3])),(log2(l049.2[,3])),(log2(l059.2[,3])),(log2(l061.2[,3])),
        (log2(l066.2[,3])),(log2(l069.2[,3])),(log2(l076.2[,3])),(log2(l077.2[,3])),
        names=c("1","2","3","4","5","6","7","8","9","11","18","21","31","49","59","61","66","69","76","77"),
        ylab="log2(fold change)", 
        col=c("green", "green","green","cyan","green","green","cyan","cyan","green","cyan","cyan","cyan",
              "green","cyan","cyan","cyan","cyan","green","cyan","cyan"),main="GC Lines Chromsome 2",las=3)

#use non-parametric because these are not normally distributed
#run same test on all lines, even those that are not "supposed" to be aneuploid
t.test((log2(l001.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l002.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l003.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l004.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l005.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l006.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l007.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l008.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l009.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l011.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
#p-value 0.05396
t.test((log2(l018.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l021.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
#p-value  0.03724
wilcox.test((log2(l021.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
#p-value = 0.003565 
t.test((log2(l031.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l049.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l059.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l061.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l066.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l069.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l076.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)
t.test((log2(l077.2[,3])), (log2(euRatio.c2)), paired=FALSE, var.equal=TRUE)

abline(h=0,lty=3)
#p3 <- 0.0006221
#mylabel = bquote(italic(p) == .(format(p3, digits = 9)))
#text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
#text(x=4, y=2.5, "*", pos=3, cex=1.5, col="red")
#text(x=14, y=2.5, "*", pos=3, cex=1.5, col="red")
###############################################################################################
#CHROMOSOME 3

#RATIOS FOR CHROMOSOME 3

euRatio.c3 <- c(l002.3[,3],l003.3[,3],l005.3[,3],l006.3[,3],l008.3[,3],l009.3[,3],l069.3[,3])


boxplot((log2(l001.3[,3])),(log2(l002.3[,3])),(log2(l003.3[,3])),(log2(l004.3[,3])),(log2(l005.3[,3])),
        (log2(l006.3[,3])),(log2(l007.3[,3])),(log2(l008.3[,3])),(log2(l009.3[,3])),(log2(l011.3[,3])),
        (log2(l018.3[,3])),(log2(l021.3[,3])),(log2(l031.3[,3])),(log2(l049.3[,3])),(log2(l059.3[,3])),
        (log2(l061.3[,3])),(log2(l066.3[,3])),(log2(l069.3[,3])),(log2(l076.3[,3])),(log2(l077.3[,3])),
        names=c( "1","2","3","4","5","6","7","8","9","11","18","21","31","49","59","61","66","69","76","77"),
        ylab="log2(fold change)", col=c("green", "green","green","cyan","green","green","cyan","cyan",
      "green","cyan","cyan","cyan","green","cyan","cyan","cyan","cyan","green","cyan","cyan"),
      main="GC Lines Chromsome 3",las=3)

#use non-parametric because these are not normally distributed
#run same test on all lines, even those that are not "supposed" to be aneuploid
t.test((log2(l001.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
t.test((log2(l002.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
t.test((log2(l003.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
t.test((log2(l004.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
t.test((log2(l005.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
t.test((log2(l006.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
t.test((log2(l007.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
t.test((log2(l008.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
t.test((log2(l009.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
t.test((log2(l011.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
#p=0.01514
t.test((log2(l018.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
#p-value = 0.007955
t.test((log2(l021.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)

t.test((log2(l031.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
#p-value = 0.04687
t.test((log2(l049.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
t.test((log2(l059.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
t.test((log2(l061.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
#p-value = 0.0404
t.test((log2(l066.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
t.test((log2(l069.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
t.test((log2(l076.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
#p-value = 0.0006532
t.test((log2(l077.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)

abline(h=0,lty=3)
#p3 <- 0.0006221
#mylabel = bquote(italic(p) == .(format(p3, digits = 9)))
#text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
#text(x=4, y=2.5, "*", pos=3, cex=1.5, col="red")
#text(x=14, y=2.5, "*", pos=3, cex=1.5, col="red")

###############################################################################################
#CHROMOSOME 4

#RATIOS FOR CHROMOSOME 4

euRatio.c4 <- c(l002.4[,3],l003.4[,3],l005.4[,3],l006.4[,3],l008.4[,3],l009.4[,3],l069.4[,3])


boxplot((log2(l001.4[,3])),(log2(l002.4[,3])),(log2(l003.4[,3])),(log2(l004.4[,3])),(log2(l005.4[,3])),
        (log2(l006.4[,3])),(log2(l007.4[,3])),(log2(l008.4[,3])),(log2(l009.4[,3])),(log2(l011.4[,3])),
        (log2(l018.4[,3])),(log2(l021.4[,3])),(log2(l031.4[,3])),(log2(l049.4[,3])),(log2(l059.4[,3])),
        (log2(l061.4[,3])),(log2(l066.4[,3])),(log2(l069.4[,3])),(log2(l076.4[,3])),(log2(l077.4[,3])),
        names=c( "1","2","3","4","5","6","7","8","9","11","18","21","31","49","59","61","66","69","76","77"),
        ylab="log2(fold change)", col=c("green", "green","green","cyan","green","green","cyan","cyan",
                                        "green","cyan","cyan","cyan","green","cyan","cyan","cyan","cyan","green","cyan","cyan"),
        main="GC Lines Chromsome 4",las=3)
#use non-parametric because these are not normally distributed
#run same test on all lines, even those that are not "supposed" to be aneuploid
t.test((log2(l001.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
t.test((log2(l002.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
t.test((log2(l003.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
t.test((log2(l004.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
t.test((log2(l005.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
t.test((log2(l006.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
t.test((log2(l007.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
t.test((log2(l008.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
t.test((log2(l009.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
t.test((log2(l011.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
#p=0.01514
t.test((log2(l018.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
#p-value = 0.007955
t.test((log2(l021.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
t.test((log2(l031.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
#p-value = 0.04687
t.test((log2(l049.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
t.test((log2(l059.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
t.test((log2(l061.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
#p-value = 0.0404
t.test((log2(l066.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
t.test((log2(l069.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
t.test((log2(l076.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)
#p-value = 0.0006532
t.test((log2(l077.4[,3])), (log2(euRatio.c4)), paired=FALSE, var.equal=TRUE)

abline(h=0,lty=3)
#p3 <- 0.0006221
#mylabel = bquote(italic(p) == .(format(p3, digits = 9)))
#text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
#text(x=4, y=2.5, "*", pos=3, cex=1.5, col="red")
#text(x=14, y=2.5, "*", pos=3, cex=1.5, col="red")
###############################################################################################
#CHROMOSOME 5 

#RATIOS FOR CHROMOSOME 5

euRatio.c5 <- c(l002.5[,3],l003.5[,3],l005.5[,3],l006.5[,3],l008.5[,3],l009.5[,3],l069.5[,3])


boxplot((log2(l001.5[,3])),(log2(l002.5[,3])),(log2(l003.5[,3])),(log2(l004.5[,3])),(log2(l005.5[,3])),
        (log2(l006.5[,3])),(log2(l007.5[,3])),(log2(l008.5[,3])),(log2(l009.5[,3])),(log2(l011.5[,3])),
        (log2(l018.5[,3])),(log2(l021.5[,3])),(log2(l031.5[,3])),(log2(l049.5[,3])),(log2(l059.5[,3])),
        (log2(l061.5[,3])),(log2(l066.5[,3])),(log2(l069.5[,3])),(log2(l076.5[,3])),(log2(l077.5[,3])),
        names=c( "1","2","3","4","5","6","7","8","9","11","18","21","31","49","59","61","66","69","76","77"),
        ylab="log2(fold change)", col=c("green", "green","green","deeppink1","green","green","cyan","cyan",
                  "green","cyan","cyan","cyan","green","deeppink1","cyan","cyan","cyan","green","cyan","cyan"),
        main="GC Lines Chromsome 5",las=3)
#use non-parametric because these are not normally distributed
#run same test on all lines, even those that are not "supposed" to be aneuploid
t.test((log2(l001.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l002.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l003.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
wilcox.test((log2(l004.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
t.test((log2(l005.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l006.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l007.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l008.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l009.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l011.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l018.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l021.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l031.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
wilcox.test((log2(l049.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
t.test((log2(l059.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l061.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l066.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l069.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l076.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l077.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)

abline(h=0,lty=3)
#p3 <- 0.0006221
#mylabel = bquote(italic(p) == .(format(p3, digits = 9)))
#text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=4, y=2.5, "*", pos=3, cex=1.5, col="red")
text(x=14, y=2.5, "*", pos=3, cex=1.5, col="red")


########################################################################################
#CHROMOSOME 6 

#RATIOS FOR CHROMOSOME 6

euRatio.c6 <- c(l002.6[,3],l003.6[,3],l005.6[,3],l006.6[,3],l008.6[,3],l009.6[,3],l069.6[,3])


boxplot((log2(l001.6[,3])),(log2(l002.6[,3])),(log2(l003.6[,3])),(log2(l004.6[,3])),(log2(l005.6[,3])),
        (log2(l006.6[,3])),(log2(l007.6[,3])),(log2(l008.6[,3])),(log2(l009.6[,3])),(log2(l011.6[,3])),
        (log2(l018.6[,3])),(log2(l021.6[,3])),(log2(l031.6[,3])),(log2(l049.6[,3])),(log2(l059.6[,3])),
        (log2(l061.6[,3])),(log2(l066.6[,3])),(log2(l069.6[,3])),(log2(l076.6[,3])),(log2(l077.6[,3])),
        names=c( "1","2","3","4","5","6","7","8","9","11","18","21","31","49","59","61","66","69","76","77"),
        ylab="log2(fold change)", col=c("green", "green","green","cyan","green","green","cyan","cyan",
                                        "green","cyan","cyan","cyan","green","cyan","cyan","cyan","cyan","green","cyan","cyan"),
        main="GC Lines Chromsome 6",las=3)
#use non-parametric because these are not normally distributed
#run same test on all lines, even those that are not "supposed" to be aneuploid
t.test((log2(l001.6[,3])), (log2(euRatio.c6)), paired=FALSE, var.equal=TRUE)
t.test((log2(l002.6[,3])), (log2(euRatio.c6)), paired=FALSE, var.equal=TRUE)
t.test((log2(l003.6[,3])), (log2(euRatio.c6)), paired=FALSE, var.equal=TRUE)
t.test((log2(l004.6[,3])), (log2(euRatio.c6)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
t.test((log2(l005.6[,3])), (log2(euRatio.c6)), paired=FALSE, var.equal=TRUE)
t.test((log2(l006.6[,3])), (log2(euRatio.c6)), paired=FALSE, var.equal=TRUE)
t.test((log2(l007.6[,3])), (log2(euRatio.c6)), paired=FALSE, var.equal=TRUE)
t.test((log2(l008.6[,3])), (log2(euRatio.c6)), paired=FALSE, var.equal=TRUE)
t.test((log2(l009.6[,3])), (log2(euRatio.c6)), paired=FALSE, var.equal=TRUE)
t.test((log2(l011.6[,3])), (log2(euRatio.c6)), paired=FALSE, var.equal=TRUE)
t.test((log2(l018.6[,3])), (log2(euRatio.c6)), paired=FALSE, var.equal=TRUE)
t.test((log2(l021.6[,3])), (log2(euRatio.c6)), paired=FALSE, var.equal=TRUE)
t.test((log2(l031.6[,3])), (log2(euRatio.c6)), paired=FALSE, var.equal=TRUE)
t.test((log2(l049.6[,3])), (log2(euRatio.c6)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
t.test((log2(l059.6[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l061.6[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l066.6[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l069.6[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l076.6[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
t.test((log2(l077.6[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)

abline(h=0,lty=3)
#p3 <- 0.0006221
#mylabel = bquote(italic(p) == .(format(p3, digits = 9)))
#text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
#text(x=4, y=2.5, "*", pos=3, cex=1.5, col="red")
#text(x=14, y=2.5, "*", pos=3, cex=1.5, col="red")
################################################################################


#CHROMOSOME 7 

#RATIOS FOR CHROMOSOME 7
euRatio.c7 <- c(l002.7[,3],l003.7[,3],l005.7[,3],l006.7[,3],l008.7[,3],l009.7[,3],l069.7[,3])


boxplot((log2(l001.7[,3])),(log2(l002.7[,3])),(log2(l003.7[,3])),(log2(l004.7[,3])),(log2(l005.7[,3])),
        (log2(l006.7[,3])),(log2(l007.7[,3])),(log2(l008.7[,3])),(log2(l009.7[,3])),(log2(l011.7[,3])),
        (log2(l018.7[,3])),(log2(l021.7[,3])),(log2(l031.7[,3])),(log2(l049.7[,3])),(log2(l059.7[,3])),
        (log2(l061.7[,3])),(log2(l066.7[,3])),(log2(l069.7[,3])),(log2(l076.7[,3])),(log2(l077.7[,3])),
        names=c( "1","2","3","4","5","6","7","8","9","11","18","21","31","49","59","61","66","69","76","77"),
        ylab="log2(fold change)", col=c("green", "green","green","cyan","green","green","cyan","cyan",
              "green","cyan","cyan","cyan","green","cyan","deeppink1","deeppink1","deeppink1","green","cyan","cyan"),
        main="GC Lines Chromsome 7",las=3)
#use non-parametric because these are not normally distributed
#run same test on all lines, even those that are not "supposed" to be aneuploid
t.test((log2(l001.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
t.test((log2(l002.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
t.test((log2(l003.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
t.test((log2(l004.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l005.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
t.test((log2(l006.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
t.test((log2(l007.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
t.test((log2(l008.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
t.test((log2(l009.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
t.test((log2(l011.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
t.test((log2(l018.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
t.test((log2(l021.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
wilcox.test((log2(l031.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
#clearly not trisomic for 31
t.test((log2(l049.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
wilcox.test((log2(l059.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
wilcox.test((log2(l061.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
wilcox.test((log2(l066.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
t.test((log2(l069.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
t.test((log2(l076.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
t.test((log2(l077.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)

abline(h=0,lty=3)
#p3 <- 0.0006221
#mylabel = bquote(italic(p) == .(format(p3, digits = 9)))
#text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=15, y=3, "*", pos=3, cex=1.5, col="red")
text(x=16, y=3, "*", pos=3, cex=1.5, col="red")
text(x=17, y=3, "*", pos=3, cex=1.5, col="red")
###############################################################################################

#CHROMOSOME 8 

#RATIOS FOR CHROMOSOME 8

euRatio.c8 <- c(l002.8[,3],l003.8[,3],l005.8[,3],l006.8[,3],l008.8[,3],l009.8[,3],l069.8[,3])


boxplot((log2(l001.8[,3])),(log2(l002.8[,3])),(log2(l003.8[,3])),(log2(l004.8[,3])),(log2(l005.8[,3])),
        (log2(l006.8[,3])),(log2(l007.8[,3])),(log2(l008.8[,3])),(log2(l009.8[,3])),(log2(l011.8[,3])),
        (log2(l018.8[,3])),(log2(l021.8[,3])),(log2(l031.8[,3])),(log2(l049.8[,3])),(log2(l059.8[,3])),
        (log2(l061.8[,3])),(log2(l066.8[,3])),(log2(l069.8[,3])),(log2(l076.8[,3])),(log2(l077.8[,3])),
        names=c( "1","2","3","4","5","6","7","8","9","11","18","21","31","49","59","61","66","69","76","77"),
        ylab="log2(fold change)", col=c("green", "green","green","cyan","green","green","cyan","cyan",
                                        "green","cyan","cyan","cyan","green","cyan","cyan","cyan","cyan","green","cyan","cyan"),
        main="GC Lines Chromsome 8",las=3)
#use non-parametric because these are not normally distributed
#run same test on all lines, even those that are not "supposed" to be aneuploid
t.test((log2(l001.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
t.test((log2(l002.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
t.test((log2(l003.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
t.test((log2(l004.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l005.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
t.test((log2(l006.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
t.test((log2(l007.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
t.test((log2(l008.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
t.test((log2(l009.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
t.test((log2(l011.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
t.test((log2(l018.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
t.test((log2(l021.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
t.test((log2(l031.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
#clearly not trisomic for 3
t.test((log2(l049.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
t.test((log2(l059.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-
t.test((log2(l061.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
t.test((log2(l066.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l069.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
t.test((log2(l076.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
t.test((log2(l077.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)

abline(h=0,lty=3)

###############################################################################################

#CHROMOSOME 8 

#RATIOS FOR CHROMOSOME 8

euRatio.c9 <- c(l002.9[,3],l003.9[,3],l005.9[,3],l006.9[,3],l008.9[,3],l009.9[,3],l069.9[,3])


boxplot((log2(l001.9[,3])),(log2(l002.9[,3])),(log2(l003.9[,3])),(log2(l004.9[,3])),(log2(l005.9[,3])),
        (log2(l006.9[,3])),(log2(l007.9[,3])),(log2(l008.9[,3])),(log2(l009.9[,3])),(log2(l011.9[,3])),
        (log2(l018.9[,3])),(log2(l021.9[,3])),(log2(l031.9[,3])),(log2(l049.9[,3])),(log2(l059.9[,3])),
        (log2(l061.9[,3])),(log2(l066.9[,3])),(log2(l069.9[,3])),(log2(l076.9[,3])),(log2(l077.9[,3])),
        names=c( "1","2","3","4","5","6","7","8","9","11","18","21","31","49","59","61","66","69","76","77"),
        ylab="log2(fold change)", col=c("green", "green","green","cyan","green","green","cyan","cyan",
               "green","cyan","cyan","cyan","green","cyan","cyan","cyan","cyan","green","deeppink1","cyan"),
        main="GC Lines Chromsome 9",las=3)
#use non-parametric because these are not normally distributed
#run same test on all lines, even those that are not "supposed" to be aneuploid
t.test((log2(l001.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
t.test((log2(l002.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
t.test((log2(l003.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
t.test((log2(l004.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l005.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
t.test((log2(l006.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
t.test((log2(l007.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
t.test((log2(l008.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
t.test((log2(l009.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
t.test((log2(l011.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
t.test((log2(l018.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
t.test((log2(l021.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
t.test((log2(l031.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
#clearly not trisoic for 3
t.test((log2(l049.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l059.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-
t.test((log2(l061.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
t.test((log2(l066.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l069.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
wilcox.test((log2(l076.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
t.test((log2(l077.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)

abline(h=0,lty=3)

text(x=19, y=2.2, "*", pos=3, cex=1.5, col="red")


###############################################################################################

#CHROMOSOME 10 

#RATIOS FOR CHROMOSOME 10

euRatio.c10 <- c(l002.10[,3],l003.10[,3],l005.10[,3],l006.10[,3],l008.10[,3],l009.10[,3],l069.10[,3])


boxplot((log2(l001.10[,3])),(log2(l002.10[,3])),(log2(l003.10[,3])),(log2(l004.10[,3])),(log2(l005.10[,3])),
        (log2(l006.10[,3])),(log2(l007.10[,3])),(log2(l008.10[,3])),(log2(l009.10[,3])),(log2(l011.10[,3])),
        (log2(l018.10[,3])),(log2(l021.10[,3])),(log2(l031.10[,3])),(log2(l049.10[,3])),(log2(l059.10[,3])),
        (log2(l061.10[,3])),(log2(l066.10[,3])),(log2(l069.10[,3])),(log2(l076.10[,3])),(log2(l077.10[,3])),
        names=c( "1","2","3","4","5","6","7","8","9","11","18","21","31","49","59","61","66","69","76","77"),
        ylab="log2(fold change)", col=c("green", "green","green","cyan","green","green","cyan","cyan",
       "green","cyan","cyan","cyan","green","cyan","cyan","cyan","cyan","green","cyan","cyan"),
        main="GC Lines Chromsome 10",las=3)
#use non-parametric because these are not normally distributed
#run same test on all lines, even those that are not "supposed" to be aneuploid
t.test((log2(l001.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
t.test((log2(l002.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
t.test((log2(l003.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
t.test((log2(l004.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l005.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
t.test((log2(l006.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
t.test((log2(l007.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
t.test((log2(l008.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
t.test((log2(l009.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
t.test((log2(l011.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
t.test((log2(l018.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
t.test((log2(l021.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
t.test((log2(l031.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
#clearly not trisoic for 3
t.test((log2(l049.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l059.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-
t.test((log2(l061.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
t.test((log2(l066.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l069.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
t.test((log2(l076.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
#p-value = 9.705e-11
wilcox.test((log2(l076.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
#p-value = 1.79e-15

t.test((log2(l077.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)

abline(h=0,lty=3)

text(x=19, y=2.2, "*", pos=3, cex=1.5, col="red")

###############################################################################################

#CHROMOSOME 11

#RATIOS FOR CHROMOSOME 11

euRatio.c11 <- c(l002.11[,3],l003.11[,3],l005.11[,3],l006.11[,3],l008.11[,3],l009.11[,3],l069.11[,3])


boxplot((log2(l001.11[,3])),(log2(l002.11[,3])),(log2(l003.11[,3])),(log2(l004.11[,3])),(log2(l005.11[,3])),
        (log2(l006.11[,3])),(log2(l007.11[,3])),(log2(l008.11[,3])),(log2(l009.11[,3])),(log2(l011.11[,3])),
        (log2(l018.11[,3])),(log2(l021.11[,3])),(log2(l031.11[,3])),(log2(l049.11[,3])),(log2(l059.11[,3])),
        (log2(l061.11[,3])),(log2(l066.11[,3])),(log2(l069.11[,3])),(log2(l076.11[,3])),(log2(l077.11[,3])),
        names=c( "1","2","3","4","5","6","7","8","9","11","18","21","31","49","59","61","66","69","76","77"),
        ylab="log2(fold change)", col=c("green", "green","green","cyan","green","green","cyan","cyan",
           "green","cyan","cyan","cyan","green","cyan","cyan","cyan","cyan","green","cyan","cyan"),
        main="GC Lines Chromsome 11",las=3)
#use non-parametric because these are not normally distributed
#run same test on all lines, even those that are not "supposed" to be aneuploid
t.test((log2(l001.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
t.test((log2(l002.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
t.test((log2(l003.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
t.test((log2(l004.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l005.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
t.test((log2(l006.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
t.test((log2(l007.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
t.test((log2(l008.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
t.test((log2(l009.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
t.test((log2(l011.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
t.test((log2(l018.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
t.test((log2(l021.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
t.test((log2(l031.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
#clearly not trisoic for 3
t.test((log2(l049.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l059.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-
t.test((log2(l061.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
t.test((log2(l066.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l069.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
t.test((log2(l076.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)
t.test((log2(l077.11[,3])), (log2(euRatio.c11)), paired=FALSE, var.equal=TRUE)

abline(h=0,lty=3)

#text(x=19, y=2.2, "*", pos=3, cex=1.5, col="red")

###############################################################################################

#CHROMOSOME 12

#RATIOS FOR CHROMOSOME 12

euRatio.c12 <- c(l002.12[,3],l003.12[,3],l005.12[,3],l006.12[,3],l008.12[,3],l009.12[,3],l069.12[,3])


boxplot((log2(l001.12[,3])),(log2(l002.12[,3])),(log2(l003.12[,3])),(log2(l004.12[,3])),(log2(l005.12[,3])),
        (log2(l006.12[,3])),(log2(l007.12[,3])),(log2(l008.12[,3])),(log2(l009.12[,3])),(log2(l011.12[,3])),
        (log2(l018.12[,3])),(log2(l021.12[,3])),(log2(l031.12[,3])),(log2(l049.12[,3])),(log2(l059.12[,3])),
        (log2(l061.12[,3])),(log2(l066.12[,3])),(log2(l069.12[,3])),(log2(l076.12[,3])),(log2(l077.12[,3])),
        names=c( "1","2","3","4","5","6","7","8","9","11","18","21","31","49","59","61","66","69","76","77"),
        ylab="log2(fold change)", col=c("green", "green","green","cyan","green","green","cyan","cyan",
       "green","cyan","deeppink1","cyan","green","cyan","cyan","cyan","cyan","green","cyan","deeppink1"),
        main="GC Lines Chromsome 12",las=3)
#use non-parametric because these are not normally distributed
#run same test on all lines, even those that are not "supposed" to be aneuploid
t.test((log2(l001.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
t.test((log2(l002.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
t.test((log2(l003.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
t.test((log2(l004.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-
t.test((log2(l005.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
t.test((log2(l006.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
t.test((log2(l007.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
t.test((log2(l008.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
t.test((log2(l009.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
t.test((log2(l011.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
wilcox.test((log2(l018.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
t.test((log2(l021.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
t.test((log2(l031.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
#clearly not trisoi for 3
t.test((log2(l049.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l059.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-
t.test((log2(l061.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
t.test((log2(l066.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l069.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
t.test((log2(l076.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
wilcox.test((log2(l077.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)

abline(h=0,lty=3)

text(x=11, y=3.9, "*", pos=3, cex=1.5, col="red")
text(x=20, y=3.9, "*", pos=3, cex=1.5, col="red")

###############################################################################################

#CHROMOSOME 13

#RATIOS FOR CHROMOSOME 13

euRatio.c13 <- c(l002.13[,3],l003.13[,3],l005.13[,3],l006.13[,3],l008.13[,3],l009.13[,3],l069.13[,3])


boxplot((log2(l001.13[,3])),(log2(l002.13[,3])),(log2(l003.13[,3])),(log2(l004.13[,3])),(log2(l005.13[,3])),
        (log2(l006.13[,3])),(log2(l007.13[,3])),(log2(l008.13[,3])),(log2(l009.13[,3])),(log2(l011.13[,3])),
        (log2(l018.13[,3])),(log2(l021.13[,3])),(log2(l031.13[,3])),(log2(l049.13[,3])),(log2(l059.13[,3])),
        (log2(l061.13[,3])),(log2(l066.13[,3])),(log2(l069.13[,3])),(log2(l076.13[,3])),(log2(l077.13[,3])),
        names=c( "1","2","3","4","5","6","7","8","9","11","18","21","31","49","59","61","66","69","76","77"),
        ylab="log2(fold change)", col=c("green", "green","green","cyan","green","green","cyan","cyan",
           "green","cyan","cyan","cyan","green","cyan","cyan","cyan","cyan","green","cyan","cyan"),
        main="GC Lines Chromsome 13",las=3)
#use non-parametric because these are not normally distributed
#run same test on all lines, even those that are not "supposed" to be aneuploid
t.test((log2(l001.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
t.test((log2(l002.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
t.test((log2(l003.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
t.test((log2(l004.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e
t.test((log2(l005.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
t.test((log2(l006.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
t.test((log2(l007.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
t.test((log2(l008.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
t.test((log2(l009.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
t.test((log2(l011.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
t.test((log2(l018.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
t.test((log2(l021.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
t.test((log2(l031.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
#clearly not trisoi for 3
t.test((log2(l049.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l059.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e
t.test((log2(l061.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
t.test((log2(l066.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l069.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
t.test((log2(l076.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)
t.test((log2(l077.13[,3])), (log2(euRatio.c13)), paired=FALSE, var.equal=TRUE)

abline(h=0,lty=3)

#text(x=11, y=3.9, "*", pos=3, cex=1.5, col="red")
#text(x=20, y=3.9, "*", pos=3, cex=1.5, col="red")

###############################################################################################

#CHROMOSOME 14

#RATIOS FOR CHROMOSOME 14

euRatio.c14 <- c(l002.14[,3],l003.14[,3],l005.14[,3],l006.14[,3],l008.14[,3],l009.14[,3],l069.14[,3])


boxplot((log2(l001.14[,3])),(log2(l002.14[,3])),(log2(l003.14[,3])),(log2(l004.14[,3])),(log2(l005.14[,3])),
        (log2(l006.14[,3])),(log2(l007.14[,3])),(log2(l008.14[,3])),(log2(l009.14[,3])),(log2(l011.14[,3])),
        (log2(l018.14[,3])),(log2(l021.14[,3])),(log2(l031.14[,3])),(log2(l049.14[,3])),(log2(l059.14[,3])),
        (log2(l061.14[,3])),(log2(l066.14[,3])),(log2(l069.14[,3])),(log2(l076.14[,3])),(log2(l077.14[,3])),
        names=c( "1","2","3","4","5","6","7","8","9","11","18","21","31","49","59","61","66","69","76","77"),
        ylab="log2(fold change)", col=c("green", "green","green","cyan","green","green","cyan","cyan",
               "green","cyan","cyan","cyan","green","cyan","cyan","cyan","cyan","green","deeppink1","cyan"),
        main="GC Lines Chromsome 14",las=3)
#use non-parametric because these are not normally distributed
#run same test on all lines, even those that are not "supposed" to be aneuploid
t.test((log2(l001.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
t.test((log2(l002.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
t.test((log2(l003.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
t.test((log2(l004.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e
t.test((log2(l005.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
t.test((log2(l006.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
t.test((log2(l007.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
t.test((log2(l008.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
t.test((log2(l009.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
t.test((log2(l011.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
t.test((log2(l018.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
t.test((log2(l021.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
t.test((log2(l031.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
#clearly not trisoi for 3
t.test((log2(l049.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l059.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e
t.test((log2(l061.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
t.test((log2(l066.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l069.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
wilcox.test((log2(l076.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
t.test((log2(l077.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)

abline(h=0,lty=3)

text(x=19, y=3, "*", pos=3, cex=1.5, col="red")
#text(x=20, y=3.9, "*", pos=3, cex=1.5, col="red")

###############################################################################################

#CHROMOSOME 15

#RATIOS FOR CHROMOSOME 15

euRatio.c15 <- c(l002.15[,3],l003.15[,3],l005.15[,3],l006.15[,3],l008.15[,3],l009.15[,3],l069.15[,3])


boxplot((log2(l001.15[,3])),(log2(l002.15[,3])),(log2(l003.15[,3])),(log2(l004.15[,3])),(log2(l005.15[,3])),
        (log2(l006.15[,3])),(log2(l007.15[,3])),(log2(l008.15[,3])),(log2(l009.15[,3])),(log2(l011.15[,3])),
        (log2(l018.15[,3])),(log2(l021.15[,3])),(log2(l031.15[,3])),(log2(l049.15[,3])),(log2(l059.15[,3])),
        (log2(l061.15[,3])),(log2(l066.15[,3])),(log2(l069.15[,3])),(log2(l076.15[,3])),(log2(l077.15[,3])),
        names=c( "1","2","3","4","5","6","7","8","9","11","18","21","31","49","59","61","66","69","76","77"),
        ylab="log2(fold change)", col=c("green", "green","green","cyan","green","green","cyan","cyan",
           "green","deeppink1","cyan","cyan","green","cyan","cyan","cyan","cyan","green","cyan","cyan"),
        main="GC Lines Chromsome 15",las=3)
#use non-parametric because these are not normally distributed
#run same test on all lines, even those that are not "supposed" to be aneuploid
t.test((log2(l001.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
t.test((log2(l002.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
t.test((log2(l003.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
t.test((log2(l004.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
#p-value < 2.25
t.test((log2(l005.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
t.test((log2(l006.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
t.test((log2(l007.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
t.test((log2(l008.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
t.test((log2(l009.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
wilcox.test((log2(l011.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
t.test((log2(l018.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
t.test((log2(l021.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
t.test((log2(l031.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
#clearly not trisoi5for 3
t.test((log2(l049.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l059.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e
t.test((log2(l061.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
t.test((log2(l066.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l069.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
t.test((log2(l076.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
t.test((log2(l077.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)

abline(h=0,lty=3)

text(x=10, y=2.5, "*", pos=3, cex=1.5, col="red")
#text(x=20, y=3.9, "*", pos=3, cex=1.5, col="red")

###############################################################################################

#CHROMOSOME 16

#RATIOS FOR CHROMOSOME 16

euRatio.c16 <- c(l002.16[,3],l003.16[,3],l005.16[,3],l006.16[,3],l008.16[,3],l009.16[,3],l069.16[,3])


boxplot((log2(l001.16[,3])),(log2(l002.16[,3])),(log2(l003.16[,3])),(log2(l004.16[,3])),(log2(l005.16[,3])),
        (log2(l006.16[,3])),(log2(l007.16[,3])),(log2(l008.16[,3])),(log2(l009.16[,3])),(log2(l011.16[,3])),
        (log2(l018.16[,3])),(log2(l021.16[,3])),(log2(l031.16[,3])),(log2(l049.16[,3])),(log2(l059.16[,3])),
        (log2(l061.16[,3])),(log2(l066.16[,3])),(log2(l069.16[,3])),(log2(l076.16[,3])),(log2(l077.16[,3])),
        names=c( "1","2","3","4","5","6","7","8","9","11","18","21","31","49","59","61","66","69","76","77"),
        ylab="log2(fold change)", col=c("green", "green","green","cyan","green","green","cyan","deeppink1",
              "green","cyan","cyan","cyan","green","cyan","cyan","cyan","cyan","green","cyan","cyan"),
        main="GC Lines Chromsome 16",las=3)
#use non-parametric because these are not normally distributed
#run same test on all lines, even those that are not "supposed" to be aneuploid
t.test((log2(l001.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
t.test((log2(l002.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
t.test((log2(l003.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
t.test((log2(l004.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
#p-value < 2.26
t.test((log2(l005.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
t.test((log2(l006.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
t.test((log2(l007.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
wilcox.test((log2(l008.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
t.test((log2(l009.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
t.test((log2(l011.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
t.test((log2(l018.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
t.test((log2(l021.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
t.test((log2(l031.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
#clearly not trisoi5for 3
t.test((log2(l049.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l059.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e
t.test((log2(l061.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-16
t.test((log2(l066.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
#p-value < 2.2e-1
t.test((log2(l069.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
t.test((log2(l076.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
t.test((log2(l077.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)

abline(h=0,lty=3)

text(x=8, y=3.9, "*", pos=3, cex=1.5, col="red")
#text(x=20, y=3.9, "*", pos=3, cex=1.5, col="red")
###############################################################################################
##boxplot and t-test for chromosome 5 trisomic
#boxplot((log2(anRatio.c5)), (log2(euRatio.c5)),names=c("Trisomic", "Euploid"),ylab="log2(fold change)", col=c("limegreen", "purple"))
#t.test((log2(anRatio.c5)), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
#t.test((log2(anRatio.c5.no50)), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
#t.test((log2(anRatio.c5.no4)), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
###p value is 0.253 whaaat

#boxplot((log2(euRatio.c5)),(log2(anRatio.c5)),names=c("Disomic","Trisomic"),ylab="log2(fold change)", col=c("turquoise", "hotpink"), main="Chr 5 Relative Gene Expression Levels",par(mfrow=c(1,1))
#abline(h=0,lty=3)
#t.test((log2(anRatio.c5)), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
##p value for monosomic: 0.253
#wilcox.test((log2(anRatio.c5)), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
##p value wilcox tets: 0.644
#p <- 0.0002234
#mylabel = bquote(italic(p) == .(format(p, digits = 9)))
##text(x=3.2,y = 3.7, labels = mylabel,cex=.8)
#text(x=2, y=3, "*", pos=3, cex=1.5, col="red")
##going to need these for all of the chromosomes that have a CNV, but can come back to this later
###################################################################################################
#CHROMOSOME 3


#RATIOS FOR CHROMOSOME 3
euRatio.c3 <- c(l002.3[,3],l003.3[,3],l005.3[,3],l006.3[,3],l009.3[,3])



l115.3 <- getChrmRatio("115",3)
boxplot( (log2(euRatio.c3)),(log2(l115.3[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine", "mediumvioletred"),main="Line 115 Chr 3",las=3)
t.test((log2(l115.3[,3])), (log2(euRatio.c3)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <- 0.00544
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.4,y = 2.5, labels = mylabel,cex=.8)

#####################################################################################################
###going to look at individual lines 
par(mfrow=c(1,1))
##Line 117
boxplot((log2(euRatio.c5)), (log2(l117.5[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine", "mediumvioletred"),main="Line 117",las=3)
t.test((log2(l117.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p2 <- 1.026e-13
mylabel = bquote(italic(p) == .(format(p2, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=3, "*", pos=3, cex=1.5, col="red")

##significantly different
boxplot((log2(euRatio.c5)),(log2(l049.5[,3])),names=c( "Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine", "mediumvioletred"),main="Line 49",las=3)
t.test((log2(l049.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p3 <- 9.394e-14
mylabel = bquote(italic(p) == .(format(p3, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=2.5, "*", pos=3, cex=1.5, col="red")

##significantly different
boxplot( (log2(euRatio.c5)),(log2(l004.5[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine", "mediumvioletred"),main="Line 4",las=3)
t.test((log2(l004.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p4 <- 1.414e-14
mylabel = bquote(italic(p) == .(format(p4, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=2.5, "*", pos=3, cex=1.5, col="red")

##significantly different
boxplot( (log2(euRatio.c5)),(log2(l050.5[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine", "mediumvioletred"),main="Line 50 Chr 5",las=3)
t.test((log2(l050.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <- 0.4559
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

l123.5 <- getChrmRatio("123",5)
boxplot( (log2(euRatio.c5)),(log2(l123.5[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine", "mediumvioletred"),main="Line 123 Chr 5",las=3)
t.test((log2(l123.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <- 4e-06
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.4,y = 2.5, labels = mylabel,cex=.8)

##significantly different 

##############################################################################################

###CHROMOSOME 7 

l115.7 <- getChrmRatio("115",7)
l031.7 <- getChrmRatio("31",7)
l059.7 <- getChrmRatio("59",7)
l061.7 <- getChrmRatio("61",7)
l066.7 <- getChrmRatio("66",7)

#RATIOS FOR CHROMOSOME 7
anRatio.c7 <- c(l115.7[,3],l031.7[,3],l059.7[,3],l061.7[,3],l066.7[,3])
euRatio.c7 <- c(l002.7[,3],l003.7[,3],l005.7[,3],l006.7[,3],l009.7[,3])

##boxplot and t-test for chromosome 7 trisomic
#boxplot((log2(anRatio.c7)), (log2(euRatio.c7)),names=c("Trisomic", "Euploid"),ylab="log2(fold change)", col=c("blue", "red"),main="Chr 7 Gene Expression")
#t.test((log2(anRatio.c7)), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
#abline(h=0,lty=3)
#p5 <- 2.2e-16
#mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
#text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

#Line 66
boxplot( (log2(euRatio.c7)),(log2(l066.7[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan", "violetred1"),main="Chr 7 Gene Expression Line 66", las=3)
t.test((log2(l066.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <- 1.337e-06
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=3, "*", pos=3, cex=1.5, col="red")

#Line 31
boxplot( (log2(euRatio.c7)),(log2(l031.7[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan", "violetred1"),main="Chr 7 Gene Expression Line 31", las=3)
t.test((log2(l031.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
p5 <- 0.7816
abline(h=0,lty=3)
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

#line 59
boxplot( (log2(euRatio.c7)),(log2(l059.7[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan", "violetred1"),main="Chr 7 Gene Expression Line 59", las=3)
t.test((log2(l059.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
p5 <- "< 2.2e-16"
abline(h=0,lty=3)
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=5, "*", pos=3, cex=1.5, col="red")
#line 61
boxplot( (log2(euRatio.c7)),(log2(l061.7[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan", "violetred1"),main="Chr 7 Gene Expression Line 61", las=3)
t.test((log2(l061.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
p5 <- "< 2.2e-16"
abline(h=0,lty=3)
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=2.5, "*", pos=3, cex=1.5, col="red")
#line 115
boxplot( (log2(euRatio.c7)),(log2(l115.7[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan", "violetred1"),main="Chr 7 Gene Expression Line 115", las=3)
t.test((log2(l115.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
p5 <- 0.1624
abline(h=0,lty=3)
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

##############################################################################################

###CHROMOSOME 8
l152.8 <- getChrmRatio("152",8)
l029.8 <- getChrmRatio("29",8)
##ratios for chromosome 8 
euRatio.c8 <- c(l002.8[,3],l003.8[,3],l005.8[,3],l006.8[,3],l009.8[,3])

#Line 152
boxplot((log2(euRatio.c8)),(log2(l152.8[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("mediumorchid3", "deeppink1"),main="Chr 8 Gene Expression Line 152", las=3)
t.test((log2(l152.8[,3])), (log2(euRatio.c8)), paired=FALSE, var.equal=TRUE)
p5 <- 3.599e-16
abline(h=0,lty=3)
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=3, "*", pos=3, cex=1.5, col="red")


#############################################################################################

###CHROMOSOME 9

l029.9 <- getChrmRatio("29",9)
l076.9 <- getChrmRatio("76",9)

#RATIOS FOR CHROMOSOME 9
anRatio.c9 <- c(l029.9[,3],l076.9[,3])
euRatio.c9 <- c(l002.9[,3],l003.9[,3],l005.9[,3],l006.9[,3],l009.9[,3])

#Line 29
boxplot((log2(euRatio.c9)),(log2(l029.9[,3])),names=c("Disomic", "Monosomic"),ylab="log2(fold change)", col=c("cyan3", "mediumvioletred"),main="Chr 9 Gene Expression Line 29", las=3)
t.test((log2(l029.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
p5 <- "< 2.2e-16"
abline(h=0,lty=3)
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=3, "*", pos=3, cex=1.5, col="red")
#Line76
boxplot((log2(euRatio.c9)),(log2(l076.9[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan3", "mediumvioletred"),main="Chr 9 Gene Expression Line 76", las=3)
t.test((log2(l076.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
p5 <- "< 2.2e-16"
abline(h=0,lty=3)
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=2.7, "*", pos=3, cex=1.5, col="red")
##############################################################################################

###Chromosome 10 

l066.10 <- getChrmRatio("66",10)
l076.10 <- getChrmRatio("76",10)


#RATIOS FOR CHROMOSOME 10
anRatio.c10 <- c(l066.10[,3])
euRatio.c10 <- c(l002.10[,3],l003.10[,3],l005.10[,3],l006.10[,3],l009.10[,3])

##boxplot and t-test for chromosome 10 line 76 trisomic for one arm of 10
boxplot( (log2(euRatio.c10)),(log2(l066.10[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan2", "deeppink1"),main="Chr 10 Gene Expression: Line 66")
t.test((log2(l066.10[,3])), (log2(euRatio.c10)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-0.4013
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

##############################################################################################

###Chromosome 12

l001.12 <- getChrmRatio("1",12)
l018.12 <- getChrmRatio("18",12)
l077.12 <- getChrmRatio("77",12)
l123.12 <- getChrmRatio("123",12)

#RATIOS FOR CHROMOSOME 10
anRatio.c12 <- c(l001.12[,3],l018.12[,3],l077.12[,3],l123.12[,3])
euRatio.c12 <- c(l002.12[,3],l003.12[,3],l005.12[,3],l006.12[,3],l009.12[,3])

#line1 
boxplot( (log2(euRatio.c12)),(log2(l001.12[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine2", "mediumvioletred"),main="Chr 12 Gene Expression: Line 1",las=3)
t.test((log2(l001.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <- 0.1625
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

#line18 
boxplot( (log2(euRatio.c12)),(log2(l018.12[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine2", "mediumvioletred"),main="Chr 12 Gene Expression: Line 18",las=3)
t.test((log2(l018.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-"< 2.2e-16"
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=7, "*", pos=3, cex=1.5, col="red")
#line77 
boxplot( (log2(euRatio.c12)),(log2(l077.12[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine2", "mediumvioletred"),main="Chr 12 Gene Expression: Line 77",las=3)
t.test((log2(l077.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-"< 2.2e-16"
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=3, "*", pos=3, cex=1.5, col="red")
#line123 
boxplot( (log2(euRatio.c12)),(log2(l123.12[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine2", "mediumvioletred"),main="Chr 12 Gene Expression: Line 123",las=3)
t.test((log2(l123.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-0.2541
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

##############################################################################################

####Chromosome 14

l076.14 <- getChrmRatio("76",14)


#RATIOS FOR CHROMOSOME 14
anRatio.c14 <- c(l076.14[,3])
euRatio.c14 <- c(l002.14[,3],l003.14[,3],l005.14[,3],l006.14[,3],l009.14[,3])

boxplot( (log2(euRatio.c14)),(log2(l076.14[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine3", "deeppink1"),main="Chr 14 Gene Expression: Line 76",las=3)
t.test((log2(l076.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-"< 2.2e-16"
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=2.5, "*", pos=3, cex=1.5, col="red")

##############################################################################################

##Chromosome 15
l011.15 <- getChrmRatio("11",15)


#RATIOS FOR CHROMOSOME 10
anRatio.c15 <- c(l011.15[,3])
euRatio.c15 <- c(l002.15[,3],l003.15[,3],l005.15[,3],l006.15[,3],l009.15[,3])

boxplot( (log2(euRatio.c15)),(log2(l011.15[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine1", "violetred1"),main="Chr 15 Gene Expression: Line 11",las=3)
t.test((log2(l011.15[,3])), (log2(euRatio.c15)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-"< 2.2e-16"
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=5.3, "*", pos=3, cex=1.5, col="red")
##############################################################################################

###Chromosome 16 
l008.16 <- getChrmRatio("8",16)
l031.16 <- getChrmRatio("31",16)
l069.16 <- getChrmRatio("69",16)
l112.16 <- getChrmRatio("112",16)
l141.16 <- getChrmRatio("141",16)

#RATIOS FOR CHROMOSOME 16
anRatio.c16 <- c(l008.16[,3],l031.16[,3],l069.16[,3],l112.16[,3],l141.16[,3])
euRatio.c16 <- c(l002.16[,3],l003.16[,3],l005.16[,3],l006.16[,3],l009.16[,3])
#line 8
boxplot( (log2(euRatio.c16)),(log2(l008.16[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan3", "mediumorchid3"),main="Chr 16 Gene Expression: Line 8",las=3)
t.test((log2(l008.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-"< 2.2e-16"
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=3.3, "*", pos=3, cex=1.5, col="red")
#line31
boxplot( (log2(euRatio.c16)),(log2(l031.16[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan3", "mediumorchid3"),main="Chr 16 Gene Expression: Line 31",las=3)
t.test((log2(l031.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-0.1757
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

#line69
boxplot( (log2(euRatio.c16)),(log2(l069.16[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan3", "mediumorchid3"),main="Chr 16 Gene Expression: Line 69",las=3)
t.test((log2(l069.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-0.0951
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

#line112
boxplot( (log2(euRatio.c16)),(log2(l112.16[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan3", "mediumorchid3"),main="Chr 16 Gene Expression: Line 112",las=3)
t.test((log2(l112.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-"< 2.2e-16"
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=5, "*", pos=3, cex=1.5, col="red")
#line141
boxplot( (log2(euRatio.c16)),(log2(l141.16[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan3", "mediumorchid3"),main="Chr 16 Gene Expression: Line 141",las=3)
t.test((log2(l141.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-"< 2.2e-16"
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=3, "*", pos=3, cex=1.5, col="red")
##############################################################################################

#make a histogram of log ratio of all genes across 
#euploid lines
##for chromosome one 
##can't do this for euploid ratios
#has to be individual euploid lines 

###CHROMOSOME I 
par(mfrow=c(5,2))
par(mar = rep(2, 4))

l002.1 <- getChrmRatio("2",1)
l003.1 <- getChrmRatio("3",1)
l005.1 <- getChrmRatio("5",1)
l006.1 <- getChrmRatio("6",1)
l009.1 <- getChrmRatio("9",1)
l152.7 <- getChrmRatio("152",7)
l152.1 <- getChrmRatio("152",1)


hist(log2(l002.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 1 Line 2")
skewness(log2(l002.1[,3]))
#distribution is left skewed
##skewness = -3.589928
mean(log2(l002.1[,3]))
#mean is -0.1290089
sd(log2(l002.1[,3]))
#sd is 1.136425
lines(density(log2(l002.1[,3])),col="green")
mc11 <- mean(log2(l002.1[,3]))
mc11
abline(v = mc11, col = "blue")
#make a histogram of log ratio of all genes across 
#euploid lines
s2 <- -3.589928
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .75, labels = mylabel,cex=.8)
sd2 <- sd(log2(l002.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .8, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .7, labels = mylabel3,cex=.8)

hist(log2(l003.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 1 Line 3")
#distribution is left skewed
##skewness = -3.589928
mean(log2(l003.1[,3]))
#mean is -0.1290089
sd(log2(l003.1[,3]))
#sd is 1.136425
lines(density(log2(l003.1[,3])),col="green")
mc11 <- mean(log2(l003.1[,3]))
mc11
abline(v = mc11, col = "blue")
#make a histogram of log ratio of all genes across 
#euploid lines
s2 <- skewness(log2(l003.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .85, labels = mylabel,cex=.8)
sd2 <- sd(log2(l003.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .8, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .75, labels = mylabel3,cex=.8)

hist(log2(l005.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 1 Line 5")
#distribution is left skewed
##skewness = -3.589928
mean(log2(l005.1[,3]))
#mean is -0.1290089
sd(log2(l005.1[,3]))
#sd is 1.136425
lines(density(log2(l005.1[,3])),col="green")
mc11 <- mean(log2(l005.1[,3]))
mc11
abline(v = mc11, col = "blue")
#make a histogram of log ratio of all genes across 
#euploid lines
s2 <- skewness(log2(l005.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l005.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .5, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .45, labels = mylabel3,cex=.8)

hist(log2(l006.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 1 Line 6")
#distribution is left skewed
##skewness = -3.589928
mean(log2(l006.1[,3]))
#mean is -0.1290089
sd(log2(l006.1[,3]))
#sd is 1.136425
lines(density(log2(l006.1[,3])),col="green")
mc11 <- mean(log2(l006.1[,3]))
mc11
abline(v = mc11, col = "blue")
#make a histogram of log ratio of all genes across 
#euploid lines
s2 <- skewness(log2(l006.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .55, labels = mylabel,cex=.8)
sd2 <- sd(log2(l006.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .5, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .45, labels = mylabel3,cex=.8)

hist(log2(l009.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 1 Line 9")
#distribution is left skewed
##skewness = -3.589928
mean(log2(l009.1[,3]))
#mean is -0.1290089
sd(log2(l009.1[,3]))
#sd is 1.136425
lines(density(log2(l009.1[,3])),col="green")
mc11 <- mean(log2(l009.1[,3]))
mc11
abline(v = mc11, col = "blue")
#make a histogram of log ratio of all genes across 
#euploid lines
s2 <- skewness(log2(l009.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .75, labels = mylabel,cex=.8)
sd2 <- sd(log2(l009.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .65, labels = mylabel3,cex=.8)


hist(log2(l007.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Trisomic Chr 1 Line 7")
#distribution is left skewed
##skewness = -3.589928
mean(log2(l007.1[,3]))
#mean is -0.1290089
sd(log2(l007.1[,3]))
#sd is 1.136425
lines(density(log2(l007.1[,3])),col="green")
mc11 <- mean(log2(l007.1[,3]))
mc11
abline(v = mc11, col = "blue")
#make a histogram of log ratio of all genes across 
#euploid lines
s2 <- skewness(log2(l007.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .75, labels = mylabel,cex=.8)
sd2 <- sd(log2(l007.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .65, labels = mylabel3,cex=.8)

hist(log2(l021.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Trisomic Chr 1 Line 21")
#distribution is left skewed
##skewness = -3.589928
mean(log2(l021.1[,3]))
#mean is -0.1290089
sd(log2(l021.1[,3]))
#sd is 1.136425
lines(density(log2(l021.1[,3])),col="green")
mc11 <- mean(log2(l021.1[,3]))
mc11
abline(v = mc11, col = "blue")
#make a histogram of log ratio of all genes across 
#euploid lines
s2 <- skewness(log2(l021.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .75, labels = mylabel,cex=.8)
sd2 <- sd(log2(l021.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .7, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .65, labels = mylabel3,cex=.8)

hist(log2(l011.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Monosomic Chr 1 Line 11")
#distribution is left skewed
##skewness = -3.589928
mean(log2(l011.1[,3]))
#mean is -0.1290089
sd(log2(l011.1[,3]))
#sd is 1.136425
lines(density(log2(l011.1[,3])),col="green")
mc11 <- mean(log2(l011.1[,3]))
mc11
abline(v = mc11, col = "blue")
#make a histogram of log ratio of all genes across 
#euploid lines
s2 <- skewness(log2(l011.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .9, labels = mylabel,cex=.8)
sd2 <- sd(log2(l011.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .8, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .7, labels = mylabel3,cex=.8)

hist(log2(l018.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Trisomic Chr 1 Line 18")
#distribution is left skewed
##skewness = -3.589928
mean(log2(l018.1[,3]))
#mean is -0.1290089
sd(log2(l018.1[,3]))
#sd is 1.136425
lines(density(log2(l018.1[,3])),col="green")
mc11 <- mean(log2(l018.1[,3]))
mc11
abline(v = mc11, col = "blue")
#make a histogram of log ratio of all genes across 
#euploid lines
s2 <- skewness(log2(l018.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .9, labels = mylabel,cex=.8)
sd2 <- sd(log2(l018.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .8, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .7, labels = mylabel3,cex=.8)

hist(log2(l152.1[,3]),freq=FALSE,breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Trisomic Chr 1 Line 152")
#distribution is left skewed
##skewness = -3.589928
mean(log2(l152.1[,3]))
#mean is -0.1290089
sd(log2(l152.1[,3]))
#sd is 1.136425
lines(density(log2(l152.1[,3])),col="green")
mc11 <- mean(log2(l152.1[,3]))
mc11
abline(v = mc11, col = "blue")
#make a histogram of log ratio of all genes across 
#euploid lines
s2 <- skewness(log2(l152.1[,3]))
mylabel = bquote(italic("Skewness") == .(format(s2, digits = 3)))
text(x=2.2,y = .9, labels = mylabel,cex=.8)
sd2 <- sd(log2(l152.1[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd2, digits = 3)))
text(x=2.2,y = .8, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = .7, labels = mylabel3,cex=.8)


###CHROMOSOME 5 

##to make a histogram that has the density curve included and the mean plotted
##line 11, monosomic for chromosome 1 
hist(log2(l002.5[,3]),freq=FALSE,breaks=40, main="Euploid Line 2 Chr 5", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l002.5[,3])),col="green")
mc11 <- mean(log2(l002.5[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l002.5[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l002.5[,3]))
#mean is -0.8662616
sd(log2(l002.5[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l002.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l002.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l003.5[,3]),freq=FALSE,breaks=40, main="Euploid Line 3 Chr 5", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l003.5[,3])),col="green")
mc11 <- mean(log2(l003.5[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l003.5[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l003.5[,3]))
#mean is -0.8662616
sd(log2(l003.5[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l003.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <- sd(log2(l003.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l005.5[,3]),freq=FALSE,breaks=40, main="Euploid Line 5 Chr 5", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l005.5[,3])),col="green")
mc11 <- mean(log2(l005.5[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l005.5[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l005.5[,3]))
#mean is -0.8662616
sd(log2(l005.5[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l005.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l005.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l006.5[,3]),freq=FALSE,breaks=40, main="Euploid Line 6 Chr 5", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l006.5[,3])),col="green")
mc11 <- mean(log2(l006.5[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l006.5[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l006.5[,3]))
#mean is -0.8662616
sd(log2(l006.5[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l006.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l006.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l009.5[,3]),freq=FALSE,breaks=40, main="Euploid Line 9 Chr 5", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l009.5[,3])),col="green")
mc11 <- mean(log2(l009.5[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l009.5[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l009.5[,3]))
#mean is -0.8662616
sd(log2(l009.5[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l009.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l009.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l004.5[,3]),freq=FALSE,breaks=40, main="Trisomic Line 4 Chr 5", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l004.5[,3])),col="green")
mc11 <- mean(log2(l004.5[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l004.5[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l004.5[,3]))
#mean is -0.8662616
sd(log2(l004.5[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l004.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l004.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l049.5[,3]),freq=FALSE,breaks=40, main="Trisomic Line 49 Chr 5", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l049.5[,3])),col="green")
mc11 <- mean(log2(l049.5[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l049.5[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l049.5[,3]))
#mean is -0.8662616
sd(log2(l049.5[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l049.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l049.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l050.5[,3]),freq=FALSE,breaks=40, main="Trisomic Line 50 Chr 5", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l050.5[,3])),col="green")
mc11 <- mean(log2(l050.5[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l050.5[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l050.5[,3]))
#mean is -0.8662616
sd(log2(l050.5[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l050.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l050.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l117.5[,3]),freq=FALSE,breaks=40, main="Trisomic Line 117 Chr 5", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l117.5[,3])),col="green")
mc11 <- mean(log2(l117.5[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l117.5[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l117.5[,3]))
#mean is -0.8662616
sd(log2(l117.5[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l117.5[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l117.5[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

###CHROMOSOME 7

hist(log2(l002.7[,3]),freq=FALSE,breaks=40, main="Euploid Line 2 Chr 7", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l002.7[,3])),col="green")
mc11 <- mean(log2(l002.7[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l002.7[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l002.7[,3]))
#mean is -0.8662616
sd(log2(l002.7[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l002.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l002.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l003.7[,3]),freq=FALSE,breaks=40, main="Euploid Line 3 Chr 7", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l003.7[,3])),col="green")
mc11 <- mean(log2(l003.7[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l003.7[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l003.7[,3]))
#mean is -0.8662616
sd(log2(l003.7[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l003.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l003.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l005.7[,3]),freq=FALSE,breaks=40, main="Euploid Line 5 Chr 7", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l005.7[,3])),col="green")
mc11 <- mean(log2(l005.7[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l005.7[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l005.7[,3]))
#mean is -0.8662616
sd(log2(l005.7[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l005.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l005.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l006.7[,3]),freq=FALSE,breaks=40, main="Euploid Line 6 Chr 7", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l006.7[,3])),col="green")
mc11 <- mean(log2(l006.7[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l006.7[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l006.7[,3]))
#mean is -0.8662616
sd(log2(l006.7[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l006.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l006.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l009.7[,3]),freq=FALSE,breaks=40, main="Euploid Line 9 Chr 7", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l009.7[,3])),col="green")
mc11 <- mean(log2(l009.7[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l009.7[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l009.7[,3]))
#mean is -0.8662616
sd(log2(l009.7[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l009.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l009.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l031.7[,3]),freq=FALSE,breaks=40, main="Trisomic Line 31 Chr 7", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.5))
lines(density(log2(l031.7[,3])),col="green")
mc11 <- mean(log2(l031.7[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l031.7[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l031.7[,3]))
#mean is -0.8662616
sd(log2(l031.7[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l031.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l031.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l059.7[,3]),freq=FALSE,breaks=40, main="Trisomic Line 59 Chr 7", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.5))
lines(density(log2(l059.7[,3])),col="green")
mc11 <- mean(log2(l059.7[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l059.7[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l059.7[,3]))
#mean is -0.8662616
sd(log2(l059.7[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l059.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l059.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l066.7[,3]),freq=FALSE,breaks=40, main="Trisomic Line 66 Chr 7", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.5))
lines(density(log2(l066.7[,3])),col="green")
mc11 <- mean(log2(l066.7[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l066.7[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l066.7[,3]))
#mean is -0.8662616
sd(log2(l066.7[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l066.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l066.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l061.7[,3]),freq=FALSE,breaks=40, main="Trisomic Line 61 Chr 7", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.5))
lines(density(log2(l061.7[,3])),col="green")
mc11 <- mean(log2(l061.7[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l061.7[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l061.7[,3]))
#mean is -0.8662616
sd(log2(l061.7[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l061.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l061.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l115.7[,3]),freq=FALSE,breaks=40, main="Trisomic Line 115 Chr 7", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l115.7[,3])),col="green")
mc11 <- mean(log2(l115.7[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l115.7[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l115.7[,3]))
#mean is -0.8662616
sd(log2(l115.7[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l115.7[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l115.7[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

###CHROMOSOME 8 
par(mfrow=c(3,2))

hist(log2(l002.8[,3]),freq=FALSE,breaks=40, main="Euploid Line 2 Chr 8", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l002.8[,3])),col="green")
mc11 <- mean(log2(l002.8[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l002.8[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l002.8[,3]))
#mean is -0.8662616
sd(log2(l002.8[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l002.8[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l002.8[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l003.8[,3]),freq=FALSE,breaks=40, main="Euploid Line 3 Chr 8", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l003.8[,3])),col="green")
mc11 <- mean(log2(l003.8[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l003.8[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l003.8[,3]))
#mean is -0.8662616
sd(log2(l003.8[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l003.8[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l003.8[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l005.8[,3]),freq=FALSE,breaks=40, main="Euploid Line 5 Chr 8", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l005.8[,3])),col="green")
mc11 <- mean(log2(l005.8[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l005.8[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l005.8[,3]))
#mean is -0.8662616
sd(log2(l005.8[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l005.8[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l005.8[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l006.8[,3]),freq=FALSE,breaks=40, main="Euploid Line 6 Chr 8", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l006.8[,3])),col="green")
mc11 <- mean(log2(l006.8[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l006.8[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l006.8[,3]))
#mean is -0.8662616
sd(log2(l006.8[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l006.8[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l006.8[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l009.8[,3]),freq=FALSE,breaks=40, main="Euploid Line 9 Chr 8", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l009.8[,3])),col="green")
mc11 <- mean(log2(l009.8[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l009.8[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l009.8[,3]))
#mean is -0.8662616
sd(log2(l009.8[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l009.8[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l009.8[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l152.8[,3]),freq=FALSE,breaks=40, main="Trisomic Line 152 Chr 8", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l152.8[,3])),col="green")
mc11 <- mean(log2(l152.8[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l152.8[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l152.8[,3]))
#mean is -0.8662616
sd(log2(l152.8[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l152.8[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l152.8[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

###CHROMOSOME 9 
par(mfrow=c(4,2))
hist(log2(l002.9[,3]),freq=FALSE,breaks=40, main="Euploid Line 2 Chr 9", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l002.9[,3])),col="green")
mc11 <- mean(log2(l002.9[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l002.9[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l002.9[,3]))
#mean is -0.8662616
sd(log2(l002.9[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l002.9[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l002.9[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l003.9[,3]),freq=FALSE,breaks=40, main="Euploid Line 3 Chr 9", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l003.9[,3])),col="green")
mc11 <- mean(log2(l003.9[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l003.9[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l003.9[,3]))
#mean is -0.8662616
sd(log2(l003.9[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l003.9[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l003.9[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l005.9[,3]),freq=FALSE,breaks=40, main="Euploid Line 5 Chr 9", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l005.9[,3])),col="green")
mc11 <- mean(log2(l005.9[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l005.9[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l005.9[,3]))
#mean is -0.8662616
sd(log2(l005.9[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l005.9[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l005.9[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l006.9[,3]),freq=FALSE,breaks=40, main="Euploid Line 6 Chr 9", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l006.9[,3])),col="green")
mc11 <- mean(log2(l006.9[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l006.9[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l006.9[,3]))
#mean is -0.8662616
sd(log2(l006.9[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l006.9[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l006.9[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l009.9[,3]),freq=FALSE,breaks=40, main="Euploid Line 9 Chr 9", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l009.9[,3])),col="green")
mc11 <- mean(log2(l009.9[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l009.9[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l009.9[,3]))
#mean is -0.8662616
sd(log2(l009.9[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l009.9[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l009.9[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l076.9[,3]),freq=FALSE,breaks=40, main="Trisomic Line 76 Chr 9", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l076.9[,3])),col="green")
mc11 <- mean(log2(l076.9[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l076.9[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l076.9[,3]))
#mean is -0.8662616
sd(log2(l076.9[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l076.9[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l076.9[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l029.9[,3]),freq=FALSE,breaks=40, main="Monosomic Line 29 Chr 9", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l029.9[,3])),col="green")
mc11 <- mean(log2(l029.9[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l029.9[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l029.9[,3]))
#mean is -0.8662616
sd(log2(l029.9[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l029.9[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l029.9[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)


###CHROMOSOME 12 
par(mfrow=c(5,2))
hist(log2(l002.12[,3]),freq=FALSE,breaks=40, main="Euploid Line 2 Chr 12", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l002.12[,3])),col="green")
mc11 <- mean(log2(l002.12[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l002.12[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l002.12[,3]))
#mean is -0.8662616
sd(log2(l002.12[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l002.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l002.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l003.12[,3]),freq=FALSE,breaks=40, main="Euploid Line 3 Chr 12", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l003.12[,3])),col="green")
mc11 <- mean(log2(l003.12[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l003.12[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l003.12[,3]))
#mean is -0.8662616
sd(log2(l003.12[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l003.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l003.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l005.12[,3]),freq=FALSE,breaks=40, main="Euploid Line 5 Chr 12", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l005.12[,3])),col="green")
mc11 <- mean(log2(l005.12[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l005.12[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l005.12[,3]))
#mean is -0.8662616
sd(log2(l005.12[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l005.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l005.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l006.12[,3]),freq=FALSE,breaks=40, main="Euploid Line 6 Chr 12", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l006.12[,3])),col="green")
mc11 <- mean(log2(l006.12[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l006.12[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l006.12[,3]))
#mean is -0.8662616
sd(log2(l006.12[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l006.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l006.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l009.12[,3]),freq=FALSE,breaks=40, main="Euploid Line 9 Chr 12", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l009.12[,3])),col="green")
mc11 <- mean(log2(l009.12[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l009.12[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l009.12[,3]))
#mean is -0.8662616
sd(log2(l009.12[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l009.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l009.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l001.12[,3]),freq=FALSE,breaks=40, main="Trisomic Line 1 Chr 12", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l001.12[,3])),col="green")
mc11 <- mean(log2(l001.12[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l001.12[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l001.12[,3]))
#mean is -0.8662616
sd(log2(l001.12[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l001.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l001.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l018.12[,3]),freq=FALSE,breaks=40, main="Trisomic Line 18 Chr 12", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l018.12[,3])),col="green")
mc11 <- mean(log2(l018.12[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l018.12[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l018.12[,3]))
#mean is -0.8662616
sd(log2(l018.12[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l018.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l018.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)


hist(log2(l077.12[,3]),freq=FALSE,breaks=40, main="Trisomic Line 77 Chr 12", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l077.12[,3])),col="green")
mc11 <- mean(log2(l077.12[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l077.12[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l077.12[,3]))
#mean is -0.8662616
sd(log2(l077.12[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l077.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l077.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l123.12[,3]),freq=FALSE,breaks=40, main="Trisomic Line 123 Chr 12", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l123.12[,3])),col="green")
mc11 <- mean(log2(l123.12[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l123.12[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l123.12[,3]))
#mean is -0.8662616
sd(log2(l123.12[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l123.12[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l123.12[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)



###CHROMOSOME 14 
par(mfrow=c(3,2))
hist(log2(l002.14[,3]),freq=FALSE,breaks=40, main="Euploid Line 2 Chr 14", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l002.14[,3])),col="green")
mc11 <- mean(log2(l002.14[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l002.14[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l002.14[,3]))
#mean is -0.8662616
sd(log2(l002.14[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l002.14[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l002.14[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l003.14[,3]),freq=FALSE,breaks=40, main="Euploid Line 3 Chr 14", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l003.14[,3])),col="green")
mc11 <- mean(log2(l003.14[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l003.14[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l003.14[,3]))
#mean is -0.8662616
sd(log2(l003.14[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l003.14[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l003.14[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l005.14[,3]),freq=FALSE,breaks=40, main="Euploid Line 5 Chr 14", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l005.14[,3])),col="green")
mc11 <- mean(log2(l005.14[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l005.14[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l005.14[,3]))
#mean is -0.8662616
sd(log2(l005.14[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l005.14[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l005.14[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l006.14[,3]),freq=FALSE,breaks=40, main="Euploid Line 6 Chr 14", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l006.14[,3])),col="green")
mc11 <- mean(log2(l006.14[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l006.14[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l006.14[,3]))
#mean is -0.8662616
sd(log2(l006.14[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l006.14[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l006.14[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l009.14[,3]),freq=FALSE,breaks=40, main="Euploid Line 9 Chr 14", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l009.14[,3])),col="green")
mc11 <- mean(log2(l009.14[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l009.14[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l009.14[,3]))
#mean is -0.8662616
sd(log2(l009.14[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l009.14[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l009.14[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l076.14[,3]),freq=FALSE,breaks=40, main="Trisomic Line 76 Chr 14", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l076.14[,3])),col="green")
mc11 <- mean(log2(l076.14[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l076.14[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l076.14[,3]))
#mean is -0.8662616
sd(log2(l076.14[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l076.14[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l076.14[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

###CHROMOSOME 10 

hist(log2(l002.10[,3]),freq=FALSE,breaks=40, main="Euploid Line 2 Chr 10", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l002.10[,3])),col="green")
mc11 <- mean(log2(l002.10[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l002.10[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l002.10[,3]))
#mean is -0.8662616
sd(log2(l002.10[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l002.10[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l002.10[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l003.10[,3]),freq=FALSE,breaks=40, main="Euploid Line 3 Chr 10", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l003.10[,3])),col="green")
mc11 <- mean(log2(l003.10[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l003.10[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l003.10[,3]))
#mean is -0.8662616
sd(log2(l003.10[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l003.10[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l003.10[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l005.10[,3]),freq=FALSE,breaks=40, main="Euploid Line 5 Chr 10", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l005.10[,3])),col="green")
mc11 <- mean(log2(l005.10[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l005.10[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l005.10[,3]))
#mean is -0.8662616
sd(log2(l005.10[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l005.10[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l005.10[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l006.10[,3]),freq=FALSE,breaks=40, main="Euploid Line 6 Chr 10", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l006.10[,3])),col="green")
mc11 <- mean(log2(l006.10[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l006.10[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l006.10[,3]))
#mean is -0.8662616
sd(log2(l006.10[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l006.10[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l006.10[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l009.10[,3]),freq=FALSE,breaks=40, main="Euploid Line 9 Chr 10", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l009.10[,3])),col="green")
mc11 <- mean(log2(l009.10[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l009.10[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l009.10[,3]))
#mean is -0.8662616
sd(log2(l009.10[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l009.10[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l009.10[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l066.10[,3]),freq=FALSE,breaks=40, main="Trisomic Line 66 Chr 10", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l066.10[,3])),col="green")
mc11 <- mean(log2(l066.10[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l066.10[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l066.10[,3]))
#mean is -0.8662616
sd(log2(l066.10[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l066.10[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l066.10[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)


###CHROMOSOME 15

hist(log2(l002.15[,3]),freq=FALSE,breaks=40, main="Euploid Line 2 Chr 15", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l002.15[,3])),col="green")
mc11 <- mean(log2(l002.15[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l002.15[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l002.15[,3]))
#mean is -0.8662616
sd(log2(l002.15[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l002.15[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l002.15[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l003.15[,3]),freq=FALSE,breaks=40, main="Euploid Line 3 Chr 15", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l003.15[,3])),col="green")
mc11 <- mean(log2(l003.15[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l003.15[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l003.15[,3]))
#mean is -0.8662616
sd(log2(l003.15[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l003.15[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l003.15[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l005.15[,3]),freq=FALSE,breaks=40, main="Euploid Line 5 Chr 15", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l005.15[,3])),col="green")
mc11 <- mean(log2(l005.15[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l005.15[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l005.15[,3]))
#mean is -0.8662616
sd(log2(l005.15[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l005.15[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l005.15[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l006.15[,3]),freq=FALSE,breaks=40, main="Euploid Line 6 Chr 15", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l006.15[,3])),col="green")
mc11 <- mean(log2(l006.15[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l006.15[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l006.15[,3]))
#mean is -0.8662616
sd(log2(l006.15[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l006.15[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l006.15[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l009.15[,3]),freq=FALSE,breaks=40, main="Euploid Line 9 Chr 15", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l009.15[,3])),col="green")
mc11 <- mean(log2(l009.15[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l009.15[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l009.15[,3]))
#mean is -0.8662616
sd(log2(l009.15[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l009.15[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l009.15[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l011.15[,3]),freq=FALSE,breaks=40, main="Trisomic Line 11 Chr 15", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l011.15[,3])),col="green")
mc11 <- mean(log2(l011.15[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l011.15[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l011.15[,3]))
#mean is -0.8662616
sd(log2(l011.15[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l011.15[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l011.15[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

###CHROMOSOME 16
par(mfrow=c(5,2))
hist(log2(l002.16[,3]),freq=FALSE,breaks=40, main="Euploid Line 2 Chr 16", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l002.16[,3])),col="green")
mc11 <- mean(log2(l002.16[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l002.16[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l002.16[,3]))
#mean is -0.8662616
sd(log2(l002.16[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l002.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l002.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l003.16[,3]),freq=FALSE,breaks=40, main="Euploid Line 3 Chr 16", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l003.16[,3])),col="green")
mc11 <- mean(log2(l003.16[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l003.16[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l003.16[,3]))
#mean is -0.8662616
sd(log2(l003.16[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l003.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l003.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l005.16[,3]),freq=FALSE,breaks=40, main="Euploid Line 5 Chr 16", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l005.16[,3])),col="green")
mc11 <- mean(log2(l005.16[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l005.16[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l005.16[,3]))
#mean is -0.8662616
sd(log2(l005.16[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l005.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l005.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l006.16[,3]),freq=FALSE,breaks=40, main="Euploid Line 6 Chr 16", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l006.16[,3])),col="green")
mc11 <- mean(log2(l006.16[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l006.16[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l006.16[,3]))
#mean is -0.8662616
sd(log2(l006.16[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l006.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l006.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l009.16[,3]),freq=FALSE,breaks=40, main="Euploid Line 9 Chr 16", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l009.16[,3])),col="green")
mc11 <- mean(log2(l009.16[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l009.16[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l009.16[,3]))
#mean is -0.8662616
sd(log2(l009.16[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l009.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l009.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l008.16[,3]),freq=FALSE,breaks=40, main="Trisomic Line 8 Chr 16", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l008.16[,3])),col="green")
mc11 <- mean(log2(l008.16[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l008.16[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l008.16[,3]))
#mean is -0.8662616
sd(log2(l008.16[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l008.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l008.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l031.16[,3]),freq=FALSE,breaks=40, main="Trisomic Line 31 Chr 16", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l031.16[,3])),col="green")
mc11 <- mean(log2(l031.16[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l031.16[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l031.16[,3]))
#mean is -0.8662616
sd(log2(l031.16[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l031.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l031.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l069.16[,3]),freq=FALSE,breaks=40, main="Trisomic Line 69 Chr 16", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l069.16[,3])),col="green")
mc11 <- mean(log2(l069.16[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l069.16[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l069.16[,3]))
#mean is -0.8662616
sd(log2(l069.16[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l069.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l069.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l112.16[,3]),freq=FALSE,breaks=40, main="Trisomic Line 112 Chr 16", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l112.16[,3])),col="green")
mc11 <- mean(log2(l112.16[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l112.16[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l112.16[,3]))
#mean is -0.8662616
sd(log2(l112.16[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l112.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l112.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

hist(log2(l141.16[,3]),freq=FALSE,breaks=40, main="Trisomic Line 141 Chr 16", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1.2))
lines(density(log2(l141.16[,3])),col="green")
mc11 <- mean(log2(l141.16[,3]))
mc11
abline(v = mc11, col = "blue")
skewness(log2(l141.16[,3]))
#distribution is left skewed
##skewness = -4.146241
mean(log2(l141.16[,3]))
#mean is -0.8662616
sd(log2(l141.16[,3]))
#sd is 0.9248915
#for line 152, which is trisomic for chr 1 
s11 <- skewness(log2(l141.16[,3]))
mylabel = bquote(italic("Skewness") == .(format(s11, digits = 3)))
text(x=2.2,y = 1.1, labels = mylabel,cex=.8)
sd11 <-sd(log2(l141.16[,3]))
mylabel2 = bquote(italic("SD") == .(format(sd11, digits = 3)))
text(x=2.2,y = .9, labels = mylabel2,cex=.8)
mylabel3 = bquote(italic("Mean") == .(format(mc11, digits = 3)))
text(x=2.2,y = 1, labels = mylabel3,cex=.8)

##boxplot and t-test for chromosome 1 trisomic
boxplot((log2(l018.1[,3])), (log2(euRatio.c1)),names=c("Trisomic", "Euploid"),ylab="log2(fold change)", col=c("turquoise", "purple"))
t.test((log2(l018.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
##


##this is for all trisomic lines merged together
hist(log2(anRatio.c1),freq=FALSE,breaks=40, main="Trisomic Chr 1", xlim=range(-3,3),xlab="log2(fold change)", ylim=range(0,1))
lines(density(log2(anRatio.c1)),na.rm=TRUE,col="green")
abline(v = mc1, col = "blue")


###library(ggplot2)
##gplot(anRatio.c1, aes(x = X)) + geom_histogram(aes(y = ..density..)) + geom_density()


hist(log2(anRatio.c5),breaks=40,main="Triploid Chr 5",xlim=range(-3,3),xlab="log2(fold change)")
skewness(log2(anRatio.c5))
#left skewed
#value is -2.89
mean(log2(anRatio.c5))
mean(anRatio.c5)
#mean is 1.14
sd(log2(anRatio.c5))

##NEED TO COME BACK AND FINISH THE REST OF THE CHROMOSOMES


##############################################################################################
#now lets look at line 11, which is 1n for c1

anRatio.1n.c1 <- l011.1[,3]
mean(log2(anRatio.1n.c1))
mean(log2(anRatio.1n.c1*1.05))
hist(log2(anRatio.1n.c1),breaks=40,main="Monosomic Chr 1",xlim=range(-3,3),xlab="log2(fold change)")
skewness(log2(anRatio.1n.c1))
##skewed left
##value is -4.146241

##############################################################################################
# let's look at the bottom X% of genes in the samples
#############
go1 <- NULL;go2 <- NULL;go3 <- NULL;go4 <- NULL;go5 <- NULL;go6 <- NULL;go7 <- NULL;go8 <- NULL;
go9 <- NULL;go10 <- NULL;go11 <- NULL;go12 <- NULL;go13 <- NULL;go14 <- NULL;go15 <- NULL;go16 <- NULL


########################
for(i in 1:16){ #for each go1 go2 go3 etc
  c <- NULL
  for (line in rep){ #for each line
    c <- cbind(c,as.character(getChrmRatio(line,i)[,3]))
    #assign(paste("go",i,sep=""),value=as.character(getChrmRatio(line,i)[,2]))
  }
  assign(paste("go",i,sep=""),value=c)
} #this gets the gene order on each line for each chromsome

#########################
colnames(go1) <- rep
colnames(go2) <- rep
colnames(go3) <- rep
colnames(go4) <- rep
colnames(go5) <- rep
colnames(go6) <- rep
colnames(go7) <- rep
colnames(go8) <- rep
colnames(go9) <- rep
colnames(go10) <- rep
colnames(go11) <- rep
colnames(go12) <- rep
colnames(go13) <- rep
colnames(go14) <- rep
colnames(go15) <- rep
colnames(go16) <- rep
##########################
#makes colnames X001, X002, etc
go1 <- data.frame(go1)
length(intersect(intersect(intersect(intersect(go1$X11[1:50],go1$X152[1:50]),go1$X18[1:50]),go1$X21[1:50]),go1$X7[1:50]))

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
e



########################################################################################################
#Cuffdiff P-values plotted against the ancestral FPKM values
########################################################################################################

#first need to read in the Cuffdiff data 
#first do this for one line, then make a loop to do it for each line
setwd = "/Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Redo/Cuffdiff_files/Diff"
samples <- c(1,2,3,4,5,6,7,8,9,11,18,21,29,31,49,50,59,61,66,69,76,77,112,115,123,141,152)

sample1 <- read.table("/Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Redo/Cuffdiff_files/gene_exp_S1.diff",header=TRUE)

plot(sample1$q_value,sample1$value_2,ylab="ancestor FPKM",xlab="FDR-adj p-value",main="Sample1")

list_of_titles <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6","Sample7","Sample8",
                    "Sample9","Sample11","Sample18","Sample21","Sample29","Sample31","Sample49","Sample50",
                    "Sample59","Sample61","Sample66","Sample69","Sample76","Sample77","Sample112","Sample115",
                    "Sample123","Sample141","Sample152")


# Plot separate ggplot figures in a loop.
library(ggplot2)

# Make list of variable names to loop over.
myFiles <- list.files(path="/Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Redo/Cuffdiff_files/Diff",pattern = "*.diff")
str(myFiles)
length(myFiles)
View(myFiles)
# Make plots.
plot_list = list()

#plot each sample in myFiles, then save to a .pdf in the working directory 
for (i in 1:length(myFiles)) {
  pdf(paste("FPKM_vs_p_val_",list_of_titles[i],".pdf",sep=""))
  sample <- read.table(myFiles[i],header=TRUE)
  sample <- data.frame(sample)
  qval <- sample$q_value
  val2 <- sample$value_2
  plot(qval,val2,ylab="ancestor FPKM",xlab="FDR-adj p-value", 
       xlim=c(0,1),ylim=c(0,10000),main=list_of_titles[i])
  dev.off() 
}


