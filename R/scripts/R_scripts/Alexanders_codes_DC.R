#DC with the library-size normalized counts

library(moments)
##########read the data
sc_counts <- read.table("C:/Users/John/Desktop/binf/cuffnorm_out_all/gene_fpkm_chrm.tsv",header=T)
colnames(sc_counts)[1] <- "chr"; colnames(sc_counts)[2] <- "gene"

sc_counts.z <- sc_counts 
sc_counts.z[,3:50] <- sc_counts.z[,3:50]+.5 #add .5 to dataset so the zeroes dont mess up

#we want to exclude genes that are very lowly expressed
gene.means <- rowMeans(sc_counts.z[,3:50]) #average across the sample
anc.means <- rowMeans(sc_counts.z[,47:50]) #average in the ancestor
#if the average fpkm is less than 5 in the ancestor or less than 5 in the overall, remove the gene
sc_counts.final <- sc_counts.z[-which(gene.means<5 | anc.means<5),]

#final matrix for ancestor and descendans
sc_anc <- sc_counts.final[,c(1,2,48,49,50)]
sc_des <- sc_counts.final[,1:47]

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
ancAvg <- rowMeans(sc_anc[,3:5],na.rm=T)
ancAvg <- cbind(sc_anc[,1:2],ancAvg)

#first write some functions that help automate the process
##########
#this gets the average FPKM for a given line from the sc_des matrix
getAvgFPKM <- function(line) { #for line as "###" string
  r1 <- paste(c("X",line,"_0"),collapse="") #triplicates are named X[LINE]_[rep]
  r2 <- paste(c("X",line,"_1"),collapse="")
  r3 <- paste(c("X",line,"_2"),collapse="")
  reps <- sc_des[,c("chr","gene",r1,r2,r3)] #get those 4 columns out of the descendant dataset
  reps[reps==0] <- NA #replace 0 values with NA as before
  desAvg <- rowMeans(reps[,c(-1,-2)],na.rm=T)
  desAvg <- cbind(sc_des[,c(1,2)],desAvg)
  #colnames(desAvg) <- c("chr","gene",paste(c(line,"_Avg"),collapse=""))
  return(desAvg) #returned is a matrix equivalent to ancAvg: chrm, gene, and line average FPKM 
}
##########
#this gets the fold change for the whole line
getFPKMRatio <- function(line) { #line must be a 3 digit character
  desAvg <- getAvgFPKM(line) #first call avg fpkm script, which returns a 3 column matrix
  ratioByGene <- desAvg[,3]/ancAvg[,3]
  rat <- cbind(sc_des[,1:2],ratioByGene) #add gene and chrm values
  colnames(rat) <- c("chrm","gene","ratio")
  return(rat)
}
##########
#this gets the fold change for just the specified chromosome
getChrmRatio <- function(line, chr){ #input line as "XXX", chr as a number
  rat <- getFPKMRatio(line) #get list of ratios which includes chrm and gene
  rat.c <- rat[(rat[,1]==chr),] #get just the rows for the appropriate chrom
  rat.cs <- rat.c[order(rat.c[,3]),]#order by column 3
  return(rat.cs)
}
##########


#now we can look at the fold changes for each chromosome and each line

#make a matrix of median ratio:
y=NULL
rep <- c("001","002","003","004","005","006","007","008","009","011","015","028","088","108","119")
for(line in rep){
  for(chrm in 1:16){
    y <- c(y,median(getChrmRatio(line,chrm)[,3],na.rm=T))
  }
}

#make a matrix out of it
ratMatMed <- matrix(y,nrow=16,ncol=15)
colnames(ratMatMed) <- rep
ratMatMed <- cbind(c(1:16),ratMatMed)

#remove the ratios that are about 1
ratMatMedRed <- ratMatMed
ratMatMedRed[ratMatMed<1.35] <- NA
ratMatMedRed[9,15] <- ratMatMed[9,15]


#make a matrix that is the LOG2 RATIO
yl=NULL
for(line in rep){
  for(chrm in 1:16){
    yl <- c(yl,median(log2(getChrmRatio(line,chrm)[,3])))
  }
}

logMat <- matrix(yl,nrow=16,ncol=15)
colnames(logMat) <- c("001","002","003","004","005","006","007","008","009","011","015","028","088","108","119")
logMat <- data.frame(logMat)


########### compare to non amplified chromosome ratios
l001.9 <- getChrmRatio("001",9)
l002.9 <- getChrmRatio("002",9)
l003.9 <- getChrmRatio("003",9)
l004.9 <- getChrmRatio("004",9)
l005.9 <- getChrmRatio("005",9)
l006.9 <- getChrmRatio("006",9)
l007.9 <- getChrmRatio("007",9)
l008.9 <- getChrmRatio("008",9)
l011.9 <- getChrmRatio("011",9)
l028.9 <- getChrmRatio("028",9)

#line 9 is 3n for c14
l009.14 <- getChrmRatio("009",14)
#line 15, 88, 119 are 3n for c9
l015.9 <- getChrmRatio("015",9)
l088.9 <- getChrmRatio("088",9)
l119.9 <- getChrmRatio("119",9)

anRatio.c9 <- c(l015.9[,3],l088.9[,3],l119.9[,3])
euRatio.c9 <- c(l001.9[,3],l002.9[,3],l003.9[,3],l004.9[,3],l005.9[,3],l006.9[,3],
                l007.9[,3],l008.9[,3],l011.9[,3],l028.9[,3])

#make a histogram of log ratio of all genes across 
#euploid lines
hist(log2(euRatio.c9),breaks=40,xlab="log2(Fold Change)",xlim=range(-3,3),main="Euploid Chr 9")
skewness(log2(euRatio.c9))
#distribution is left skewed

mean(log2(euRatio.c9))
#mean is .026, close to 0
sd(log2(euRatio.c9))

#do the same for the aneuploid lines
hist(log2(anRatio.c9),breaks=20,main="Triploid Chr 9",xlim=range(-3,3),xlab="log2(fold change)")
skewness(log2(anRatio.c9))
#right skewed
mean(log2(anRatio.c9))
mean(anRatio.c9)
#mean is .548, a bit less than .58 expectation
sd(log2(anRatio.c9))
#sd in line with euploid lines


#now lets look at line 9, which is 3n for c14

anRatio.c14 <- l009.14[,3]
mean(log2(anRatio.c14))
mean(log2(anRatio.c14*1.05))

# let's look at the bottom X% of genes in the samples
#############
go1 <- NULL;go2 <- NULL;go3 <- NULL;go4 <- NULL;go5 <- NULL;go6 <- NULL;go7 <- NULL;go8 <- NULL;
go9 <- NULL;go10 <- NULL;go11 <- NULL;go12 <- NULL;go13 <- NULL;go14 <- NULL;go15 <- NULL;go16 <- NULL


########################
for(i in 1:16){ #for each go1 go2 go3 etc
  c <- NULL
  for (line in rep){ #for each line
    c <- cbind(c,as.character(getChrmRatio(line,i)[,2]))
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
go9 <- data.frame(go9)
length(intersect(intersect(go9$X015[1:50],go9$X088[1:50]),go9$X119[1:50]))
#we get 32 shared genes in the bottom 75 genes of the three triploid
#this is way more than expected, but not atypical for pairings of 3 euploid lines


#if we compare sets of 3 of the euploids, what do we get?

#make the list of reps we should choose from
xrep <- c(1:8,10,12)

#now make a loop to randomly sample these
overlap <- NULL
for(i in 1:1000){
triple <- sample(xrep,3,replace=F)
bounds <- 1:75
overlap[i] <- length(intersect(intersect(go9[bounds,triple[1]],go9[bounds,triple[2]]),go9[bounds,triple[3]]))
}
mean(overlap)
max(overlap)
min(overlap)
