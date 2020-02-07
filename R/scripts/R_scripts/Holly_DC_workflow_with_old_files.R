#Holly's DC scripts, modified from DCfinal.R
#DC with the library-size normalized counts

###Rate of aneuploidy for GC and MA data 
#are they statistically significantly different?
GC_rate <- 1.9969e-4 
MA_rate <- 9.7e-5
rate <- c(GC_rate,MA_rate)
chisq.test(rate)

library(moments)
##########read the data
sc_counts <- read.csv("/Users/hollymcqueary/Dropbox/McQueary/Dosage-Compensation/Old/genes_fpkm.chrs.csv",header=T)
###colnames(sc_counts) <- as.character(unlist(sc_counts[1,]))
###sc_counts = sc_counts[-1,]
class(sc_counts)
colnames(sc_counts)[1] <- "chr"; colnames(sc_counts)[2] <- "gene"
sc_counts.z <- sc_counts 
sc_counts.z[,3:92] <- sc_counts.z[,3:92]+.5 #add .5 to dataset so the zeroes dont mess up


warnings()
#we want to exclude genes that are very lowly expressed
gene.means <- rowMeans(sc_counts.z[,3:92]) #average across the sample
anc.means <- rowMeans(sc_counts.z[,87:92]) #average in the ancestor
#if the average fpkm is less than 5 in the ancestor or less than 5 in the overall, remove the gene
sc_counts.final <- sc_counts.z[-which(gene.means<5 | anc.means<5),]

#final matrix for ancestor and descendans
sc_anc <- sc_counts.final[,c(1,2,87,88,89,90,91,92)]
sc_des <- sc_counts.final[,1:86]

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
ancAvg <- rowMeans(sc_anc[,3:6],na.rm=T)
ancAvg <- cbind(sc_anc[,1:2],ancAvg)
#ancMin <- min(sc_anc[,3:6],na.rm=T) 
#ancMin <- cbind(sc_anc[,1:2],ancMin)
#ancMax <- max(sc_anc[,3:6],na.rm=T)
#ancMax <- cbind(sc_anc[,1:2],ancMax)

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
getFPKMRatio <- function(line) { #line must be a 3 digit character
  desAvg <- getAvgFPKM(line) #first call avg fpkm script, which returns a 3 column matrix
  ratioByGene <- desAvg[,3]/ancAvg[,3]
  rat <- cbind(sc_des[,1:2],ratioByGene) #add gene and chrm values
  colnames(rat) <- c("chrm","gene","ratio")
  return(rat)
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
getChrmRatio <- function(line, chr){ #input line as "XXX", chr as a number
  rat <- getFPKMRatio(line) #get list of ratios which includes chrm and gene
  rat.c <- rat[(rat[,1]==chr),] #get just the rows for the appropriate chrom
  rat.cs <- rat.c[order(rat.c[,3]),]#order by column 3
  return(rat.cs)
}
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

y=NULL
rep <- c("112","115","117","11","123","141","152","18","1","21","29","2","31","3","49","4","50","59","5","61","66","69","6","76","77","7","8","9")
for(line in rep){
  for(chrm in 1:16){
    y <- c(y,median(getChrmRatio(line,chrm)[,3],na.rm=T))
  }
}


#make a matrix out of it
ratMatMed <- matrix(y,nrow=16,ncol=28)
colnames(ratMatMed) <- rep
ratMatMed <- cbind(c(1:16),ratMatMed)

#remove the ratios that are about 1
ratMatMedRed <- ratMatMed
ratMatMedRed[ratMatMed<1.35] <- NA
ratMatMedRed[9,15] <- ratMatMed[9,15]

#make a matrix out of it- min 
#ratMatMedMin <- matrix(ymin,nrow=16,ncol=28)
#colnames(ratMatMedMin) <- rep
#ratMatMedMin <- cbind(c(1:16),ratMatMedMin)


#make a matrix out of it- min 
#ratMatMedMax <- matrix(ymax,nrow=16,ncol=28)
#colnames(ratMatMedMax) <- rep
#ratMatMedMax <- cbind(c(1:16),ratMatMedMax)



#make a matrix that is the LOG2 RATIO
yl=NULL
for(line in rep){
  for(chrm in 1:16){
    yl <- c(yl,median(log2(getChrmRatio(line,chrm)[,3])))
  }
}

logMat <- matrix(yl,nrow=16,ncol=28)
colnames(logMat) <- c("112","115","117","11","123","141","152","18","1","21","29","2","31","3","49","4","50","59","5","61","66","69","6","76","77","7","8","9")
logMat <- data.frame(logMat)
chrms <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
logMat2 <-cbind(chrms,logMat)

##make a matrix of the standard deviations for each line and each chromosome
ysd=NULL
for(line in rep){
  for(chrm in 1:16){
    ysd <- c(ysd,sd(log2(getChrmRatio(line,chrm)[,3])))
  }
}

logMatsd <- matrix(ysd,nrow=16,ncol=28)
colnames(logMatsd) <- c("112","115","117","11","123","141","152","18","1","21","29","2","31","3","49","4","50","59","5","61","66","69","6","76","77","7","8","9")
logMatsd <- data.frame(logMatsd)
chrms <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
logMatsd2 <-cbind(chrms,logMatsd)




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
X112 <- subset(logMat, select=c("X112"))
#make a table with the chromosome numer in there too 
X112c <- cbind(chrms,X112)
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
plot(X112~chrms, data=X112,ylim=c(-1.5,1.5),xlab="Chromosome", ylab="log2(fold change)",type="p",main="MA Line 112",pch=16,col="darkturquoise")
#this is just an example plot script
####plot(1:2, y.bar, ylim=c(0, 25), xlim=c(0.5, 2.5), xaxt="n", pch=16,
     ###xlab="", ylab="Tree density", cex=1.5)
##this makes the x axis have the chromosome # labels
axis(1, 1:16, c("1", "2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"))
##next I would like to figure out how to put error bars, but it is not completely necessary right now 
##add error bars
X112.sd <- subset(logMatsd2, select=c("X112"))
X112.sd1 <- X112.sd[1,]
arrows(1:2,X112-X112.sd1, 1:2, X112+X112.sd1, code=3,angle=90,length=0.05)
arrows(x, avg-sdev, x, avg+sdev, length=0.05, angle=90, code=3)
###arrows(1:2, X112-X112min, 1:2, X112+X112max, code=3, angle=90, length=0.05)
x <- 1:16
plot(X112~chrms, data=X112,
     ylim=range(c(X112-X112.sd, X112+X112.sd)),
     pch=16, xlab="Chromosome", ylab="log2(fold change",
     main="MA Line 112", col="darkturquoise"
)

###arrows(x, X112-X112.sd, x, X112+X112.sd, length=0.05, angle=90, code=3)

########### compare to non amplified chromosome ratios
##starting with chromosome I, going in order
##going to need these for all of the chromosomes that have a CNV, but can come back to this later
###probably an easier way to do this, but will come back to this later

l002.1 <- getChrmRatio("2",1)
l003.1 <- getChrmRatio("3",1)
l005.1 <- getChrmRatio("5",1)
l006.1 <- getChrmRatio("6",1)
l009.1 <- getChrmRatio("9",1)

l002.2 <- getChrmRatio("2",2)
l003.2 <- getChrmRatio("3",2)
l005.2 <- getChrmRatio("5",2)
l006.2 <- getChrmRatio("6",2)
l009.2 <- getChrmRatio("9",2)

l002.3 <- getChrmRatio("2",3)
l003.3 <- getChrmRatio("3",3)
l005.3 <- getChrmRatio("5",3)
l006.3 <- getChrmRatio("6",3)
l009.3 <- getChrmRatio("9",3)

l002.4 <- getChrmRatio("2",4)
l003.4 <- getChrmRatio("3",4)
l005.4 <- getChrmRatio("5",4)
l006.4 <- getChrmRatio("6",4)
l009.4 <- getChrmRatio("9",4)

l002.5 <- getChrmRatio("2",5)
l003.5 <- getChrmRatio("3",5)
l005.5 <- getChrmRatio("5",5)
l006.5 <- getChrmRatio("6",5)
l009.5 <- getChrmRatio("9",5)

l002.6 <- getChrmRatio("2",6)
l003.6 <- getChrmRatio("3",6)
l005.6 <- getChrmRatio("5",6)
l006.6 <- getChrmRatio("6",6)
l009.6 <- getChrmRatio("9",6)

l002.7 <- getChrmRatio("2",7)
l003.7 <- getChrmRatio("3",7)
l005.7 <- getChrmRatio("5",7)
l006.7 <- getChrmRatio("6",7)
l009.7 <- getChrmRatio("9",7)

l002.8 <- getChrmRatio("2",8)
l003.8 <- getChrmRatio("3",8)
l005.8 <- getChrmRatio("5",8)
l006.8 <- getChrmRatio("6",8)
l009.8 <- getChrmRatio("9",8)

l002.9 <- getChrmRatio("2",9)
l003.9 <- getChrmRatio("3",9)
l005.9 <- getChrmRatio("5",9)
l006.9 <- getChrmRatio("6",9)
l009.9 <- getChrmRatio("9",9)

l002.10 <- getChrmRatio("2",10)
l003.10 <- getChrmRatio("3",10)
l005.10 <- getChrmRatio("5",10)
l006.10 <- getChrmRatio("6",10)
l009.10 <- getChrmRatio("9",10)

l002.11 <- getChrmRatio("2",11)
l003.11 <- getChrmRatio("3",11)
l005.11 <- getChrmRatio("5",11)
l006.11 <- getChrmRatio("6",11)
l009.11 <- getChrmRatio("9",11)

l002.12 <- getChrmRatio("2",12)
l003.12 <- getChrmRatio("3",12)
l005.12 <- getChrmRatio("5",12)
l006.12 <- getChrmRatio("6",12)
l009.12 <- getChrmRatio("9",12)

l002.13 <- getChrmRatio("2",13)
l003.13 <- getChrmRatio("3",13)
l005.13 <- getChrmRatio("5",13)
l006.13 <- getChrmRatio("6",13)
l009.13 <- getChrmRatio("9",13)

l002.14 <- getChrmRatio("2",14)
l003.14 <- getChrmRatio("3",14)
l005.14 <- getChrmRatio("5",14)
l006.14 <- getChrmRatio("6",14)
l009.14 <- getChrmRatio("9",14)

l002.15 <- getChrmRatio("2",15)
l003.15 <- getChrmRatio("3",15)
l005.15 <- getChrmRatio("5",15)
l006.15 <- getChrmRatio("6",15)
l009.15 <- getChrmRatio("9",15)

l002.16 <- getChrmRatio("2",16)
l003.16 <- getChrmRatio("3",16)
l005.16 <- getChrmRatio("5",16)
l006.16 <- getChrmRatio("6",16)
l009.16 <- getChrmRatio("9",16)

#this works, just gives an interesting curve 
plot((log2(l002.1[,3])),ylab="log2(fold change)", xlab="Chromosome")
library(ggplot2)
#autoplot((log2(l002.1[,3])),ylab="log2(fold change)", xlab="Chromosome")

#ancestor ratios
ancAvg <- rowMeans(sc_anc[,3:6],na.rm=T)
ancAvg <- cbind(sc_anc[,1:2],ancAvg)

##############################################################################################

#CHROMOSOME 1 
#line 11 is 1n for c1
l011.1 <- getChrmRatio("11",1)
#line 152,18,21,7 are 3n for c1
l152.1 <- getChrmRatio("152",1)
l018.1 <- getChrmRatio("18",1)
l021.1 <- getChrmRatio("21",1)
l007.1 <- getChrmRatio("7",1)

#RATIOS FOR CHROMOSOME 1 
anRatio.c1 <- c(l152.1[,3],l018.1[,3],l021.1[,3],l007.1[,3])
euRatio.c1 <- c(l002.1[,3],l003.1[,3],l005.1[,3],l006.1[,3],l009.1[,3])

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

boxplot((log2(euRatio.c1)),(log2(l007.1[,3])),names=c( "Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan", "deeppink1"),main="Line 7",las=3)
t.test((log2(l007.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p3 <- 0.0006221
mylabel = bquote(italic(p) == .(format(p3, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=2.4, "*", pos=3, cex=1.5, col="red")

boxplot((log2(euRatio.c1)),(log2(l018.1[,3])),names=c( "Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan", "deeppink1"),main="Line 18",las=3)
t.test((log2(l018.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p3 <- 0.004569
mylabel = bquote(italic(p) == .(format(p3, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=2.4, "*", pos=3, cex=1.5, col="red")

boxplot((log2(euRatio.c1)),(log2(l011.1[,3])),names=c( "Disomic", "Monosomic"),ylab="log2(fold change)", col=c("cyan", "deeppink1"),main="Line 11",las=3)
t.test((log2(l011.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p3 <-  0.002302
mylabel = bquote(italic(p) == .(format(p3, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=2.4, "*", pos=3, cex=1.5, col="red")

boxplot((log2(euRatio.c1)),(log2(l021.1[,3])),names=c( "Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan", "deeppink1"),main="Line 21",las=3)
t.test((log2(l021.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p3 <-  0.1214
mylabel = bquote(italic(p) == .(format(p3, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)


boxplot((log2(euRatio.c1)),(log2(l152.1[,3])),names=c( "Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan", "deeppink1"),main="Line 152",las=3)
t.test((log2(l152.1[,3])), (log2(euRatio.c1)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p3 <-  1.546e-06
mylabel = bquote(italic(p) == .(format(p3, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
text(x=2, y=2.55, "*", pos=3, cex=1.5, col="red")

##############################################################################################
#CHROMOSOME 5 
#line 117 is 3n for c5
l117.5 <- getChrmRatio("117",5)
l049.5 <- getChrmRatio("49",5)
l004.5 <- getChrmRatio("4",5)
l050.5 <- getChrmRatio("50",5)

#RATIOS FOR CHROMOSOME 5
anRatio.c5 <- c(l117.5[,3],l049.5[,3],l004.5[,3],l050.5[,3])
euRatio.c5 <- c(l002.5[,3],l003.5[,3],l005.5[,3],l006.5[,3],l009.5[,3])
anRatio.c5.no50 <- c(l117.5[,3],l049.5[,3],l004.5[,3])
anRatio.c5.no4 <- c(l117.5[,3],l049.5[,3],l050.5[,3])

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


##significantly different
boxplot( (log2(euRatio.c5)),(log2(l004.5[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine", "mediumvioletred"),main="Line 4",las=3)
t.test((log2(l004.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p4 <- 1.414e-14
mylabel = bquote(italic(p) == .(format(p4, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

##significantly different
boxplot( (log2(euRatio.c5)),(log2(l050.5[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine", "mediumvioletred"),main="Line 50 Chr 5",las=3)
t.test((log2(l050.5[,3])), (log2(euRatio.c5)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <- 0.4638
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
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

#Line 31
boxplot( (log2(euRatio.c7)),(log2(l031.7[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan", "violetred1"),main="Chr 7 Gene Expression Line 31", las=3)
t.test((log2(l031.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
p5 <- 0.8172
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

#line 61
boxplot( (log2(euRatio.c7)),(log2(l061.7[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan", "violetred1"),main="Chr 7 Gene Expression Line 61", las=3)
t.test((log2(l061.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
p5 <- "< 2.2e-16"
abline(h=0,lty=3)
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

#line 115
boxplot( (log2(euRatio.c7)),(log2(l115.7[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan", "violetred1"),main="Chr 7 Gene Expression Line 115", las=3)
t.test((log2(l115.7[,3])), (log2(euRatio.c7)), paired=FALSE, var.equal=TRUE)
p5 <- 0.1992
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

#Line76
boxplot((log2(euRatio.c9)),(log2(l076.9[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan3", "mediumvioletred"),main="Chr 9 Gene Expression Line 76", las=3)
t.test((log2(l076.9[,3])), (log2(euRatio.c9)), paired=FALSE, var.equal=TRUE)
p5 <- "< 2.2e-16"
abline(h=0,lty=3)
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

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
p5 <-0.4511
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
p5 <- 0.245
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

#line18 
boxplot( (log2(euRatio.c12)),(log2(l018.12[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine2", "mediumvioletred"),main="Chr 12 Gene Expression: Line 18",las=3)
t.test((log2(l018.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-"< 2.2e-16"
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

#line77 
boxplot( (log2(euRatio.c12)),(log2(l077.12[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine2", "mediumvioletred"),main="Chr 12 Gene Expression: Line 77",las=3)
t.test((log2(l077.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-"< 2.2e-16"
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

#line123 
boxplot( (log2(euRatio.c12)),(log2(l123.12[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine2", "mediumvioletred"),main="Chr 12 Gene Expression: Line 123",las=3)
t.test((log2(l123.12[,3])), (log2(euRatio.c12)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-0.3206
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

##############################################################################################

####Chromosome 14

l076.14 <- getChrmRatio("76",14)


#RATIOS FOR CHROMOSOME 10
anRatio.c14 <- c(l076.14[,3])
euRatio.c14 <- c(l002.14[,3],l003.14[,3],l005.14[,3],l006.14[,3],l009.14[,3])

boxplot( (log2(euRatio.c14)),(log2(l076.14[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("aquamarine3", "deeppink1"),main="Chr 14 Gene Expression: Line 76",las=3)
t.test((log2(l076.14[,3])), (log2(euRatio.c14)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-"< 2.2e-16"
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)


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

##############################################################################################

###Chromosome 16 
l008.16 <- getChrmRatio("8",16)
l031.16 <- getChrmRatio("31",16)
l069.16 <- getChrmRatio("69",16)
l112.16 <- getChrmRatio("112",16)
l141.16 <- getChrmRatio("141",16)

#RATIOS FOR CHROMOSOME 10
anRatio.c16 <- c(l008.16[,3],l031.16[,3],l069.16[,3],l112.16[,3],l141.16[,3])
euRatio.c16 <- c(l002.16[,3],l003.16[,3],l005.16[,3],l006.16[,3],l009.16[,3])
#line 8
boxplot( (log2(euRatio.c16)),(log2(l008.16[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan3", "mediumorchid3"),main="Chr 16 Gene Expression: Line 8",las=3)
t.test((log2(l008.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-"< 2.2e-16"
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

#line31
boxplot( (log2(euRatio.c16)),(log2(l031.16[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan3", "mediumorchid3"),main="Chr 16 Gene Expression: Line 31",las=3)
t.test((log2(l031.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-0.212
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

#line69
boxplot( (log2(euRatio.c16)),(log2(l069.16[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan3", "mediumorchid3"),main="Chr 16 Gene Expression: Line 69",las=3)
t.test((log2(l069.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-0.1216
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

#line112
boxplot( (log2(euRatio.c16)),(log2(l112.16[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan3", "mediumorchid3"),main="Chr 16 Gene Expression: Line 112",las=3)
t.test((log2(l112.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-"< 2.2e-16"
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)

#line141
boxplot( (log2(euRatio.c16)),(log2(l141.16[,3])),names=c("Disomic", "Trisomic"),ylab="log2(fold change)", col=c("cyan3", "mediumorchid3"),main="Chr 16 Gene Expression: Line 141",las=3)
t.test((log2(l141.16[,3])), (log2(euRatio.c16)), paired=FALSE, var.equal=TRUE)
abline(h=0,lty=3)
p5 <-"< 2.2e-16"
mylabel = bquote(italic(p) == .(format(p5, digits = 9)))
text(x=2.3,y = 2.5, labels = mylabel,cex=.8)
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
