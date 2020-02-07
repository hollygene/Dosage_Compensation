# need mean of log2 values for each sample 

  # make data frame with Line ---> 
                #    log2 mean
                #    chr # 


chr1Data.old <- data.frame(Line=c("01","02","03","04","05","06","07","08","09","011","015","028","088","0108","0119"), chrI=c((mean(log2(l001.1.old[,3]))),(mean(log2(l002.1.old[,3]))),(mean(log2(l003.1.old[,3]))),(mean(log2(l004.1.old[,3]))),(mean(log2(l005.1.old[,3]))),
                               (mean(log2(l006.1.old[,3]))),(mean(log2(l007.1.old[,3]))),(mean(log2(l008.1.old[,3]))),(mean(log2(l009.1.old[,3])))
                               ,(mean(log2(l011.1.old[,3]))),(mean(log2(l015.1[,3]))),(mean(log2(l028.1.old[,3]))),(mean(log2(l088.1[,3]))),(mean(log2(l108.1[,3]))),(mean(log2(l119.1[,3])))))


chr1Data <- data.frame(Line=c("1","2","3","4","5","7","8","9","11","18","21","31","49","59","61","66","69","76","77","29","50","112","115","117","123","141","152"), 
                       chrI=c((mean((log2(l001.1[,3])))),(mean((log2(l002.1[,3])))),(mean((log2(l003.1[,3])))),(mean((log2(l004.1[,3])))),(mean((log2(l005.1[,3])))),
                      (mean((log2(l007.1[,3])))),(mean((log2(l008.1[,3])))),(mean((log2(l009.1[,3])))),(mean((log2(l011.1[,3])))),
                      (mean((log2(l018.1[,3])))),(mean((log2(l021.1[,3])))),(mean((log2(l031.1[,3])))),(mean((log2(l049.1[,3])))),(mean((log2(l059.1[,3])))),
                      (mean((log2(l061.1[,3])))),(mean((log2(l066.1[,3])))),(mean((log2(l069.1[,3])))),(mean((log2(l076.1[,3])))),(mean((log2(l077.1[,3]))))
                    ,(mean((log2(l029.1[,3])))),(mean((log2(l050.1[,3])))),(mean((log2(l112.1[,3])))),(mean((log2(l115.1[,3])))),(mean((log2(l117.1[,3])))),
                                                            (mean((log2(l123.1[,3])))),(mean((log2(l141.1[,3])))),(mean((log2(l152.1[,3]))))))


chr1DataCombo <- rbind(chr1Data,chr1Data.old)



# attempt at making a loop for this

df <- data.frame(matrix(ncol=42, nrow=16, data=NA))
getMeans <- function{
  
chrData <- chrI=c((mean(log2(l001.1.old[,3]))),(mean(log2(l002.1.old[,3]))),(mean(log2(l003.1.old[,3]))),
                                  (mean(log2(l004.1.old[,3]))),(mean(log2(l005.1.old[,3]))),
                (mean(log2(l006.1.old[,3]))),(mean(log2(l007.1.old[,3]))),(mean(log2(l008.1.old[,3]))),(mean(log2(l009.1.old[,3])))
                ,(mean(log2(l011.1.old[,3]))),(mean(log2(l015.1[,3]))),(mean(log2(l028.1.old[,3]))),(mean(log2(l088.1[,3]))),
                (mean(log2(l108.1[,3]))),(mean(log2(l119.1[,3]))),(mean((log2(l001.1[,3])))),(mean((log2(l002.1[,3])))),(mean((log2(l003.1[,3])))),(mean((log2(l004.1[,3])))),(mean((log2(l005.1[,3])))),
                (mean((log2(l007.1[,3])))),(mean((log2(l008.1[,3])))),(mean((log2(l009.1[,3])))),(mean((log2(l011.1[,3])))),
                (mean((log2(l018.1[,3])))),(mean((log2(l021.1[,3])))),(mean((log2(l031.1[,3])))),(mean((log2(l049.1[,3])))),(mean((log2(l059.1[,3])))),
                (mean((log2(l061.1[,3])))),(mean((log2(l066.1[,3])))),(mean((log2(l069.1[,3])))),(mean((log2(l076.1[,3])))),(mean((log2(l077.1[,3]))))
                ,(mean((log2(l029.1[,3])))),(mean((log2(l050.1[,3])))),(mean((log2(l112.1[,3])))),(mean((log2(l115.1[,3])))),(mean((log2(l117.1[,3])))),
                (mean((log2(l123.1[,3])))),(mean((log2(l141.1[,3])))),(mean((log2(l152.1[,3]))))))

View(chrData)
chr1Data <- data.frame(Line=c("1","2","3","4","5","7","8","9","11","18","21","31","49","59","61","66","69","76","77","29","50","112","115","117","123","141","152"), 
                       chrI=c((mean((log2(l001.1[,3])))),(mean((log2(l002.1[,3])))),(mean((log2(l003.1[,3])))),(mean((log2(l004.1[,3])))),(mean((log2(l005.1[,3])))),
                              (mean((log2(l007.1[,3])))),(mean((log2(l008.1[,3])))),(mean((log2(l009.1[,3])))),(mean((log2(l011.1[,3])))),
                              (mean((log2(l018.1[,3])))),(mean((log2(l021.1[,3])))),(mean((log2(l031.1[,3])))),(mean((log2(l049.1[,3])))),(mean((log2(l059.1[,3])))),
                              (mean((log2(l061.1[,3])))),(mean((log2(l066.1[,3])))),(mean((log2(l069.1[,3])))),(mean((log2(l076.1[,3])))),(mean((log2(l077.1[,3]))))
                              ,(mean((log2(l029.1[,3])))),(mean((log2(l050.1[,3])))),(mean((log2(l112.1[,3])))),(mean((log2(l115.1[,3])))),(mean((log2(l117.1[,3])))),
                              (mean((log2(l123.1[,3])))),(mean((log2(l141.1[,3])))),(mean((log2(l152.1[,3]))))))


chr1DataCombo <- rbind(chr1Data,chr1Data.old)







# chr1DataCombo$Line <- as.factor(chr1DataCombo$Line)
View(chr1DataCombo)
chr1DataCombo$copy_number <- as.factor(chr1DataCombo$copy_number)
# chr 1 has 87 genes on it that I looked at
View(chr1DataCombo)
chr1DataCombo$Line=factor(chr1DataCombo$Line, c("1","2","3","4","5","7","8","9","11","18","21","31","49","59","61","66","69","76","77","01","02","03","04","05","06","07","08","09","011","015","028","29","50","088","0108","112","115","117","0119","123","141","152"))


p <- ggplot(chr1DataCombo, aes(x=Line, y=y,color=copy_number)) + 
  geom_boxplot() + theme_classic() + scale_color_manual(values=c("blue","gray","red")) 







### have this data! need to get mean not median 
y.mean.MA=NULL
rep.MA <- c("112","115","117","123","141","152","29","50")
#View(rep.MA)
for(line in rep.MA){
  for(chr in 1:16){
    y.mean.MA <- c(y.mean.MA,mean(getChrmRatio.MA(line,chr)[,3],na.rm=T))
  }
}

ratMatMed.mean.MA <- matrix(y.mean.MA,nrow=16,ncol=8)
colnames(ratMatMed.mean.MA) <- rep.MA
ratMatMed.mean.MA <- cbind(c(1:16),ratMatMed.mean.MA)
colnames(ratMatMed.mean.MA)[1] <- "chr"

#remove the ratios that are about 1
ratMatMedRed.mean.MA <- ratMatMed.mean.MA
ratMatMedRed.mean.MA[ratMatMedRed.mean.MA<1.35] <- NA

y.MA.mean.old=NULL
rep.MA.old <- c("001","002","003","004","005","006","007","008","009","011","015","028","088","108")
#View(rep.MA)
for(line in rep.MA.old){
  for(chr in 1:16){
    y.MA.mean.old <- c(y.MA.mean.old,mean(getChrmRatio.MA.old(line,chr)[,3],na.rm=T))
  }
}

ratMatMed.MA.mean.old <- matrix(y.MA.mean.old,nrow=16,ncol=14)
colnames(ratMatMed.MA.mean.old) <- rep.MA.old
ratMatMed.MA.mean.old <- cbind(c(1:16),ratMatMed.MA.mean.old)
colnames(ratMatMed.MA.mean.old)[1] <- "chr"

#remove the ratios that are about 1
ratMatMedRed.MA.mean.old <- ratMatMed.MA.mean.old
ratMatMedRed.MA.mean.old[ratMatMedRed.MA.mean.old<1.35] <- NA
#ratMatMedRed.MA[9,15] <- ratMatMedRed.MA[9,15]


#y.GC=NULL
#rep.GC <- c("11","18","1","21","2","31","3","49","4","59","5","61","66","69","6","76","77","7","8","9")
#for(line in rep.GC){
#  for(chrm in 1:16){
#    y.GC <- c(y.GC,median(getChrmRatio.GC(line,chrm)[,3],na.rm=T))
#  }
#}

#ratMatMed.GC <- matrix(y.GC,nrow=16,ncol=20)
#colnames(ratMatMed.GC) <- rep.GC
#ratMatMed.GC <- cbind(c(1:16),ratMatMed.GC)
#colnames(ratMatMed.GC)[1] <- "chr"


y.GC=NULL
colnames(sc_des.GC)
rep.GC <- c("11","18","1","2","31","3","49","4","5","61","69","76","77","7","8","9")
for(line in rep.GC){
  for(chrm in 1:16){
    y.GC<- c(y.GC,mean(getChrmRatio.GC(line,chrm)[,3],na.rm=T))
  }
}

ratMatMed.GC <- matrix(y.GC,nrow=16,ncol=16)
colnames(ratMatMed.GC) <- rep.GC
ratMatMed.GC <- cbind(c(1:16),ratMatMed.GC)
colnames(ratMatMed.GC)[1] <- "chr"
#remove the ratios that are about 1
ratMatMedRed.GC <- ratMatMed.GC
ratMatMedRed.GC[ratMatMedRed.GC<1.35] <- NA
#ratMatMedRed.MA[9,15] <- ratMatMedRed.MA[9,15]



#make a matrix that is the LOG2 RATIO
# log2 ratios 
yl.MA=NULL
for(line in rep.MA){
  for(chrm in 1:16){
    yl.MA <- c(yl.MA,mean(log2(getChrmRatio.MA(line,chrm)[,3])))
  }
}

warnings()
logMat.MA <- matrix(yl.MA,nrow=16,ncol=8)
colnames(logMat.MA) <- rep.MA
logMat.MA <- data.frame(logMat.MA)
chrms <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
logMat2.MA <-cbind(chrms,logMat.MA)
#write.csv(logMat2.MA, file="logMat.MA.csv")

#log 2 ratios 
yl.MA.old=NULL
for(line in rep.MA.old){
  for(chrm in 1:16){
    yl.MA.old <- c(yl.MA.old,mean(log2(getChrmRatio.MA.old(line,chrm)[,3])))
  }
}
length(rep.MA.old)
warnings()
logMat.MA.old <- matrix(yl.MA.old,nrow=16,ncol=14)
colnames(logMat.MA.old) <- rep.MA.old
logMat.MA.old <- data.frame(logMat.MA.old)
chrms <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
logMat2.MA.old <-cbind(chrms,logMat.MA.old)

# log2 ratios 
yl.GC=NULL
for(line in rep.GC){
  for(chrm in 1:16){
    yl.GC <- c(yl.GC,mean(log2(getChrmRatio.GC(line,chrm)[,3])))
  }
}
warnings()
length(rep.GC)
logMat.GC <- matrix(yl.GC,nrow=16,ncol=16)
colnames(logMat.GC) <- rep.GC
logMat.GC <- data.frame(logMat.GC)
chrms <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
logMat2.GC <-cbind(chrms,logMat.GC)

logMat2.GC
logMat2.MA
logMat2.MA.old
dim(logMatAll)
logMatAll <- cbind(logMat2.GC, logMat2.MA[,2:8], logMat2.MA.old[2:14])
# View(logMatAll)
# cols <- colnames(logMatAll[,2:37])
# cols
# library(reshape)
# logMatNew <- reshape(logMatAll,direction="long")
# 
# aes(x=chrms,color=copy_number))


p <- ggplot(logMatAll,aes(x=2:37,y=chrms)) 

m=matrix(1:12,3,4)
m

as.vector(m)

as.vector(t(m))


logMatTest <- as.vector(t(logMatAll[,2:37]))

length(logMatTest)
View(logMatTest)

chrms = rep(logMatAll[,1],each=36)
View(chrms)


(t(logMatAll[,1]))


plot(chrms,logMatTest)

p <- ggplot(chrms,logMatTest, aes(x=chrms,y=logMatTest))

write.csv(logMatAll, file="logMatAll.csv")




