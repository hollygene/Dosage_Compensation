file.names <- list.files(path = ".", pattern="*.txt")

for(i in 1:length(file.names)){

 file <- read.table(file.names[i], sep= "\t", header=TRUE)
 pdf(file= paste("Graphs/",file.names[i],".HS_Graph.pdf", sep=""))
 plot(file$GCEvents, col=as.factor(file$Chromosome), ylab="GCEvents",ylim=c(0,10), xlab="Bins (10kb)")
 dev.off()
 }