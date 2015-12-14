
args <- commandArgs(trailingOnly = TRUE)
print(args)


library(DNAcopy)

#sapply(1:npatients, function(i){
  
  print(paste("Processing patient", args[1], sep=" "))
  Data <- read.table(paste("TB" ,args[1],".txt",sep=""),
                     header=TRUE, 
                     sep="\t")
  
  Data$Chr <- gsub("X",23,Data$Chr)
  Data$Chr <- gsub("Y",24,Data$Chr)
  Data$Chr <- as.numeric(Data$Chr) 
  
  cna_object <- CNA(as.numeric(as.character(Data[,5])),chrom=as.numeric(as.character(Data$Chr)),maploc=as.numeric(as.character(Data$Position)),data.type="logratio")
  
  segment_cna <- segment(cna_object,verbose=2)
  
  write.table(segment_cna$out, 
              file=paste("segment_TB",args[1],".txt",sep=""),
              quote=FALSE, 
              row.names=FALSE, 
              sep='\t')
  
  save(segment_cna,file=paste("segment_object_TB",args[1],".RData",sep=""))
  
#})
