setwd("~/FinalCodeCNVs")

gs_file_names <- list.files("~/FinalCodeCNVs/Output_ProstateTrio_GenomeStudio")

gs_file_names <- gs_file_names[-length(gs_file_names)] # Delete the genomestudio_guide.doc


for(i in 1:length(gs_file_names)){
  
  cc <- rep("NULL",7) #7 = numebr of columns of the files
  cc[c(2,3,4,5,6,7)]<-NA
    
    
  output_gs_data <- read.table(paste("~/FinalCodeCNVs/Output_ProstateTrio_GenomeStudio/",gs_file_names[i],sep=""),
                               skip=9,#+(j-1)*num_each,
                               header=TRUE,
                               sep="\t",
                               colClasses=cc)

  sample_name <- unique(output_gs_data$Sample.Name)
  
  header <- read.table(paste("~/FinalCodeCNVs/Output_ProstateTrio_GenomeStudio/",gs_file_names[i],sep=""),
                       header=FALSE,
                       sep="\t",
                       nrow=8, fill=TRUE)
  
  for(j in 1:length(as.character(sample_name))){
    extra_line1 <- paste("File",j,"of",length(as.character(sample_name)), sep=" ")
    extra_line2 <- "[Data]"
    
    write.table(header, paste("~/FinalCodeCNVs/DataFiles/",as.character(sample_name[j]),".txt",sep=""),row.names=FALSE,col.names=FALSE, quote=FALSE)
    write.table(extra_line1, paste("~/FinalCodeCNVs/DataFiles/",as.character(sample_name[j]),".txt",sep=""),append=TRUE,row.names=FALSE,col.names=FALSE, quote=FALSE)
    write.table(extra_line2, paste("~/FinalCodeCNVs/DataFiles/",as.character(sample_name[j]),".txt",sep=""),append=TRUE,row.names=FALSE,col.names=FALSE, quote=FALSE)
        
    output_sample_data <- output_gs_data[which(as.character(output_gs_data$Sample.Name)==as.character(sample_name[j])),]
    
    write.table(output_sample_data, paste("~/FinalCodeCNVs/DataFiles/",as.character(sample_name[j]),".txt",sep=""),append=TRUE,row.names=FALSE, quote=FALSE)
    
  }

}











