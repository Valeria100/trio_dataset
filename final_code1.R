Arraystolookat <- read.csv("Arraystolookat.csv", header=TRUE, sep=",")

#----------------------------------------------------------------------------------------------------------------------#
# First check if the patients name is the same and is in the same order
pnames <- matrix(0,nrow(Arraystolookat),6)
pnames[,1] <- sapply(strsplit(as.character(Arraystolookat$CC1), "_"), function(i) ifelse(length(i)!=0, i[1], NA)) 
pnames[,2] <- sapply(strsplit(as.character(Arraystolookat$CC2), "_"), function(i) ifelse(length(i)!=0, i[1], NA))
pnames[,3] <- sapply(strsplit(as.character(Arraystolookat$CC3), "_"), function(i) ifelse(length(i)!=0, i[1], NA))
pnames[,4] <- sapply(strsplit(as.character(Arraystolookat$CC4), "_"), function(i) ifelse(length(i)!=0, i[1], NA))
pnames[,5] <- sapply(strsplit(as.character(Arraystolookat$ICGC1), " "), function(i) ifelse(length(i)!=0, i[1], NA)) 
pnames[,6] <- sapply(strsplit(as.character(Arraystolookat$ICGC2), " "), function(i) ifelse(length(i)!=0, i[1], NA)) 

pnames[,1] <- gsub("[.]0",".",pnames[,1])
pnames[,2] <- gsub("[.]0",".",pnames[,2])
pnames[,3] <- gsub("[.]0",".",pnames[,3])
pnames[,4] <- gsub("[.]0",".",pnames[,4])

all.equal(pnames[,1],pnames[,2])#TRUE
all.equal(pnames[!is.na(pnames[,3]),1], pnames[!is.na(pnames[,3]),3])#TRUE
all.equal(pnames[!is.na(pnames[,4]),1], pnames[!is.na(pnames[,4]),4])#TRUE
all.equal(pnames[!is.na(pnames[,5]),1], pnames[!is.na(pnames[,5]),5])#TRUE
all.equal(pnames[!is.na(pnames[,6]),1], pnames[!is.na(pnames[,6]),6])#TRUE
#----------------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------------#
#CamCap B, T or G in a matrix
CC <- matrix(0,nrow(Arraystolookat),4)
CC[,1] <- sapply(strsplit(as.character(Arraystolookat$CC1), "_"), function(i) ifelse(length(i)!=0, i[length(i)], NA)) 
# length(i) - because I want the last element of the splitting (not always the 3rd one!)
# ifelse - because some cells are empty 
CC[,2] <- sapply(strsplit(as.character(Arraystolookat$CC2), "_"), function(i) ifelse(length(i)!=0, i[length(i)], NA))
CC[,3] <- sapply(strsplit(as.character(Arraystolookat$CC3), "_"), function(i) ifelse(length(i)!=0, i[length(i)], NA))
CC[,4] <- sapply(strsplit(as.character(Arraystolookat$CC4), "_"), function(i) ifelse(length(i)!=0, i[length(i)], NA))

# ICGC Blood or v* in a matrix
ICGC <- matrix(0,nrow(Arraystolookat),2)

ICGC[,1] <- sapply(strsplit(as.character(Arraystolookat$ICGC1), " "), function(i) ifelse(length(i)!=0, tolower(i[length(i)]), NA)) 
ICGC[,2] <- sapply(strsplit(as.character(Arraystolookat$ICGC2), " "), function(i) ifelse(length(i)!=0, tolower(i[length(i)]), NA)) 

ICGC[56,2] <- "blood" # here there was a typo - blood is missing in the name

ICGC <- gsub("^v[0-9]{1,2}","tumour",ICGC)


# CamCam file names in a matrix
CCID <- matrix(0,nrow(Arraystolookat),4)
CCID[,1] <- sapply(strsplit(as.character(Arraystolookat$CCID1), "/"), function(i) ifelse(length(i)!=0, i[1], NA)) 
# length(i) - because I want the last element of the splitting (not always the 3rd one!)
# ifelse - because some cells are empty 
CCID[,2] <- sapply(strsplit(as.character(Arraystolookat$CCID1), "/"), function(i) ifelse(length(i)>=2, i[2], NA))
CCID[,3] <- sapply(strsplit(as.character(Arraystolookat$CCID1), "/"), function(i) ifelse(length(i)>=3, i[3], NA))
CCID[,4] <- sapply(strsplit(as.character(Arraystolookat$CCID1), "/"), function(i) ifelse(length(i)>=4, i[4], NA))

#Replace those files that have been rerun and you can see in the other column
CCID2 <- strsplit(gsub("N/A", "NA", as.character(Arraystolookat$CCID2)),"/")
CCID2replace <- sapply(CCID2, function(i) which(i!="NA"))
ccid2_replace <- sapply(CCID2replace, function(i) ifelse(length(i)==0, 0, i ))
ccid2_replace_ind <- which(ccid2_replace!=0)

for(i in 1:length(ccid2_replace_ind)){
  CCID[ccid2_replace_ind[i],ccid2_replace[ccid2_replace_ind[i]]] <- CCID2[[ccid2_replace_ind[i]]][ccid2_replace[ccid2_replace_ind[i]]]
}


ICGCID <- matrix(0,nrow(Arraystolookat),2)
ICGCID[,1] <- sapply(strsplit(as.character(Arraystolookat$ICGCID1), "/"), function(i) ifelse(length(i)!=0, i[1], NA)) 
ICGCID[,2] <- sapply(strsplit(as.character(Arraystolookat$ICGCID1), "/"), function(i) ifelse(length(i)>=2, i[2], NA))
ICGCID[ICGCID=="NA"] <- NA
#----------------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------------#
# List the filenames for tumour, benign and germline for both camcap and icgc
source("functions_list.R")

tumour_filenames_camcap <- list_filenames("T",CC[,-4],CCID)
benign_filenames_camcap <- list_filenames("B",CC[,-4],CCID)
germline_filenames_camcap <- list_filenames("G",CC[, -4],CCID) # Ignore the 4th column. Only one repeated sample of germline

tumour_filenames_icgc <- list_filenames("tumour",ICGC,ICGCID)
germline_filenames_icgc <- list_filenames("blood",ICGC,ICGCID)


table_complete_names <- matrix(NA, nrow(Arraystolookat), 6)
table_complete_names[,1] <- pnames[,1]
table_complete_names[,2] <- tumour_filenames_camcap
table_complete_names[,3] <- benign_filenames_camcap
table_complete_names[,4] <- germline_filenames_camcap
table_complete_names[,5] <- tumour_filenames_icgc
table_complete_names[,6] <- germline_filenames_icgc

colnames(table_complete_names) <- c("Patients", "TumourCC", "BenignCC", "GermlineCC", "TumourICGC", "GermlineICGC")

table_complete_names[table_complete_names=="NONE"] <- NA
#Delete the patients for which there's no benign -> is.na(BenignCC)    (Arraytolookat$Trio==No)
table_complete_names <- table_complete_names[-which(is.na(table_complete_names[,3])),]

# Delete those patient s that don't have germline files nor in CC nor in ICGC
table_complete_names <- table_complete_names[-which(is.na(table_complete_names[,4]) & is.na(table_complete_names[,6])),]

write.table(table_complete_names, "table_complete_names.csv", col.names=TRUE, row.names=FALSE, sep=",")
#----------------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------------#
#Merge all the logR and BAF in one file for each tumour, benign and germline.

#I copied all TRIOS files in the folder DataFiles where I also wrote the LP600...files

#The data comes from 3 differnt manifest :25M-8v1_1_B, 2.5-8v1_A and 2.5-8v1_C.
# I use two examples of these 3 manifests (1 each) and select the SNPs that exists in both


cc <- rep("NULL",6) #7 = numebr of columns of the files
cc[c(2,3,4)]<-NA

file_temp <- read.table("DataFiles/LP6005185-DNA_A05.txt",skip=10,header=TRUE,sep=" ",colClasses=cc) #HumanOmni2.5-8v1_A.bpm
file_temp2 <- rbind(file_temp,read.table("DataFiles/Neal193.FinalReport.txt",skip=10,header=TRUE,sep="\t",colClasses=cc)) #HumanOmni2.5-8v1_C.bpm

file_temp3 <- file_temp2[duplicated(file_temp2[,c(1,2)])==TRUE,]

# not_dup <- tail(which(duplicated(file_temp2)==FALSE))[-c(1,2)]
# 
# file_temp[not_dup,]
# SNP.Name Chr Position
# rs35437602  XY  1303007
# rs311094  XY  2604460
# rs311093  XY  2603913
# rs7876632  XY   722984
# 
# > file_temp2[file_temp2$SNP.Name=="rs35437602",]
# SNP.Name Chr Position
# rs35437602   0        0
# > file_temp2[file_temp2$SNP.Name=="rs311094",]
# SNP.Name Chr Position
# rs311094  XY  2654460
# > file_temp2[file_temp2$SNP.Name=="rs311093",]
# SNP.Name Chr Position
# rs311093  XY  2653913
# > file_temp2[file_temp2$SNP.Name=="rs7876632",]
# SNP.Name Chr Position
# rs7876632  XY   772984

file_temp4 <- read.table(paste("DataFiles/LP6005858-DNA_D01.txt",sep=""),skip=10,header=TRUE,sep=" ",colClasses=cc) #HumanOmni25M-8v1-1_B.bpm

file_temp5 <- file_temp4[which(file_temp4$SNP.Name%in%file_temp3$SNP.Name & file_temp4$Chr%in%file_temp3$Chr),] # for some reasons duplicated doesn't recognize all of them.

#Delete all the chr=0, MT, XY

all_snps <- file_temp5[which(file_temp5$Chr!=0),]
all_snps <- all_snps[which(all_snps$Chr!="MT"),]
all_snps <- all_snps[which(all_snps$Chr!="XY"),]

#----------------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------------#
# Order the snps

all_snps <- all_snps[order(all_snps[,2],all_snps[,3]),]

#----------------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------------#
#Read files and create logr and baf files

source("functions_list.R")

tumour_logr <- create_logr_baf_files(table_complete_names[,2], table_complete_names[,5],all_snps,"logr")
rownames(tumour_logr) <- tumour_logr[,1]
tumour_logr <- tumour_logr[,-1]
write.table(tumour_logr, "tumour_logr.txt", quote=FALSE, row.names=TRUE)
save(tumour_logr, file="tumour_logr.RData")

benign_logr <- create_logr_baf_files(table_complete_names[,3], NULL, all_snps,"logr")
rownames(benign_logr) <- benign_logr[,1]
benign_logr <- benign_logr[,-1]
write.table(benign_logr, "benign_logr.txt", quote=FALSE, row.names=TRUE)
save(benign_logr, file="benign_logr.RData")

germline_logr <- create_logr_baf_files(table_complete_names[,4], table_complete_names[,6],all_snps,"logr")
rownames(germline_logr) <- germline_logr[,1]
germline_logr <- germline_logr[,-1]
write.table(germline_logr, "germline_logr.txt", quote=FALSE, row.names=TRUE)
save(germline_logr, file="germline_logr.RData")

#-------------

tumour_baf <- create_logr_baf_files(table_complete_names[,2], table_complete_names[,5],all_snps,"baf")
rownames(tumour_baf) <- tumour_baf[,1]
tumour_baf <- tumour_baf[,-1]
write.table(tumour_baf, "tumour_baf.txt", quote=FALSE, row.names=FALSE)
save(tumour_baf, file="tumour_baf.RData")

benign_baf <- create_logr_baf_files(table_complete_names[,3], NULL, all_snps,"baf")
rownames(benign_baf) <- benign_baf[,1]
benign_baf <- benign_baf[,-1]
write.table(benign_baf, "benign_baf.txt", quote=FALSE, row.names=FALSE)
# benign_baf <- read.table("benign_baf.txt", header=TRUE, sep=" ")
save(benign_baf, file="benign_baf.RData")

germline_baf <- create_logr_baf_files(table_complete_names[,4], table_complete_names[,6],all_snps,"baf")
rownames(germline_baf) <- germline_baf[,1]
germline_baf <- germline_baf[,-1]
write.table(germline_baf, "germline_baf.txt", quote=FALSE, row.names=FALSE)
germline_baf <- read.table("germline_baf.txt", header=TRUE, sep=" ")
save(germline_baf, file="germline_baf.RData")

#----------------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------------#
#Delete NAs

source("functions_list.R")

# - They are not
clean_benign_logr <- discard_nas(benign_logr,3)
clean_germline_logr <- discard_nas(germline_logr,3)
clean_tumour_logr <- discard_nas(tumour_logr,3)

clean_benign_baf <- discard_nas(benign_baf,3)
clean_germline_baf <- discard_nas(germline_baf,3)
clean_tumour_baf <- discard_nas(tumour_baf,3)

chrs <- c(1:22,"X","Y")#,"XY","MT")

source("functions_list.R")

benign50 <- add_50_before_and_after(clean_benign_logr,chrs)
germline50 <- add_50_before_and_after(clean_germline_logr,chrs)
tumour50 <- add_50_before_and_after(clean_tumour_logr,chrs)

span=0.3
correct_single_benign <- apply_loess_to_single_chr(benign50,chrs,span)
correct_single_germline <- apply_loess_to_single_chr(germline50,chrs,span)
correct_single_tumour <- apply_loess_to_single_chr(tumour50,chrs,span)

benign <- lapply(correct_single_benign, function(i) i[-c(1:50,((nrow(i)-49):nrow(i))),])
germline <- lapply(correct_single_germline, function(i) i[-c(1:50,((nrow(i)-49):nrow(i))),])
tumour <- lapply(correct_single_tumour, function(i) i[-c(1:50,((nrow(i)-49):nrow(i))),])

wind <- 5
runmed_benign <- apply_runmed(benign,wind)
runmed_germline <- apply_runmed(germline,wind)
runmed_tumour <- apply_runmed(tumour,wind)

#---------------------------------------------------------------------------------------------------------------------------#
clean_benign_baf_chr <- lapply(chrs,function(i) as.matrix(clean_benign_baf[clean_benign_baf[,1]==i,-c(1,2)]))
clean_germline_baf_chr <- lapply(chrs,function(i) as.matrix(clean_germline_baf[clean_germline_baf[,1]==i,-c(1,2)]))
clean_tumour_baf_chr <- lapply(chrs,function(i) as.matrix(clean_tumour_baf[clean_tumour_baf[,1]==i,-c(1,2)]))

wind <- 5
runmed_benign_baf <- apply_runmed(clean_benign_baf_chr,wind)
runmed_germline_baf <- apply_runmed(clean_germline_baf_chr,wind)
runmed_tumour_baf <- apply_runmed(clean_tumour_baf_chr,wind)

save(runmed_tumour,file="runmed_tumour.RData")
save(runmed_benign, file="runmed_benign.RData")
save(runmed_germline,file="runmed_germline.RData")

load("runmed_benign.RData")
load("runmed_germline.RData")
load("runmed_tumour.RData")

save(runmed_tumour_baf,file="runmed_tumour_baf.RData")
save(runmed_benign_baf, file="runmed_benign_baf.RData")
save(runmed_germline_baf,file="runmed_germline_baf.RData")

load("runmed_benign_baf.RData")
load("runmed_germline_baf.RData")
load("runmed_tumour_baf.RData")


#---------------------------------------------------------------------------------------------------------------------------#
#  Plot the values way higher than the rest and see if they are close together - enough to be a CNV
#---------------------------------------------------------------------------------------------------------------------------#
load("benign_logr.RData")
BG <- mapply('-',runmed_benign, runmed_germline, SIMPLIFY=FALSE)
BGbaf <- mapply('-',runmed_benign_baf, runmed_germline_baf, SIMPLIFY=FALSE)

chrs <- c(1:22,"X","Y")
positions <- lapply(chrs,function(i) benign_logr[benign_logr$Chr==i,2])
snps <- lapply(chrs,function(i) rownames(benign_logr[benign_logr$Chr==i,]))

# All those highr thnan 0.5 and their neighbourhood

indeces_higher <- vector("list",62)
indeces_higher <- vector("list", length(BG))
for(i in 1:length(BG)){
  indeces_higher[[i]] <- vector("list", ncol(BG[[i]]))
  for(j in 1:ncol(BG[[i]])){
    indeces_higher[[i]][[j]] <- which(BG[[i]][,j]>=0.5)
  }
}

neighb_indeces_higher <- vector("list", length(indeces_higher))
for(i in 1:length(indeces_higher)){
  neighb_indeces_higher[[i]] <- vector("list", length(indeces_higher[[i]]))
  for(j in 1:length(indeces_higher[[i]])){
    neighb <- NULL
    if(length(indeces_higher[[i]][[j]])!=0){
      for(k in 1:length(indeces_higher[[i]][[j]])){
        neighb <- c(neighb, rep((indeces_higher[[i]][[j]][k]-10):(indeces_higher[[i]][[j]][k]+10)))
        neighb2 <- neighb[which(neighb>0 & neighb<=length(BG[[i]][,j]))]
      }
      neighb_indeces_higher[[i]][[j]] <- unique(neighb2)
    }
  }
}


dir.create(file.path(getwd(),"/Plots/plot_original_05"), showWarnings = FALSE)

for(i in 1:length(neighb_indeces_higher)){
  print(i)
  for(j in 1:length(neighb_indeces_higher[[i]])){
    print(j)
    if(length(neighb_indeces_higher[[i]][[j]])!=0){
      tiff(file=paste(getwd(), "/Plots/plot_original_05/Patient_",j,"_chr_",i,".tiff",sep=""),width = 3000, height = 1500)
      par(mfrow=c(2,1))
      plot(positions[[i]][neighb_indeces_higher[[i]][[j]]],BG[[i]][neighb_indeces_higher[[i]][[j]],j],
           main=paste("LOGR - Patient_",j,"_chr_",i,sep=""), pch=20,ylim=c(-1,1),
           xlab="positions",ylab="LogR",cex.main=3,cex.lab=2,cex.axis=2)
      points(positions[[i]][indeces_higher[[i]][[j]]],BG[[i]][indeces_higher[[i]][[j]],j],pch=20,ylim=c(-1,1),col="red")
      
      plot(positions[[i]][neighb_indeces_higher[[i]][[j]]],BGbaf[[i]][neighb_indeces_higher[[i]][[j]],j],
           main=paste("BAF - Patient_",j,"_chr_",i,sep=""), pch=20,ylim=c(-1,1),
           xlab="positions",ylab="LogR",cex.main=3,cex.lab=2,cex.axis=2)
      points(positions[[i]][indeces_higher[[i]][[j]]],BGbaf[[i]][indeces_higher[[i]][[j]],j],pch=20,ylim=c(-1,1),col="red")
      dev.off()
    }
  }
}


#---------------------------------------------------------------------------------------------------------------------------#
#  Create file for PennCNV and DNAcopy
#---------------------------------------------------------------------------------------------------------------------------#

source("functions_list.R")

dir.create(file.path(getwd(),"Input_penncnv_dnacopy"), showWarnings = FALSE)

#Benign-Germline
create_input_files(benign_logr, runmed_benign, runmed_germline, runmed_benign_baf, runmed_germline_baf, "BG")
#Tumour-Benign
create_input_files(benign_logr, runmed_tumour, runmed_benign, runmed_tumour_baf, runmed_benign_baf, "TB")
#Tumour-Germline
create_input_files(benign_logr, runmed_tumour, runmed_germline, runmed_tumour_baf, runmed_germline_baf, "TG")
