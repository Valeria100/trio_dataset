#---------------------------------------------------------------------------------------------------------------------------#
#  Read and create a table for the penncnv files
#---------------------------------------------------------------------------------------------------------------------------#

source("functions_list.R")

file_path <- paste(getwd(),"/output_penncnv_files",sep="")
output_penncnv_bg <- read_penncnv(file_path,"BGall.rawcnv", "BG") 
output_penncnv_tb <- read_penncnv(file_path,"TBall.rawcnv", "TB") 
output_penncnv_tg <- read_penncnv(file_path,"TGall.rawcnv", "TG") 


write.table(output_penncnv_bg, file=paste(file_path,"/output_penncnv_bg.txt",sep=""), quote=FALSE, row.names=FALSE)
write.table(output_penncnv_tb, file=paste(file_path,"/output_penncnv_tb.txt",sep=""), quote=FALSE, row.names=FALSE)
write.table(output_penncnv_tg, file=paste(file_path,"/output_penncnv_tg.txt",sep=""), quote=FALSE, row.names=FALSE)


#---------------------------------------------------------------------------------------------------------------------------#
#  Read and plot the BG files
#---------------------------------------------------------------------------------------------------------------------------#

#ordered by chr
BG <- mapply('-',runmed_benign, runmed_germline, SIMPLIFY=FALSE)
BGbaf <- mapply('-',runmed_benign_baf, runmed_germline_baf, SIMPLIFY=FALSE)

chrs <- c(1:22,"X","Y")
positions <- lapply(chrs,function(i) benign_logr[benign_logr$Chr==i,2])
snps <- lapply(chrs,function(i) rownames(benign_logr[benign_logr$Chr==i,]))

#ordered by patient
penncnv_bg <- lapply(1:62, function(i) output_penncnv_bg[output_penncnv_bg$patient==i,])

penncnv_bg_chrs1 <- lapply(1:62, function(i){lapply(chrs[-c(23,24)], function(j) penncnv_bg[[i]][penncnv_bg[[i]]$chr==j,])})

penncnv_bg_chrs <- lapply(penncnv_bg_chrs1, function(i){lapply(i, function(j)  j[order(j$start_pos),] )})


indeces_penncnv <- vector("list",62)
# i-patients j-chrs
for(i in 1:62){
  indeces_penncnv[[i]] <- vector("list",22)
  for(j in 1:22){#ignore X and Y - penncnv doesn't use them
    for(k in 1:length(penncnv_bg_chrs[[i]][[j]]))
    indeces_penncnv[[i]][[j]] <- c(indeces_penncnv[[i]][[j]],which(positions[[j]]>=penncnv_bg_chrs[[i]][[j]][k,]$start_pos & positions[[j]]<=penncnv_bg_chrs[[i]][[j]][k,]$end_pos))
  } 
}

#---------------------------------------------------------------------------------------------------------------------------#

#Plot BG logr and baf
for(i in 1:62){
  for(j in 1:22){
    if(length(indeces_penncnv[[i]][[j]])!=0){
      tiff(file=paste("~/FinalCodeCNVs/Plots/only_bg/Patient_",i,"_chr_",j,".tiff",sep=""),width = 3000, height = 1500)
      par(mfrow=c(2,1))
      
      plot(positions[[j]],BG[[j]][,i], main=paste("LOGR - Patient_",i,"_chr_",j,sep=""), pch=20,ylim=c(-1,1),
           xlab="positions",ylab="LogR",cex.main=3,cex.lab=2,cex.axis=2)
      points(positions[[j]][indeces_penncnv[[i]][[j]]],BG[[j]][indeces_penncnv[[i]][[j]],i], pch=20, col="red",ylim=c(-1,1))
      
      plot(positions[[j]],BGbaf[[j]][,i], main=paste("BAF - Patient_",i,"_chr_",j,sep=""), pch=20,ylim=c(-1,1),
           xlab="positions",ylab="BAF",cex.main=3,cex.lab=2,cex.axis=2)
      points(positions[[j]][indeces_penncnv[[i]][[j]]],BGbaf[[j]][indeces_penncnv[[i]][[j]],i], pch=20, col="red",ylim=c(-1,1))
      dev.off()
    }
  }  
}

#---------------------------------------------------------------------------------------------------------------------------#

#Plot BG and compare it with Benign and Germline alone
for(i in 1:62){
  for(j in 1:22){
    if(length(indeces_penncnv[[i]][[j]])!=0){
      tiff(file=paste("~/FinalCodeCNVs/Plots/bg_with_single_signals/Patient_",i,"_chr_",j,".tiff",sep=""),width = 3000, height = 1500)
      
      par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
      plot(positions[[j]],BG[[j]][,i], main="LOGR BG", pch=20,ylim=c(-1,1),
           xlab="Positions",ylab="LogR",cex.main=3,cex.lab=2,cex.axis=2)
      points(positions[[j]][indeces_penncnv[[i]][[j]]],BG[[j]][indeces_penncnv[[i]][[j]],i], pch=20, col="red",ylim=c(-1,1))
      
      plot(positions[[j]],BGbaf[[j]][,i], main="BAF BG", pch=20,ylim=c(-1,1),
           xlab="",ylab="BAF",cex.main=2,cex.lab=2,cex.axis=2)
      points(positions[[j]][indeces_penncnv[[i]][[j]]],BGbaf[[j]][indeces_penncnv[[i]][[j]],i], pch=20, col="red",ylim=c(-1,1))
      
      plot(positions[[j]],runmed_benign[[j]][,i], main="BENIGN", pch=20,ylim=c(-1,1),
           xlab="",ylab="LogR",cex.main=2,cex.lab=2,cex.axis=2)
      points(positions[[j]][indeces_penncnv[[i]][[j]]],runmed_benign[[j]][indeces_penncnv[[i]][[j]],i], pch=20, col="red",ylim=c(-1,1))
      
      plot(positions[[j]],runmed_germline[[j]][,i], main="GERMLINE", pch=20,ylim=c(-1,1),
           xlab="",ylab="LogR",cex.main=2,cex.lab=2,cex.axis=2)
      points(positions[[j]][indeces_penncnv[[i]][[j]]],runmed_germline[[j]][indeces_penncnv[[i]][[j]],i], pch=20, col="red",ylim=c(-1,1))
      
      mtext(paste("Patient",i," chr ",j,sep=""), outer = TRUE, cex = 3)
      dev.off()
    }
  }  
}

#---------------------------------------------------------------------------------------------------------------------------#
# Plot BG penncnv and the neighbourhood (10 before and after)

neighb_indeces_penncnv <- vector("list",62)
# i-patients j-chrs
for(i in 1:62){
  neighb_indeces_penncnv[[i]] <- vector("list",22)
  for(j in 1:22){#ignore Xand Y - penncnv doesn't use them
#     print(paste(i,j, sep=" "))
    for(k in 1:length(penncnv_bg_chrs[[i]][[j]])){
      indeces <- which(positions[[j]]>=penncnv_bg_chrs[[i]][[j]][k,]$start_pos & positions[[j]]<=penncnv_bg_chrs[[i]][[j]][k,]$end_pos)
      
      if(length(indeces)!=0){
        Start <- ifelse(indeces[1]>=20,indeces[1]-20,1)
        
        ifelse(indeces[1]==1,Start_segment <- 1,Start_segment <- rep(Start:(indeces[1]-1))) 
        
        End <- ifelse(indeces[length(indeces)]<=(nrow(BG[[j]])-20),indeces[length(indeces)]+20,nrow(BG[[j]]))
        ifelse(indeces[length(indeces)]==nrow(BG[[j]]),End_segment <- nrow(BG[[j]]),End_segment <- rep((indeces[length(indeces)]+1):End)) 
          
        new_indeces <- unique(c(Start_segment,indeces,End_segment))
        neighb_indeces_penncnv[[i]][[j]] <- c(neighb_indeces_penncnv[[i]][[j]],new_indeces)
      }
      
    }
  } 
}

#---------------------------------------------------------------------------------------------------------------------------#

#Plot BG logr and baf with neighbours
for(i in 1:62){
  for(j in 1:22){
    
    ind <- neighb_indeces_penncnv[[i]][[j]]
    
    if(length(indeces_penncnv[[i]][[j]])!=0){
      tiff(file=paste("~/FinalCodeCNVs/Plots/plot_neighbourhood/Patient_",i,"_chr_",j,".tiff",sep=""),
           width = 3000, height = 1500)
      par(mfrow=c(2,1))
      
      plot(positions[[j]][ind],BG[[j]][ind,i], main=paste("LOGR - Patient_",i,"_chr_",j,sep=""), pch=20,ylim=c(-1,1),
           xlab="positions",ylab="LogR",cex.main=3,cex.lab=2,cex.axis=2)
      points(positions[[j]][indeces_penncnv[[i]][[j]]],BG[[j]][indeces_penncnv[[i]][[j]],i], pch=20, col="red",ylim=c(-1,1))
      
      plot(positions[[j]][ind],BGbaf[[j]][ind,i], main=paste("BAF - Patient_",i,"_chr_",j,sep=""), pch=20,ylim=c(-1,1),
           xlab="positions",ylab="BAF",cex.main=3,cex.lab=2,cex.axis=2)
      points(positions[[j]][indeces_penncnv[[i]][[j]]],BGbaf[[j]][indeces_penncnv[[i]][[j]],i], pch=20, col="red",ylim=c(-1,1))
      dev.off()
    }
  }  
}

#---------------------------------------------------------------------------------------------------------------------------#

#Plot BG logr and baf with neighbours and original signal higher than 0.5
dir.create(paste(getwd(),"/Plots/plot_neighbourhood_pcnv_05",sep=""))

for(i in 1:62){
  for(j in 1:22){
    
    ind <- neighb_indeces_penncnv[[i]][[j]]
    
    if(length(indeces_penncnv[[i]][[j]])!=0){
      tiff(file=paste("~/FinalCodeCNVs/Plots/plot_neighbourhood_pcnv_05/Patient_",i,"_chr_",j,".tiff",sep=""),
           width = 3000, height = 1500)
      par(mfrow=c(2,1))
      
      plot(positions[[j]][ind],BG[[j]][ind,i], 
           main=paste("LOGR - Patient_",i,"_chr_",j,sep=""), pch=20,ylim=c(-1,1),
           xlab="positions",ylab="LogR",cex.main=3,cex.lab=2,cex.axis=2)
      points(positions[[j]][indeces_penncnv[[i]][[j]]],BG[[j]][indeces_penncnv[[i]][[j]],i], pch=20, col="red",ylim=c(-1,1))
      
      points(positions[[j]][neighb_indeces_higher[[j]][[i]]],BG[[j]][neighb_indeces_higher[[j]][[i]],i], pch=20, col="blue",ylim=c(-1,1))
      points(positions[[j]][indeces_higher[[j]][[i]]],BG[[j]][indeces_higher[[j]][[i]],i],pch=20,ylim=c(-1,1),col="green")
      
      legend("topright",c("neighbourhood_penncnv","penncnv","higher_05","neighbourhood_higher05"),col=c("black","red","green","blue"), cex=2,lty=2)
      
      plot(positions[[j]][ind],BGbaf[[j]][ind,i], 
           main=paste("BAF - Patient_",i,"_chr_",j,sep=""), pch=20,ylim=c(-1,1),
           xlab="positions",ylab="BAF",cex.main=3,cex.lab=2,cex.axis=2)
      points(positions[[j]][indeces_penncnv[[i]][[j]]],BGbaf[[j]][indeces_penncnv[[i]][[j]],i], pch=20, col="red",ylim=c(-1,1))
      
      points(positions[[j]][neighb_indeces_higher[[j]][[i]]],BGbaf[[j]][neighb_indeces_higher[[j]][[i]],i], pch=20, col="blue",ylim=c(-1,1))
      points(positions[[j]][indeces_higher[[j]][[i]]],BGbaf[[j]][indeces_higher[[j]][[i]],i],pch=20,ylim=c(-1,1),col="green")
      
      dev.off()
    }
  }  
}


#---------------------------------------------------------------------------------------------------------------------------#
# Plot BG penncnv and the neighbourhood (10 before and after) 
# SINGLE SEGMENTS!!!

neighb_indeces_penncnv_single <- indeces_penncnv_single <- vector("list",62)
# i-patients j-chrs
for(i in 1:62){
  neighb_indeces_penncnv_single[[i]] <- indeces_penncnv_single[[i]] <- vector("list",22)
  for(j in 1:22){#ignore Xand Y - penncnv doesn't use them
    #     print(paste(i,j, sep=" "))
    neighb_indeces_penncnv_single[[i]][[j]] <- indeces_penncnv_single[[i]][[j]] <- vector("list",nrow(penncnv_bg_chrs[[i]][[j]]))
    for(k in 1:length(penncnv_bg_chrs[[i]][[j]])){
      
       indeces <- which(positions[[j]]>=penncnv_bg_chrs[[i]][[j]][k,]$start_pos & positions[[j]]<=penncnv_bg_chrs[[i]][[j]][k,]$end_pos)
      
      if(length(indeces)!=0){
        Start <- ifelse(indeces[1]>=20,indeces[1]-20,1)
        
        ifelse(indeces[1]==1,Start_segment <- 1,Start_segment <- rep(Start:(indeces[1]-1))) 
        
        End <- ifelse(indeces[length(indeces)]<=(nrow(BG[[j]])-20),indeces[length(indeces)]+20,nrow(BG[[j]]))
        ifelse(indeces[length(indeces)]==nrow(BG[[j]]),End_segment <- nrow(BG[[j]]),End_segment <- rep((indeces[length(indeces)]+1):End)) 
        
        new_indeces <- unique(c(Start_segment,indeces,End_segment))
        neighb_indeces_penncnv_single[[i]][[j]][[k]] <- new_indeces
        indeces_penncnv_single[[i]][[j]][[k]] <- indeces
                                                 
      }
      
    }
  } 
}



#---------------------------------------------------------------------------------------------------------------------------#

#Plot BG logr and baf with neighbours
for(i in 1:62){
  for(j in 1:22){
    
    if(length(indeces_penncnv_single[[i]][[j]])!=0){
      
      for(k in 1:length(indeces_penncnv_single[[i]][[j]])){
      #print(paste(i,j,k, sep=" "))
      
        
        real <- indeces_penncnv_single[[i]][[j]][[k]]
        neigh <- neighb_indeces_penncnv_single[[i]][[j]][[k]]
      if(length(real)!=0){
        tiff(file=paste("~/FinalCodeCNVs/Plots/plot_single_neighbourhood/Patient_",i,"_chr_",j,"_segment_",k,".tiff",sep=""),width = 3000, height = 1500)
        par(mfrow=c(2,1))        
        plot(neigh,BG[[j]][neigh,i], main=paste("LOGR - Patient_",i,"_chr_",j,"_segment_",k,sep=""), pch=20,ylim=c(-1,1),
             xlab="positions",ylab="LogR",cex.main=3,cex.lab=2,cex.axis=2)
        points(real,BG[[j]][real,i], pch=20, col="red",ylim=c(-1,1))
        
        plot(neigh,BGbaf[[j]][neigh,i], main=paste("BAF - Patient_",i,"_chr_",j,sep=""), pch=20,ylim=c(-1,1),
             xlab="positions",ylab="BAF",cex.main=3,cex.lab=2,cex.axis=2)
        points(real,BGbaf[[j]][real,i], pch=20, col="red",ylim=c(-1,1))
        dev.off()
      }
        
      }
    }
  }  
}


#---------------------------------------------------------------------------------------------------------------------------#
#  Read the DNAcopy files
#---------------------------------------------------------------------------------------------------------------------------#
source("functions_list.R")

create_segment(table_complete_names,"BG/segment_object_BG","segment_dnac_BG.RData","segment_gain_loss_BG.RData")
create_segment(table_complete_names,"TG/segment_object_TG","segment_dnac_TG.RData","segment_gain_loss_TG.RData")
create_segment(table_complete_names,"TB/segment_object_TB","segment_dnac_TB.RData","segment_gain_loss_TB.RData")

#I moved the files inside /output_DNAcopy_files/ folder


#---------------------------------------------------------------------------------------------------------------------------#
# Plot the DNAcopy results
load("segment_dnac_BG.RData")
load("segment_gain_loss_BG.RData")
dir.create(paste(getwd(),"/Plots/GainLoss_correspondence_files",sep=""), showWarnings = FALSE)


for(i in 1:length(segment_gain_loss)){
  for(j in 1:length(chrs)){
    gainloss <- segment_gain_loss[[i]][segment_gain_loss[[i]]$chrom==chrs[j],"gain"]+segment_gain_loss[[i]][segment_gain_loss[[i]]$chrom==chrs[j],"loss"]
    indices <- which(gainloss!=0)
    tiff(paste(getwd(),"/Plots/GainLoss_correspondence_files/Patient",i,"Chr",j,".tiff",sep=""),width = 3000, height = 1500)
    par(mfrow=c(2,1))
    plot(positions[[j]],BG[[j]][,i], main="LOGR", pch=20,ylim=c(-1,1))
    points(positions[[j]][indices],BG[[j]][indices,i],col="red")
    plot(positions[[j]],BGbaf[[j]][,i], main="BAF", pch=20,ylim=c(-1,1))    
    points(positions[[j]][indices],BGbaf[[j]][indices,i],col="red")
    
    dev.off()
  }
}




#---------------------------------------------------------------------------------------------------------------------------#
# Create bedtools files
#---------------------------------------------------------------------------------------------------------------------------#

source("functions_list.R")

dir.create(paste(getwd(),"/input_bedtools_files",sep=""), showWarnings = FALSE)
setwd(paste(getwd(),"/input_bedtools_files",sep=""))
dir.create(paste(getwd(),"/BG",sep=""), showWarnings = FALSE)
dir.create(paste(getwd(),"/TG",sep=""), showWarnings = FALSE)
dir.create(paste(getwd(),"/TB",sep=""), showWarnings = FALSE)
setwd("~/FinalCodeCNVs")

uni_pat_bg <- write_pcnv_bed(62,"~/FinalCodeCNVs/output_penncnv_files/output_penncnv_bg.txt",
                             "~/FinalCodeCNVs/input_bedtools_files/BG/bg_pcnv_bed_values_patient_",
                             "~/FinalCodeCNVs/input_bedtools_files/BG/bg_pcnv_bed_patient_")

uni_pat_tg <- write_pcnv_bed(62,"~/FinalCodeCNVs/output_penncnv_files/output_penncnv_tg.txt",
                             "~/FinalCodeCNVs/input_bedtools_files/TG/tg_pcnv_bed_values_patient_",
                             "~/FinalCodeCNVs/input_bedtools_files/TG/tg_pcnv_bed_patient_")

uni_pat_tb <- write_pcnv_bed(62,"~/FinalCodeCNVs/output_penncnv_files/output_penncnv_tb.txt",
                             "~/FinalCodeCNVs/input_bedtools_files/TB/tb_pcnv_bed_values_patient_",
                             "~/FinalCodeCNVs/input_bedtools_files/TB/tb_pcnv_bed_patient_")

save(uni_pat_bg,file="uni_pat_bg.RData")
save(uni_pat_tg,file="uni_pat_tg.RData")
save(uni_pat_tb,file="uni_pat_tb.RData")

load("uni_pat_bg.RData")
load("uni_pat_tg.RData")
load("uni_pat_tb.RData")


#segment_gain_loss

source("functions_list.R")

write_dnac_bed("~/FinalCodeCNVs/output_DNAcopy_files/segment_gain_loss_BG.RData",
               uni_pat_bg,
               "~/FinalCodeCNVs/input_bedtools_files/BG/bg_dnacopy_bed_values_patient_",
               "~/FinalCodeCNVs/input_bedtools_files/BG/bg_dnacopy_bed_patient_")

write_dnac_bed("~/FinalCodeCNVs/output_DNAcopy_files/segment_gain_loss_TG.RData",
               uni_pat_tg,
               "~/FinalCodeCNVs/input_bedtools_files/TG/tg_dnacopy_bed_values_patient_",
               "~/FinalCodeCNVs/input_bedtools_files/TG/tg_dnacopy_bed_patient_")

write_dnac_bed("~/FinalCodeCNVs/output_DNAcopy_files/segment_gain_loss_TB.RData",
               uni_pat_tb,
               "~/FinalCodeCNVs/input_bedtools_files/TB/tb_dnacopy_bed_values_patient_",
               "~/FinalCodeCNVs/input_bedtools_files/TB/tb_dnacopy_bed_patient_")


