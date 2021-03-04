## For all 33 chr finding the haplptype count in a 10 SNP window ##
# Based off haplotype div script#
library(dplyr)
library(tidyverse)
library(janitor)
library(stringdist)

setwd("H:/PHD_2ndYR/Haplotype_diversity")

for (k in 1:33){
      
  
      phased<-read.table(paste0("AlphaPeel_files/all_chr/AlphaPeel_out_chr", k, ".called_phase.0.95"),stringsAsFactors = FALSE)
      #Lucys phased data for chr 15, each row is a SNP, 1317 rows 
      markers<-read.table(paste0("AlphaPeel_files/all_chr/Deer_AlphaPeel_marker_file_chr", k, ".txt"), stringsAsFactors = FALSE)%>%
        select(V2) #reading in marker names used OR V3 = BP 
      names(phased) <- c("IID", markers$V2[-1])
      
      
      
      rum_roh<-read.table("old_files/Froh_2019.txt", header = TRUE) %>%
        select(IID) ## IDs used in ROH analysis 
      rum_roh$in_dataset<-"yes" #adding vector to filter by later
      phased_lab<-plyr::join(phased, rum_roh)%>%filter(in_dataset == "yes")%>% 
        select(-IID, -in_dataset) #filtering out IDs not used in ROH analysis 
      
      
      
      ##### counting number of haplotypes in defined region #####
      ###########################################################
      
     
      Haplotype_count <- list()
      
      for (i in seq(from = 1, to = ncol(phased_lab), by = 10)) {
        until <- ifelse((i+9) > ncol(phased_lab),  ncol(phased_lab), i+9)
        if(until==i){break}
        window_temp <- phased_lab[, c(i:(until))]
        if (ncol(window_temp)<10){break}
        window_unite <- window_temp %>% 
          unite("string_haplo", na.rm = TRUE, remove = FALSE) %>% 
          filter_all(all_vars(!grepl("9", .)))%>%
          group_by(string_haplo)%>%
          tally() %>%
          nrow()
       
        Haplotype_count <- c(Haplotype_count, window_unite)
      }

      window_haplotype_count<-as.data.frame(unlist(Haplotype_count))
      window_haplotype_count$window_number<-(1:nrow(window_haplotype_count))
      names(window_haplotype_count)[1]<-"haplotype_count"
      window_haplotype_count$CHR<-k
      save(window_haplotype_count, Haplotype_count, phased_lab, file = paste0("Haplotype_count_chr", k, ".RData"))
      
      
}





Genomewide_hap_count <- NULL

for(z in 1:33){
  
  load(paste0("Haplotype_count_chr", z, ".RData"))
  
  Genomewide_hap_count <- rbind(Genomewide_hap_count, window_haplotype_count)
  
  
}


write.table(Genomewide_hap_count,
                   file = "Genomewide_hap_COUNT_tbl.txt",
                   row.names = F, quote = F, sep = "\t")
