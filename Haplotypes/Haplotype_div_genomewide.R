### For all 33 chr finding haplotype diveristy for 10 SNP windows ###


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
      
     
      Haplotype_diversity <- list()
      
      for (i in seq(from = 1, to = ncol(phased_lab), by = 10)) {
        until <- ifelse((i+9) > ncol(phased_lab),  ncol(phased_lab), i+9) #10 SNP windows 
        if(until==i){break} #if window is only 1 snp long 
        window_temp <- phased_lab[, c(i:(until))]
        if (ncol(window_temp)<10){break}
        window_unite <- window_temp %>% 
          unite("string_haplo", na.rm = TRUE, remove = FALSE) %>% 
          filter_all(all_vars(!grepl("9", .)))%>% #removing the non-called 9s
          group_by(string_haplo)%>% #grouping by unique haplotype
          tally() #counting number of unique haplotypes
        window_unite$ns<-NA
        for (j in 1:nrow(window_unite)){
          window_unite$ns[j]<-window_unite$n[j]*((window_unite$n[j])-1)
          
        }
        
        tops<-sum(window_unite$ns)
        bottom<-sum(window_unite$n)*(sum(window_unite$n)-1)
        D<-1-(tops/bottom)
        Haplotype_diversity <- c(Haplotype_diversity, D)
      }

      window_haplotype_div<-as.data.frame(unlist(Haplotype_diversity))
      window_haplotype_div$window_number<-(1:nrow(window_haplotype_div))
      names(window_haplotype_div)[1]<-"haplotype_div"
      window_haplotype_div$CHR<-k
      save(window_haplotype_div, Haplotype_diversity, phased_lab, file = paste0("Haplotype_div_chr", k, ".RData"))
      
      
}





Genomewide_hap_diversity <- NULL

for(z in 1:33){
  
  load(paste0("Haplotype_div_chr", z, ".RData"))
  
  Genomewide_hap_diversity <- rbind(Genomewide_hap_diversity, window_haplotype_div)
  
  
}


write.table(Genomewide_hap_diversity,
                   file = "Genomewide_hap_div_tbl.txt",
                   row.names = F, quote = F, sep = "\t")
