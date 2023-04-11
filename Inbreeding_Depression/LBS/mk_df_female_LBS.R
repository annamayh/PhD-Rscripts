library(RODBC)
library(tidyverse)
library(MCMCglmm)
library(MasterBayes)
library(data.table)


db<-"C:\\Users\\s1881212\\Documents\\Deer_database_2022/RedDeer2.05.accdb" #open connection
con<-odbcConnectAccess2007(db)

female_LBS<-sqlFetch(con, "sys_LBS+LRS") #built in query - Calculate LBS and LRS for dead hinds

life<-sqlFetch(con, "tbllife") %>% 
  dplyr::select(Code, BirthYear, MumCode)

odbcClose(con)



setwd("H:/")

FROH<-read.table("PhD_4th_yr/2023_ROH_search/2021_calves_ROH_UpdatedSortedMb_032023.hom", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,CHR,KB) %>% dplyr::rename(Code=IID)%>%filter(nchar(Code)==5)

KB_perLG<-FROH%>%dplyr::group_by(Code, CHR)%>%
  dplyr::summarise(KB_chr=sum(KB))%>%
  ungroup %>% 
  complete(Code, CHR, fill = list(KB_chr = 0)) #completes for all chromosomes 

deermap <- read.csv("PhD_3rdYR/Data_files/Genome_assembly_mCerEla1.1.csv", header = T, stringsAsFactors = F)%>%
  dplyr::filter(!CHR %in% c("All","All_auto","X","unplaced"))

froh_per_chr<-plyr::join(KB_perLG,deermap)%>% mutate(chr_froh=KB_chr/length_Kb)%>%
  dplyr::select(-length_Mb,-length)%>% reshape2::dcast(Code~CHR) %>% mutate(FROHsum = rowSums(.[2:34]))%>%mutate(FROH_sum_div=FROHsum/33)

colnames(froh_per_chr) <- c("Code", paste0("FROH_chr", 1:33),"FROHsum","FROH_sum_div")





female_LBS_df=female_LBS%>%
  left_join(froh_per_chr)%>%left_join(life)%>%
  na.omit()



write.table(female_LBS_df,
            file = "PhD_4th_yr/Inbreeding_depression_models/LBS/female_LBS_df.txt",
            row.names = F, quote = F, sep = ",",na = "NA") #saving tables as txt file 


