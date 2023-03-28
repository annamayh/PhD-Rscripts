library(RODBC)
library(tidyverse)
library(data.table)



db<-"C:\\Users\\s1881212\\Documents\\Deer_database_2022/RedDeer2.05.accdb" #open connection
con<-odbcConnectAccess2007(db)
birth_wt<-sqlFetch(con, "sys_BirthWt") 
# getting trait 
life<-sqlFetch(con, "tbllife") %>% dplyr::select(Code, BirthYear, MumCode, Sex) 
# getting table life with birth yr etc
mum_stat<-sqlFetch(con, "sys_HindStatusAtConception") #
# also getting the status of the hind when calf born
#remember need the R-friendly button on!!###
ped<-sqlFetch(con, "sys_Pedigree")%>%dplyr::select(Code, MumCode,Sire)

mum_birthyr<-sqlFetch(con, "tbllife")%>%dplyr::select(Code,BirthYear)%>%dplyr::rename(mum_birthyear=BirthYear,MumCode=Code)
odbcClose(con) #close connection



mum_all_stats=mum_stat%>%dplyr::rename(MumCode=Mum, Code=Calf)%>%
  inner_join(mum_birthyr)%>%
  mutate(mum_age=CalfBirthYear-mum_birthyear)%>%
  mutate(mum_age_sq=mum_age^2)


## read in FROH info 

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
  dplyr::select(-length_Mb,-length)%>% reshape2::dcast(Code~CHR) %>% 
  mutate(FROHsum = rowSums(.[2:34]))%>%mutate(FROH_sum_div=FROHsum/33)


colnames(froh_per_chr) <- c("Code", paste0("FROH_chr", 1:33),"FROHsum","FROH_sum_div")

FROH_full<-read.table("PhD_4th_yr/2023_ROH_search/2021_calves_ROH_UpdatedSortedMb_032023.hom.indiv", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,KB) %>% dplyr::rename(Code=IID)%>%mutate(FROH=KB/2591865)%>%filter(nchar(Code)==5)




## combine to get df

birth_wt_df=birth_wt%>%left_join(mum_all_stats)%>%
  left_join(life)%>%
  left_join(froh_per_chr)%>%
  select(-CalfBirthYear)





write.table(birth_wt_df,
            file = "PhD_4th_yr/Inbreeding_depression_models/birth_weight/birth_wt_df2023.txt",
            row.names = F, quote = F, sep = ",",na = "NA") #saving tables as txt file 


write.table(ped,
            file = "PhD_4th_yr/Inbreeding_depression_models/birth_weight/pedigree.txt",
            row.names = F, quote = F, sep = ",",na = "NA")


