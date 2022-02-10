library(RODBC)
library(dplyr)
library(plyr)
library(asreml)
library(tidyr)
library(reshape2)
library("purrr")
library("brms")
library("nadiv")

###########################################################################
### reading pedigree etc from database ####################################
###########################################################################
#### when doing on server this will have to read in from the server #######

db<-"C:\\Users\\s1881212\\Documents\\Deer_database_08_2021\\RedDeer1.99.accdb" #open connection
con<-odbcConnectAccess2007(db)
birth_wt<-sqlFetch(con, "sys_BirthWt") 
# getting trait 
life<-sqlFetch(con, "tbllife") %>% select(Code, BirthYear, MumCode, Sex) 
# getting table life with birth yr etc
mum_stat<-sqlFetch(con, "sys_HindStatusAtConception") %>% select(-CalfBirthYear)%>%dplyr::rename(MumCode=Mum, Code=Calf)
# also getting the status of the hind when calf born
#remember need the R-friendly button on!!###
ped<-sqlFetch(con, "sys_Pedigree")%>%select(Code, MumCode,Sire)
ped$Code <- as.factor(ped$Code)
ped$MumCode <- as.factor(ped$MumCode)
ped$Sire <- as.factor(ped$Sire)

prep<-prepPed(ped)

ped_matrix <- as.matrix(nadiv::makeA(prep))


odbcClose(con) #close connection

################################################
### read in ROH file to work out LG FROH ###
#################################################
setwd("H:/")

###need to make this with updated map positions
FROH<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom", header=T, stringsAsFactors = F)%>%
  select(IID,CHR,KB) %>% dplyr::rename(Code=IID)

KB_perLG<-FROH%>%dplyr::group_by(Code, CHR)%>%
  dplyr::summarise(KB_chr=sum(KB))%>%
  ungroup %>% 
  complete(Code, CHR, fill = list(KB_chr = 0)) #completes for all chromosomes 


deermap <- read.table("PhD_3rdYR/Data_files/TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T, stringsAsFactors = F)
deersum <- deermap %>% 
  dplyr::group_by(CEL.LG) %>%
  dplyr::summarise(Max_EstMb = max(Estimated.Mb.Position))%>%
  dplyr::mutate(chr_length_KB=Max_EstMb/1000)%>%
  dplyr::rename(CHR=CEL.LG)


froh_per_chr<-join(KB_perLG,deersum)%>% mutate(chr_froh=KB_chr/chr_length_KB)%>%
  select(-Max_EstMb)%>% dcast(Code~CHR)


colnames(froh_per_chr) <- c("Code", paste0("FROH_chr", 1:33))

Birth_wt_df<-birth_wt%>%join(life)%>%join(mum_stat)%>%join(froh_per_chr)%>%na.omit()

#need to set factors for categorical variables#
Birth_wt_df$Sex <- factor(Birth_wt_df$Sex)
Birth_wt_df$Code <- factor(Birth_wt_df$Code)
Birth_wt_df$MotherStatus <- factor(Birth_wt_df$MotherStatus)
Birth_wt_df$MumCode <- factor(Birth_wt_df$MumCode)


#save(Birth_wt_df, file="PhD_3rdYR/Data_files/R_data_for_brms/birthwt_df.RData")


#####################################################
###### Running model in stan using brms paclage #####
#####################################################

### following code from https://juliengamartin.github.io/wam_tuto/brms-1.html to fit animal model in brms
# Another great paper with stuff for running animal models with repeated measures in MCMCglmm and brms:
# https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2656.2009.01639.x

Birth_wt_df<-read.table("PhD_3rdYR/Data_files/Files_for_brms/Birth_weight_df.txt", header=T, stringsAsFactors = F,row.names=NULL)


m1 <- brm(CaptureWt~ 1+(1|mm(FROH_chr1,FROH_chr2,FROH_chr3,FROH_chr4,FROH_chr5,FROH_chr6,FROH_chr7,FROH_chr8,
                            FROH_chr9,FROH_chr10,FROH_chr11,FROH_chr12,FROH_chr13,FROH_chr14,FROH_chr15,FROH_chr16,
                            FROH_chr17,FROH_chr18,FROH_chr19,FROH_chr20,FROH_chr21,FROH_chr22,FROH_chr23,FROH_chr24,
                            FROH_chr25,FROH_chr26,FROH_chr27,FROH_chr28,FROH_chr29,FROH_chr30,FROH_chr31,FROH_chr32,
                            FROH_chr33))+
          (1 | gr(Code, cov = ped_matrix))+ #adding in pedigree
          Sex + AgeHrs + MotherStatus + #fixed effects
          (1|BirthYear) +(1|MumCode), ##random effects
          data=Birth_wt_df,
          data2 = list(ped_matrix = ped_matrix), #data for ped info
          cores = 4, iter=1000) #default iterations is 2000, i have 4 cores on laptop, chains=4 default

#Number of cores (defaults to 1). On non-Windows systems, this argument can
#be set globally via the mc.cores option.
#When using within chain parallelization it is  advisable to use just as many threads in total as you have CPU cores


#########################################################################################################
####### SAVING DATA FILES CREATED TO RUN MODEL ON SERVER ################################################
#########################################################################################################


#write.table(Birth_wt_df,
 #           file = "PhD_3rdYR/Data_files/Files_for_brms/Birth_weight_df.txt",
  #          col.names = TRUE, row.names = FALSE, sep = "\t")

#write.table(ped,
#            file = "PhD_3rdYR/Data_files/Files_for_brms/Pedigree_df_for_brms.txt",
#            row.names = F, quote = F, sep = "\t")


