library(RODBC)
library(dplyr)
library(plyr)
library(asreml)
library(tidyr)
library(reshape2)

qgm_summary <-function(asr){
  print("CONVERGENCE")
  print(paste0(asr$converge," with logL = ",asr$loglik))
  cat("\n")
  print("VARIANCE COMPONENTS")
  print(summary(asr)$varcomp)
  cat("\n")
  print("FIXED EFFECTS")
  print(as.data.frame(summary(asr,coef=TRUE)$coef.fixed))
  cat("\n")
  print("WALD TESTS")
  print(as.data.frame(wald(asr,denDF="numeric",ssType="conditional")$Wald)[,c(1,2,4,6)])
}

###########################################################################
### reading pedigree etc from database ####################################
###########################################################################

db<-"C:\\Users\\s1881212\\Documents\\Deer_database_08_2021\\RedDeer1.99.accdb" #open connection
con<-odbcConnectAccess2007(db)
birth_wt<-sqlFetch(con, "sys_BirthWt") 
# getting trait 
life<-sqlFetch(con, "tbllife") %>% select(Code, BirthYear, MumCode, Sex) 
# getting table life with birth yr etc
mum_stat<-sqlFetch(con, "sys_HindStatusAtConception") %>% select(-CalfBirthYear)%>%dplyr::rename(MumCode=Mum, Code=Calf)
# also getting the status of the hind when calf born
#remember need the R-friendly button on!!###
ped<-sqlFetch(con, "sys_Pedigree")%>%select(Code, MumCode, Sire)
odbcClose(con) #close connection


################################################
### read in ROH file to work out FROH per chr ###
#################################################
setwd("H:/")

FROH<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_bp_092021.hom", header=T, stringsAsFactors = F)%>%
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

################################################################
##### Make list of chromosomal FROH plus FROH -focal chr ######
###### to use in models ######################################
##############################################################

for_models<-list()
for (k in 1:33){
  
  for_models[[k]]<-froh_per_chr
  for_models[[k]]$froh_minus_focal<-NA
  
  col_end_chunk_no1<-k
  col_start_chunk_no2<-k+2
 
  if(k==1){
  for_models[[k]]$froh_minus_focal<-rowSums(for_models[[k]][(3:34)])# for chr 1 only need one block to sum
  } 
  else if(k==33){
    for_models[[k]]$froh_minus_focal<-rowSums(for_models[[k]][(2:33)])# for chr 33 again only need one block to sum
  } 
  else{
      for_models[[k]]$froh_minus_focal<-rowSums(for_models[[k]][c(2:col_end_chunk_no1,col_start_chunk_no2:34)])
    # for all the rest sum across rows leaving out the focal chromosome 
    }
    
  for_models[[k]]<-mutate(for_models[[k]],froh_minus_focal_div=froh_minus_focal/32 )#Divide by total number of chr in rest of FROH 
  for_models[[k]]<-dplyr::select(for_models[[k]], Code, froh_minus_focal_div, paste0("FROH_chr", k))#just keeping inportant columns 

  }






###############################################################################################
######## NOW  SET THINGS UP TO RUN MODEL #####################################################
#############################################################################################


Birth_wt_df<-birth_wt %>%join(life)%>%join(mum_stat)%>%join(df_14)%>%na.omit()


#ped_check<-join(ped, FROH, by="Code", type="right") ## or could do whole pedigree
#ainv <- ainverse(ped)
#what to do about missing data??


Birth_wt_df$Sex <- factor(Birth_wt_df$Sex)
Birth_wt_df$Code <- factor(Birth_wt_df$Code)
Birth_wt_df$MotherStatus <- factor(Birth_wt_df$MotherStatus)
Birth_wt_df$MumCode <- factor(Birth_wt_df$MumCode)


chr_model <- asreml(fixed= CaptureWt ~ 1 + Sex + AgeHrs + froh_minus_focal_div + FROH_chr1 + MotherStatus, 
                    random= ~ BirthYear +MumCode, 
                    data=Birth_wt_df)
## add in pedigree

qgm_summary(chr_model)
