library(RODBC)
library(dplyr)
library(plyr)
library(asreml)
library(tidyr)
library(reshape2)
library("purrr")

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

ped[is.na(ped)] <- 0
ainv <- ainverse(ped)

################################################
### read in ROH file to work out LG FROH ###
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


########





for_models<-list()
for (k in 1:33){

  
Focal_LG_FROH<-FROH%>%filter(CHR==k)%>%dplyr::group_by(Code, CHR)%>%
  dplyr::summarise(KB_chr=sum(KB))%>%
  ungroup %>% 
  complete(Code, CHR, fill = list(KB_chr = 0)) %>% join(deersum)%>%
  mutate(focal_chr_froh=KB_chr/chr_length_KB)%>%
  select(Code, focal_chr_froh)


deersum_minus_focal<-deersum%>%filter(CHR!=k)%>%filter(CHR!=34)
rest_genome_length<-sum(deersum_minus_focal$chr_length_KB)

Rest_FROH<-FROH%>%filter(CHR!=k)%>%dplyr::group_by(Code)%>%
  dplyr::summarise(KB_chr=sum(KB))%>%
  ungroup %>% 
  mutate(Rest_chr_froh=KB_chr/rest_genome_length)%>%
  select(-KB_chr)
  
for_models[[k]]<-join(Focal_LG_FROH,Rest_FROH, type="right")%>%
  mutate_all(~replace(., is.na(.), 0))#replaces NAs with ) using dplyr

} 

###############################################################################################
######## NOW  SET THINGS UP TO RUN MODEL #####################################################
#############################################################################################

for (s in 1:33) {
  
  LG<-for_models[[s]]

  Birth_wt_df<-birth_wt%>%join(life)%>%join(mum_stat)%>%join(LG)%>%na.omit()

  Birth_wt_df$Sex <- factor(Birth_wt_df$Sex)
  Birth_wt_df$Code <- factor(Birth_wt_df$Code)
  Birth_wt_df$MotherStatus <- factor(Birth_wt_df$MotherStatus)
  Birth_wt_df$MumCode <- factor(Birth_wt_df$MumCode)


  chr_model_with_ped <- asreml(fixed= CaptureWt ~ 1 + Sex + AgeHrs + Rest_chr_froh + focal_chr_froh + MotherStatus, 
                             random= ~ vm(Code,ainv)+BirthYear +MumCode, 
                             residual= ~ idv(units),
                             data=Birth_wt_df)

  print(paste0("QGM SUMMARY FOR LG", s))
  qgm_summary(chr_model_with_ped)

}


