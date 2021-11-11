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

  
Focal_LG_FROH<-FROH%>%filter(CHR==k)%>%dplyr::group_by(Code, CHR)%>% #filter for focal chr and grouping by id and chr
  dplyr::summarise(KB_chr=sum(KB))%>% #getting total KB per chr per id
  ungroup %>% 
  complete(Code, CHR, fill = list(KB_chr = 0)) %>% join(deersum)%>% #joining to total length of chr
  mutate(focal_chr_froh=KB_chr/chr_length_KB)%>% # working out froh for focal chr 
  select(Code, focal_chr_froh)#selecting important columns


deersum_minus_focal<-deersum%>%filter(CHR!=k)%>%filter(CHR!=34) #getting length of chr  minus focal chr
rest_genome_length<-sum(deersum_minus_focal$chr_length_KB) #adding all chr together to get length of genome minus focal

Rest_FROH<-FROH%>%filter(CHR!=k)%>%dplyr::group_by(Code)%>%# getting total KB in ROH for rest of genome
  dplyr::summarise(KB_chr=sum(KB))%>%
  ungroup %>% 
  mutate(Rest_chr_froh=KB_chr/rest_genome_length)%>% #Kb in ROH minus focal/ genome length minus focal 
  select(-KB_chr)
  
for_models[[k]]<-join(Focal_LG_FROH,Rest_FROH, type="right")%>%
  mutate_all(~replace(., is.na(.), 0))#replaces NAs with 0 using dplyr

} 

###############################################################################################
######## NOW  SET THINGS UP TO RUN MODEL #####################################################
#############################################################################################

run_model<-function(focal_LG_df){

  Birth_wt_df<-birth_wt%>%join(life)%>%join(mum_stat)%>%join(focal_LG_df)%>%na.omit() #joining all relevant tables for analysis
#setting as factors 
  Birth_wt_df$Sex <- factor(Birth_wt_df$Sex)
  Birth_wt_df$Code <- factor(Birth_wt_df$Code)
  Birth_wt_df$MotherStatus <- factor(Birth_wt_df$MotherStatus)
  Birth_wt_df$MumCode <- factor(Birth_wt_df$MumCode)


  chr_model_with_ped <- asreml(fixed= CaptureWt ~ 1 + Sex + AgeHrs + Rest_chr_froh + focal_chr_froh + MotherStatus, 
                             random= ~ vm(Code,ainv)+BirthYear +MumCode, 
                             residual= ~ idv(units),
                             data=Birth_wt_df)
  chr_model_with_ped
}
  

all_LG_models<-map(for_models,run_model) # running model for all LG 
  
########################################################
### pulling out effect sizes for birth weight per LG  ###
##########################################################
get_sol_focal_LG<-function(model_output){

fixed_eff_table<-as.data.frame(summary(model_output,coef=TRUE)$coef.fixed)
focal_LG_solution<-fixed_eff_table[6,"solution"] # effect size is the 6th row down in the solutions column 
focal_LG_solution
}

effect_sizes<-as.data.frame(map(all_LG_models, get_sol_focal_LG)%>% #getting the effect size for each LG and making table
  unlist())%>%
  tibble::rownames_to_column()%>%
  setNames(c("CHR","effect_size_of_LG"))


size_v_length<-join(effect_sizes,deersum)%>%select(-"Max_EstMb")

library(ggplot2)
ggplot(size_v_length, aes(x=chr_length_KB,y=effect_size_of_LG))+
  geom_point()+
  theme_classic()+
  geom_text(label=size_v_length$CHR, vjust=1.25)+
  geom_smooth(method=lm)
  


