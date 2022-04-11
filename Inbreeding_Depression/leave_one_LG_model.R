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
life<-sqlFetch(con, "tbllife") %>% dplyr::select(Code, BirthYear, MumCode, Sex) 
# getting table life with birth yr etc
mum_stat<-sqlFetch(con, "sys_HindStatusAtConception") %>% dplyr::select(-CalfBirthYear)%>%dplyr::rename(MumCode=Mum, Code=Calf)
# also getting the status of the hind when calf born
#remember need the R-friendly button on!!###
ped<-sqlFetch(con, "sys_Pedigree")%>%dplyr::select(Code, MumCode, Sire)
odbcClose(con) #close connection

ped[is.na(ped)] <- 0
ainv <- ainverse(ped)

################################################
### read in ROH file to work out LG FROH ###
#################################################
setwd("H:/")

FROH<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,CHR,KB) %>% dplyr::rename(Code=IID)

KB_perLG<-FROH%>%dplyr::group_by(Code, CHR)%>%
  dplyr::summarise(KB_chr=sum(KB))%>%
  ungroup %>% 
  complete(Code, CHR, fill = list(KB_chr = 0)) #completes for all chromosomes 

deermap <- read.csv("PhD_3rdYR/Data_files/Genome_assembly_mCerEla1.1.csv", header = T, stringsAsFactors = F)%>%
  dplyr::filter(!CHR %in% c("All","All_auto","X","unplaced"))%>%select(1:4)

froh_per_chr<-join(KB_perLG,deermap)%>% mutate(chr_froh=KB_chr/length_Kb)%>%
  dplyr::select(-length_Mb,-length)%>% reshape2::dcast(Code~CHR)


colnames(froh_per_chr) <- c("Code", paste0("FROH_chr", 1:33))


########


for_models<-list()
for (k in 1:33){

  
Focal_LG_FROH<-FROH%>%filter(CHR==k)%>%dplyr::group_by(Code, CHR)%>% #filter for focal chr and grouping by id and chr
  dplyr::summarise(KB_chr=sum(KB))%>% #getting total KB per chr per id
  ungroup %>% 
  complete(Code, CHR, fill = list(KB_chr = 0)) %>% join(deermap)%>% #joining to total length of chr
  mutate(focal_chr_froh=KB_chr/length_Kb)%>% # working out froh for focal chr 
  dplyr::select(Code, focal_chr_froh)#selecting important columns


deersum_minus_focal<-deermap%>%filter(CHR!=k)%>%filter(CHR!=34) #getting length of chr  minus focal chr
rest_genome_length<-sum(deersum_minus_focal$length_Kb) #adding all chr together to get length of genome minus focal

Rest_FROH<-FROH%>%filter(CHR!=k)%>%dplyr::group_by(Code)%>%# getting total KB in ROH for rest of genome
  dplyr::summarise(KB_chr=sum(KB))%>%
  ungroup %>% 
  mutate(Rest_chr_froh=KB_chr/rest_genome_length)%>% #Kb in ROH minus focal/ genome length minus focal 
  dplyr::select(-KB_chr)
  
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


  chr_model_with_ped <- asreml(fixed= CaptureWt ~ 1 + Sex + AgeHrs + 
                                 Rest_chr_froh + focal_chr_froh + MotherStatus, 
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
  setNames(c("CHR","effect_size"))


size_v_length<-join(effect_sizes,deersum)%>%select(-"Max_EstMb")

library(ggplot2)
ggplot(size_v_length, aes(x=chr_length_KB,y=effect_size_of_LG))+
  geom_point()+
  theme_classic()+
  geom_text(label=size_v_length$CHR, vjust=1.25)+
  geom_smooth(method=lm)
  

### pulling out standard deviations ####

get_SE_focal_LG<-function(model_output){
  
  fixed_eff_table<-as.data.frame(summary(model_output,coef=TRUE)$coef.fixed)
  focal_LG_solution<-fixed_eff_table[6,"std error"] # effect size is the 6th row down in the solutions column 
  focal_LG_solution
}

SE<-as.data.frame(map(all_LG_models, get_SE_focal_LG)%>% 
                    unlist())%>%
                    tibble::rownames_to_column()%>%
                    setNames(c("CHR","SE"))


############################################
#### forest plot of effect sizes and CI ###
###########################################
effect_SE<-effect_sizes%>%join(SE)%>% 
  mutate(upper_CI=effect_size+(2*SE))%>%mutate(lower_CI=effect_size-(2*SE))#%>%order(effect_SE$chr_length_KB)

effect_SE$overlap_zero<-ifelse(effect_SE$CHR %in% c("1","10","12","14","18","19","22","26","30", "24","25"),"Yes","No")


#effect_SE$CHR<-as.factor(effect_SE$CHR)
effect_SE$CHR<-ordered(as.integer(effect_SE$CHR))

ggplot(data=effect_SE, aes(x=CHR, y=effect_size, ymin=lower_CI, ymax=upper_CI,color=overlap_zero)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Linkage group", y="Effect size on birth weight + CI", title = "Effect size of Linkage group FROH on Birthweight (kg)") +
  theme_bw()+  # use a white background                 
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey50","red"),guide=FALSE)
