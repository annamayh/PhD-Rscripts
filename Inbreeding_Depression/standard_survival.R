library(RODBC)
library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)
library("purrr")
library("MCMCglmm")
library("data.table")
library(tibble)
library(ggplot2)
library(MasterBayes)


###########################################################################
### reading pedigree etc from database ####################################
###########################################################################
# 
# 
# db<-"C:\\Users\\s1881212\\Documents\\Deer_database_08_2021\\RedDeer1.99.accdb" #open connection
# con<-odbcConnectAccess2007(db)
# 
# mum_birthyr<-sqlFetch(con, "tbllife")%>%dplyr::select(Code,BirthYear)%>%dplyr::rename(mum_birthyear=BirthYear,MumCode=Code)
# mum_stat<-sqlFetch(con, "sys_HindStatusAtConception") %>% dplyr::select(-CalfBirthYear)%>%dplyr::rename(MumCode=Mum, Code=Calf)
# birth_wt<-sqlFetch(con, "sys_BirthWt")%>%dplyr::select(BirthWt,Code)
# 
# life<-sqlFetch(con, "tbllife") %>% 
#   dplyr::select(Code, BirthDay,BirthYear, BirthMonth, BirthYear,DeathDay,DeathMonth,DeathYear, Sex, MumCode, DeathType)%>%
#   filter(DeathType!= c("S","A","D"))%>%
#   drop_na(BirthYear,DeathYear,DeathMonth,BirthMonth)%>%
#   join(mum_birthyr)%>%mutate(mum_age=BirthYear-mum_birthyear)%>%
#   join(mum_stat)%>%join(birth_wt)
# 
# ped<-sqlFetch(con, "sys_Pedigree")%>%dplyr::select(Code, MumCode,Sire)
# 
# odbcClose(con)
# 
# ########################
# ### neonatl survival ###
# #######################
# 
# #ids that die before Oct 1st in their first yr
# neo<-life%>%mutate(mum_age_sq=mum_age^2)
# 
# neo$neonatal_survival<-NA
# neo$neonatal_survival<-if_else(neo$DeathMonth<10 & neo$DeathYear==neo$BirthYear, 0, 1)
# 
# 
# ########################
# ### winter survival ####
# ########################
# winter<-neo%>%filter(neonatal_survival==1) #removing 714 ids that died as neonates
# #if ids die over winter (Oct 1st - May 1st)                       
# winter$winter_surv<-NA
# winter$winter_surv<-if_else((winter$DeathYear==winter$BirthYear|
#                                winter$DeathYear==(winter$BirthYear)+1&winter$DeathMonth<5), 0, 1)#dying between Summer of of yr born and end of yr
# 
# ########################
# ## yearling survival i.e. up to May 1st of second yr
# ###########################
# year<-winter%>%filter(winter_surv==1) #removing  ids that die over winter
# 
# year$yearling_surv<-NA
# year$yearling_surv<-if_else(year$DeathYear==(year$BirthYear)+2 & year$DeathMonth<5|
#                               year$DeathYear==(year$BirthYear)+1 & year$DeathMonth>5, 0, 1)
# 
# ##########################
# ## juvenile survival #####
# ##########################
# 
# juv<-year%>%filter(yearling_surv==1)%>%dplyr::select(-yearling_surv,-winter_surv,-neonatal_survival)## removing 864 ids that died before first yr
# juv$juvenile_surv<-1
# 
# juv_all<-life%>%join(juv)%>%mutate_at(vars("juvenile_surv"),~replace_na(.,0))%>%mutate(mum_age_sq=mum_age^2)
# 
# 
# #set up df
# setwd("H:/")
# 
# FROH<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom", header=T, stringsAsFactors = F)%>%
#   dplyr::select(IID,CHR,KB) %>% dplyr::rename(Code=IID)
# 
# KB_perLG<-FROH%>%dplyr::group_by(Code, CHR)%>%
#   dplyr::summarise(KB_chr=sum(KB))%>%
#   ungroup %>% 
#   complete(Code, CHR, fill = list(KB_chr = 0)) #completes for all chromosomes 
# 
# deermap <- read.csv("PhD_3rdYR/Data_files/Genome_assembly_mCerEla1.1.csv", header = T, stringsAsFactors = F)%>%
#   dplyr::filter(!CHR %in% c("All","All_auto","X","unplaced"))
# 
# froh_per_chr<-join(KB_perLG,deermap)%>% mutate(chr_froh=KB_chr/length_Kb)%>%
#   dplyr::select(-length_Mb,-length)%>% reshape2::dcast(Code~CHR)
# 
# colnames(froh_per_chr) <- c("Code", paste0("FROH_chr", 1:33))
# ########
# 
# for_models<-list()
# for (k in 1:33){
#  
#   Focal_LG_FROH<-FROH%>%filter(CHR==k)%>%dplyr::group_by(Code, CHR)%>% #filter for focal chr and grouping by id and chr
#     dplyr::summarise(KB_chr=sum(KB))%>% #getting total KB per chr per id
#     ungroup %>% 
#     complete(Code, CHR, fill = list(KB_chr = 0)) %>% join(deermap)%>% #joining to total length of chr
#     mutate(focal_chr_froh=KB_chr/length_Kb)%>% # working out froh for focal chr 
#     dplyr::select(Code, focal_chr_froh)#selecting important columns
#   
#   
#   deersum_minus_focal<-deermap%>%filter(CHR!=k)%>%filter(CHR!=34) #getting length of chr  minus focal chr
#   rest_genome_length<-sum(deersum_minus_focal$length_Kb) #adding all chr together to get length of genome minus focal
#   
#   Rest_FROH<-FROH%>%filter(CHR!=k)%>%dplyr::group_by(Code)%>%# getting total KB in ROH for rest of genome
#     dplyr::summarise(KB_chr=sum(KB))%>%
#     ungroup %>% 
#     mutate(Rest_chr_froh=KB_chr/rest_genome_length)%>% #Kb in ROH minus focal/ genome length minus focal 
#     dplyr::select(-KB_chr)
#   
#   for_models[[k]]<-join(Focal_LG_FROH,Rest_FROH, type="right")%>%
#     mutate_all(~replace(., is.na(.), 0))#replaces NAs with 0 using dplyr
#   
# } 

#save(for_models,juv_all, file="PhD_3rdYR/Model_inputs/df_ALLsurv_standardmodel.RData")




load("PhD_3rdYR/Model_inputs/df_ALLsurv_standardmodel.RData")

k<-100
prior<-list(R=list(V=1,fix=1),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))






run_model<-function(focal_LG_df){
  
  Birth_wt_df<-juv_all%>%join(focal_LG_df)%>%na.omit() #joining all relevant tables for analysis
  #setting as factors 
  Birth_wt_df$Sex <- factor(Birth_wt_df$Sex)
  Birth_wt_df$Code <- factor(Birth_wt_df$Code)
  Birth_wt_df$MotherStatus <- factor(Birth_wt_df$MotherStatus)
  Birth_wt_df$MumCode <- factor(Birth_wt_df$MumCode)
  
  chr_model<-MCMCglmm(juvenile_surv~1 + Sex + MotherStatus + mum_age+mum_age_sq+
                        Rest_chr_froh + focal_chr_froh, 
                   random= ~ BirthYear +MumCode,
                   family="threshold",
                   data=Birth_wt_df,
                   prior=prior, 
                   nitt=100000,burnin=20000) ## will need to run this on server as 33x. 
  chr_model
}


# test<-for_models[[1]]
# test_model<-run_model(test)
# 
# plot(test_model)

all_LG_models<-map(for_models,run_model) # running model for all LG 

save(all_LG_models, file="PhD_3rdYR/Model outputs/Juvenile_survival_0-2/standard_juvenile_noped_noBW2.RData")


test<-all_LG_models[[25]]
plot(test)

########################################################
### pulling out effect sizes for birth weight per LG  ###
##########################################################


get_focal_sol<-function(model){
    extract_focal_chr<-as.data.frame(model$Sol)%>%dplyr::select(focal_chr_froh) ## 

    sol<-apply(extract_focal_chr,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
    CI_upper<-apply(extract_focal_chr,2,quantile,probs = c(0.95)) #gets upper confidence interval for all solutions 
    CI_lower<-apply(extract_focal_chr,2,quantile,probs = c(0.05)) #gets lower CI for all solutions

    row_focal<-matrix(c(sol,CI_upper,CI_lower), 1, 3)
    row_focal

}


focal_sols_all<-map(all_LG_models,get_focal_sol)


focal_sols<-as.data.frame(do.call(rbind, focal_sols_all))%>%#binds all matricies togther in a df
   set_names(c("sols","CI_up","CI_low"))%>% add_column(CHR = 1:33)
  
focal_sols$overlap_zero<-ifelse(focal_sols$CHR %in% c("23","17","12","18"),"Yes","No")


focal_sols$CHR<-as.factor(focal_sols$CHR)

ggplot(data=focal_sols, aes(x=CHR, y=sols, ymin=CI_low, ymax=CI_up,color=overlap_zero)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="Effect on juvenile survival + 95% CI", title = "Effect on juvenile survival + 95% CI") +
  theme_bw()+  # use a white background                 
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey50","red"),guide=FALSE)






