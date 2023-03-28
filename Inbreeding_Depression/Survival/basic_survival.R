library(RODBC)
library(tidyverse)
library("MCMCglmm")
library(MasterBayes)


###########################################################################
### reading pedigree etc from database ####################################
###########################################################################
# 
# 
db<-"C:\\Users\\s1881212\\Documents\\Deer_database_08_2021\\RedDeer1.99.accdb" #open connection
con<-odbcConnectAccess2007(db)

mum_birthyr<-sqlFetch(con, "tbllife")%>%dplyr::select(Code,BirthYear)%>%dplyr::rename(mum_birthyear=BirthYear,MumCode=Code)
mum_stat<-sqlFetch(con, "sys_HindStatusAtConception") %>% dplyr::select(-CalfBirthYear)%>%dplyr::rename(MumCode=Mum, Code=Calf)
birth_wt<-sqlFetch(con, "sys_BirthWt")%>%dplyr::select(BirthWt,Code)

life<-sqlFetch(con, "tbllife") %>%
  dplyr::select(Code, BirthDay,BirthYear, BirthMonth, BirthYear,DeathDay,DeathMonth,DeathYear, Sex, MumCode, DeathType)%>%
  filter(DeathType!= c("S","A","D"))%>%
  drop_na(BirthYear,DeathYear,DeathMonth,BirthMonth)%>%
  plyr::join(mum_birthyr)%>%mutate(mum_age=BirthYear-mum_birthyear)%>%
  plyr::join(mum_stat)%>%plyr::join(birth_wt)

ped<-sqlFetch(con, "sys_Pedigree")%>%dplyr::select(Code, MumCode,Sire)

odbcClose(con)

########################
### neonatl survival ###
#######################

#ids that die before Oct 1st in their first yr
neo<-life%>%mutate(mum_age_sq=mum_age^2)

neo$neonatal_survival<-NA
neo$neonatal_survival<-if_else(neo$DeathMonth<10 & neo$DeathYear==neo$BirthYear, 0, 1)


########################
### winter survival ####
########################
winter<-neo%>%filter(neonatal_survival==1) #removing 714 ids that died as neonates
#if ids die over winter (Oct 1st - May 1st)
winter$winter_surv<-NA
winter$winter_surv<-if_else((winter$DeathYear==winter$BirthYear|
                               winter$DeathYear==(winter$BirthYear)+1&winter$DeathMonth<5), 0, 1)#dying between Summer of of yr born and end of yr

########################
## yearling survival i.e. up to May 1st of second yr
###########################
year<-winter%>%filter(winter_surv==1) #removing  ids that die over winter

year$yearling_surv<-NA
year$yearling_surv<-if_else(year$DeathYear==(year$BirthYear)+2 & year$DeathMonth<5|
                              year$DeathYear==(year$BirthYear)+1 & year$DeathMonth>5, 0, 1)

##########################
## juvenile survival #####
##########################

juv<-year%>%filter(yearling_surv==1)%>%dplyr::select(-yearling_surv,-winter_surv,-neonatal_survival)## removing 864 ids that died before first yr
juv$juvenile_surv<-1

juv_all<-life%>%plyr::join(juv)%>%mutate_at(vars("juvenile_surv"),~replace_na(.,0))%>%mutate(mum_age_sq=mum_age^2)


#set up df
setwd("H:/")

FROH_id<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom.indiv", header=T, stringsAsFactors = F)%>%
  dplyr::rename(Code=IID)%>%mutate(FROH=KB/2591865)


mum_FROH=FROH_id%>%rename(MumCode=Code, MumFROH=FROH)%>%dplyr::select(MumCode, MumFROH)

survival_df<-juv_all%>%plyr::join(FROH_id)%>%
  plyr::join(mum_FROH)%>%
  dplyr::select(Code, juvenile_surv, Sex, MotherStatus, mum_age, mum_age_sq, FROH, BirthYear, MumCode, MumFROH)%>%
  na.omit() #



k<-100
prior<-list(R=list(V=1,fix=1),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


survival_df$Sex <- factor(survival_df$Sex)
survival_df$Code <- factor(survival_df$Code)
survival_df$MotherStatus <- factor(survival_df$MotherStatus)
survival_df$MumCode <- factor(survival_df$MumCode)
survival_df$BirthYear <- factor(survival_df$BirthYear)
  
  basic_surv_m<-MCMCglmm(juvenile_surv~1 + Sex + MotherStatus + mum_age+
                           mum_age_sq+FROH+MumFROH, 
                      random= ~ BirthYear +MumCode,
                      family="threshold",
                      data=survival_df,
                      prior=prior, 
                      #nitt=100000,burnin=20000
                      ) ##

  
summary(basic_surv_m)  

# Location effects: juvenile_surv ~ 1 + Sex + MotherStatus + mum_age + mum_age_sq + FROH + MumFROH 
# 
# post.mean   l-95% CI   u-95% CI eff.samp  pMCMC    
# (Intercept)               0.040119  -0.772455   0.861504   1022.8  0.938    
# Sex2                     -0.423646  -0.571532  -0.279818    828.4 <0.001 ***
#   MotherStatusNaïve        -0.169785  -0.498500   0.123170   1000.0  0.282    
# MotherStatusSummer yeld   0.030490  -0.233860   0.306173   1000.0  0.796    
# MotherStatusTrue yeld     0.134808  -0.067933   0.326293   1000.0  0.202    
# MotherStatusWinter yeld  -0.286043  -0.582548  -0.020663   1000.0  0.040 *  
#   mum_age                   0.199064   0.056650   0.367773    908.3  0.008 ** 
#   mum_age_sq               -0.013830  -0.021905  -0.005755    888.0 <0.001 ***
#   FROH                     -7.391613 -10.519353  -4.180944   1000.0 <0.001 ***
#   MumFROH                  -1.167029  -6.020535   2.807501   1000.0  0.620    




##### now add in spatial region 



db<-"C:\\Users\\s1881212\\Documents\\Deer_database_08_2021\\RedDeer1.99.accdb" #open connection
con<-odbcConnectAccess2007(db)

census<-sqlFetch(con, "tblCensus")%>%dplyr::select(Date,Code,Northing,Easting)
census_yr=census%>%separate(Date,into=c("Year","month","Day"),sep="-")

odbcClose(con)

LocToReg6 <- function(E, N) {
  ifelse(N < 8019, "SG", #south glen
         ifelse(E < 1361, "LA", # laundry greens
                ifelse(E < 1366,
                       ifelse(N > 8033, "NG", "MG"), # north glen, mid  glen 
                       ifelse(E < 1373 , "IM", "SI")))) }

## grouped per year

ave_N_E_peryear=census_yr%>%dplyr::select(-month, -Day)%>%
  group_by(Code,Year)%>%
  summarise(yr_mean_N=mean(Northing), yr_mean_E=mean(Easting))

ave_N_E_peryear$Reg <- with(ave_N_E_peryear, LocToReg6(yr_mean_E, yr_mean_N))



surv_df_mum_loc=as.data.frame(ave_N_E_peryear%>%rename(MumCode=Code, BirthYear=Year)%>% #so we get the average location of the mum during the birth year of the calf
  right_join(survival_df)%>%na.omit()%>%dplyr::select(-yr_mean_E,-yr_mean_N))
  

prior2<-list(R=list(V=1,fix=1),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G2=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G3=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

surv_df_mum_loc$Sex <- factor(surv_df_mum_loc$Sex)
surv_df_mum_loc$Code <- factor(surv_df_mum_loc$Code)
surv_df_mum_loc$MotherStatus <- factor(surv_df_mum_loc$MotherStatus)
surv_df_mum_loc$MumCode <- factor(surv_df_mum_loc$MumCode)
surv_df_mum_loc$BirthYear <- factor(surv_df_mum_loc$BirthYear)
surv_df_mum_loc$Reg <- factor(surv_df_mum_loc$Reg)

  
basic_surv_loc<-MCMCglmm(juvenile_surv~1 + Sex + MotherStatus + mum_age+
                         mum_age_sq+FROH+MumFROH, 
                       random= ~ BirthYear +MumCode +Reg,
                       family="threshold",
                       data=surv_df_mum_loc,
                       prior=prior2, 
                       #nitt=100000,burnin=20000
) ##


summary(basic_surv_loc)  

# Location effects: juvenile_surv ~ 1 + Sex + MotherStatus + mum_age + mum_age_sq + FROH + MumFROH 
# 
# post.mean  l-95% CI  u-95% CI eff.samp  pMCMC    
# (Intercept)             -0.413905 -1.345159  0.394859   1000.0  0.370    
# Sex2                    -0.411319 -0.546802 -0.267191   1000.0 <0.001 ***
#   MotherStatusNaïve       -0.108623 -0.400096  0.201323   1000.0  0.522    
# MotherStatusSummer yeld  0.021172 -0.245838  0.271046    888.8  0.880    
# MotherStatusTrue yeld    0.114792 -0.063153  0.326877   1000.0  0.250    
# MotherStatusWinter yeld -0.274386 -0.531231  0.016776   1000.0  0.046 *  
#   mum_age                  0.257876  0.101842  0.409514   1000.0 <0.001 ***
#   mum_age_sq              -0.016289 -0.024538 -0.008715   1000.0 <0.001 ***
#   FROH                    -6.224244 -9.765107 -3.198858   1000.0 <0.001 ***
#   MumFROH                 -0.193528 -4.281787  3.783348   1000.0  0.894    
# 
# 


plot(basic_surv_loc)



## try Fgrm instead


Fgrm_full<-read.table("PhD_4th_yr/test_Fgrm_for_loeske/Fgrm_full_subset.ibc", header=T, stringsAsFactors = F)%>%
  rename(Code=IID)%>%dplyr::select(Code, Fhat3)

mum_fgrm=Fgrm_full%>%rename(MumCode=Code, Mumfgrm=Fhat3)

survival_df_fgrm=as.data.frame(Fgrm_full%>%
                                 left_join(survival_df)%>%left_join(mum_fgrm)%>%
                                 na.omit())



survival_df_fgrm$Sex <- factor(survival_df_fgrm$Sex)
survival_df_fgrm$Code <- factor(survival_df_fgrm$Code)
survival_df_fgrm$MotherStatus <- factor(survival_df_fgrm$MotherStatus)
survival_df_fgrm$MumCode <- factor(survival_df_fgrm$MumCode)
survival_df_fgrm$BirthYear <- factor(survival_df_fgrm$BirthYear)
#survival_df_fgrm$Reg <- factor(survival_df_fgrm$Reg)


basic_surv_m_fgrm<-MCMCglmm(juvenile_surv~1 + Sex + MotherStatus + mum_age+
                         mum_age_sq+Fhat3+Mumfgrm, 
                       random= ~ BirthYear +MumCode,
                       family="threshold",
                       data=survival_df_fgrm,
                       prior=prior, 
                       #nitt=100000,burnin=20000
) ##

summary(basic_surv_m_fgrm)


## now again with Region



surv_df_mum_loc_fgrm=as.data.frame(ave_N_E_peryear%>%rename(MumCode=Code, BirthYear=Year)%>% #so we get the average location of the mum during the birth year of the calf
                                right_join(survival_df_fgrm)%>%na.omit()%>%dplyr::select(-yr_mean_E,-yr_mean_N))


surv_df_mum_loc_fgrm$Sex <- factor(surv_df_mum_loc_fgrm$Sex)
surv_df_mum_loc_fgrm$Code <- factor(surv_df_mum_loc_fgrm$Code)
surv_df_mum_loc_fgrm$MotherStatus <- factor(surv_df_mum_loc_fgrm$MotherStatus)
surv_df_mum_loc_fgrm$MumCode <- factor(surv_df_mum_loc_fgrm$MumCode)
surv_df_mum_loc_fgrm$BirthYear <- factor(surv_df_mum_loc_fgrm$BirthYear)
surv_df_mum_loc_fgrm$Reg <- factor(surv_df_mum_loc_fgrm$Reg)


basic_surv_loc_fgrm<-MCMCglmm(juvenile_surv~1 + Sex + MotherStatus + mum_age+
                           mum_age_sq+Fhat3+Mumfgrm, 
                         random= ~ BirthYear +MumCode +Reg,
                         family="threshold",
                         data=surv_df_mum_loc_fgrm,
                         prior=prior2, 
                         #nitt=100000,burnin=20000
) ##


summary(basic_surv_loc_fgrm)  

plot(basic_surv_loc_fgrm)
