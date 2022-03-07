library(RODBC)
library(dplyr)
library(plyr)
library(asreml)
library(tidyr)
library(reshape2)
library("purrr")
library("brms")
library("tidybayes")




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

odbcClose(con) #close connection


setwd("H:/")

###need to make this with updated map positions
FROH<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom.indiv", header=T, stringsAsFactors = F)%>%
  select(IID,KB) %>% dplyr::rename(Code=IID)%>%mutate(FROH=KB/2495700)


Birth_wt_df<-birth_wt%>%join(life)%>%join(mum_stat)%>%join(FROH)%>%na.omit()

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

m1 <- brm(CaptureWt~ 1 +
          FROH + Sex + AgeHrs + MotherStatus + #fixed effects
          (1|BirthYear) +(1|MumCode), ##random effects
          data=Birth_wt_df,
          cores = 4, iter=5000) #default iterations is 2000, i have 4 cores on laptop, chains=4 default

#Number of cores (defaults to 1). On non-Windows systems, this argument can
#be set globally via the mc.cores option.
#When using within chain parallelization it is  advisable to use just as many threads in total as you have CPU cores


## now post stuff to play with 

summary(m1)

plot(m1)

get_variables(m1) #shows all the vairaibles in you model


