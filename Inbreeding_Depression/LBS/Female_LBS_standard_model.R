library(RODBC)
library(tidyverse)
library(MCMCglmm)
library(MasterBayes)


db<-"C:\\Users\\s1881212\\Documents\\Deer_database_08_2021\\RedDeer1.99.accdb" #open connection
con<-odbcConnectAccess2007(db)

female_LBS<-sqlFetch(con, "sys_LBS+LRS") #built in query - Calculate LBS and LRS for dead hinds

life<-sqlFetch(con, "tbllife") %>% 
  dplyr::select(Code, BirthYear, MumCode)

odbcClose(con)



setwd("H:/")

FROH_full<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom.indiv", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,KB) %>% dplyr::rename(Code=IID)%>%mutate(FROH=KB/2591865)


female_LBS_data=female_LBS%>%
  plyr::join(life)%>%plyr::join(FROH_full)%>%na.omit()



#trace plots are a bit shit
prior.1 = list(R = list(V = diag(2), nu = 0.002,fix = 2),
                   G = list(G1 = list(V = diag(2), nu = 0.002),
                            G2 = list(V = diag(2), nu = 0.002)))

##these are okay?
prior.2 = list(R = list(V = diag(2), nu = 2,fix = 2),
               G = list(G1 = list(V = diag(2), nu = 2),
                        G2 = list(V = diag(2), nu = 2)))

#again okay but not amazing for random 
prior.3 = list(R = list(V = diag(2), nu = 2,fix = 2),
               G = list(G1 = list(V = diag(2), nu = 0.002),
                        G2 = list(V = diag(2), nu = 0.002)))



#what does at.level() do? when i add it in, priors dont work - think it if you have multiple 

female_LBS_model=MCMCglmm(LBS~ trait+trait:FROH-1,
                          #first fixed term = intercept of the Poisson process 
                          #second fixed term = intercept of the zero-inflation  
                          #remaining terms are trait (Poi & Zi) specific contrasts from the intercept.
                          random=~idh(trait):MumCode+idh(trait):BirthYear, 
                          rcov = ~ idh(trait):units,
                          data = female_LBS_data, 
                          family = "zipoisson", 
                          prior = prior.4,
                          #verbose = TRUE, 
                          pr = TRUE, 
                          pl = TRUE,
                          nitt=150000, burnin=50000
                          )

plot(female_LBS_model)
summary(female_LBS_model)

