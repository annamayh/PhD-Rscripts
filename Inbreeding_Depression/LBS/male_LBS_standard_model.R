## MALE LBS ###
library(RODBC)
library(tidyverse)
library(MCMCglmm)
library(MasterBayes)


db<-"C:\\Users\\s1881212\\Documents\\Deer_database_08_2021\\RedDeer1.99.accdb" #open connection
con<-odbcConnectAccess2007(db)

ped<-sqlFetch(con, "sys_Pedigree")%>%dplyr::select(Code, BirthYear,MumCode,Sire)
  
life<-sqlFetch(con, "tbllife") %>% 
 dplyr::select(Code,DeathYear,Sex)

odbcClose(con)


sire_count=ped%>%dplyr::select(Code,Sire)%>%
  group_by(Sire)%>%tally()%>% ##grouping by sire and counting up ids sired by id
  rename(Code=Sire, LBS=n)


dead_males=life%>%
  filter(Sex=="2")%>% #select males
  na.omit() ##removing males that are still alive (i.e. no death recorded)


all_male_LBS=left_join(dead_males,sire_count)
all_male_LBS[is.na(all_male_LBS)] <- 0 #males with no recorded sired ids given 0 LBS


setwd("H:/")

FROH_full<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom.indiv", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,KB) %>% dplyr::rename(Code=IID)%>%mutate(FROH=KB/2591865)


male_LBS_df=all_male_LBS%>%
  plyr::join(FROH_full)%>%plyr::join(ped)%>%
  dplyr::select(Code,LBS,BirthYear,MumCode,FROH)%>%
  na.omit()




##these are okay?
prior.2 = list(R = list(V = diag(2), nu = 2,fix = 2),
               G = list(G1 = list(V = diag(2), nu = 2),
                        G2 = list(V = diag(2), nu = 2)))

#again okay but not amazing for random 
prior.3 = list(R = list(V = diag(2), nu = 2,fix = 2),
               G = list(G1 = list(V = diag(2), nu = 0.002),
                        G2 = list(V = diag(2), nu = 0.002)))





male_LBS_model=MCMCglmm(LBS~ trait+trait:FROH-1,
                          #first fixed term = intercept of the Poisson process 
                          #second fixed term = intercept of the zero-inflation  
                          #remaining terms are trait (Poi & Zi) specific contrasts from the intercept.
                          random=~idh(trait):MumCode+idh(trait):BirthYear, 
                          rcov = ~ idh(trait):units,
                          data = male_LBS_df, 
                          family = "zipoisson", 
                          prior = prior.2,
                          #verbose = FALSE, 
                          pr = TRUE, 
                          pl = TRUE,
                          nitt=150000, burnin=50000, #thin = 500
                            )



plot(male_LBS_model)
summary(male_LBS_model)


# 
# Iterations = 50001:149991
# Thinning interval  = 10
# Sample size  = 10000 
# 
# DIC: 2268.91 
# 
# G-structure:  ~idh(trait):MumCode
# 
# post.mean l-95% CI u-95% CI eff.samp
# traitLBS.MumCode       0.5110   0.1909   0.8814    301.4
# traitzi_LBS.MumCode    0.8749   0.1789   1.7390    119.8
# 
# ~idh(trait):BirthYear
# 
# post.mean l-95% CI u-95% CI eff.samp
# traitLBS.BirthYear       0.3339   0.1267   0.5915    595.7
# traitzi_LBS.BirthYear    2.0786   0.7338   3.6788    293.7
# 
# R-structure:  ~idh(trait):units
# 
# post.mean l-95% CI u-95% CI eff.samp
# traitLBS.units        0.668   0.3338    1.038    298.9
# traitzi_LBS.units     1.000   1.0000    1.000      0.0
# 
# Location effects: LBS ~ trait + trait:FROH - 1 
# 
# post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# traitLBS            2.3858   1.8242   2.9105    456.0 <1e-04 ***
#   traitzi_LBS         1.5270   0.6893   2.3640    426.5 0.0004 ***
#   traitLBS:FROH     -24.7843 -34.2809 -14.8155    308.7 <1e-04 ***
#   traitzi_LBS:FROH   -1.3056 -14.7156  11.8237    259.9 0.8794   


## no effect on zi ?? = not what jisca found ... maybe try using hurdle model to see if that replicates her results?
## so with hurdle model ... is showing hurdle to be significant?
#i think it should be a zi poisson .. and not a hurdle

# traithu_LBS:FROH   10.1794   1.9981  19.1477    970.5 0.0156 *  
  

