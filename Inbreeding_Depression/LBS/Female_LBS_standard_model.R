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


# 
# #trace plots are a bit shit
# prior.1 = list(R = list(V = diag(2), nu = 0.002,fix = 2),
#                    G = list(G1 = list(V = diag(2), nu = 0.002),
#                             G2 = list(V = diag(2), nu = 0.002)))

##these are okay?
prior.2 = list(R = list(V = diag(2), nu = 2,fix = 2),
               G = list(G1 = list(V = diag(2), nu = 2),
                        G2 = list(V = diag(2), nu = 2)))

#again okay but not amazing for random 
prior.3 = list(R = list(V = diag(2), nu = 2,fix = 2),
               G = list(G1 = list(V = diag(2), nu = 0.002),
                        G2 = list(V = diag(2), nu = 0.002)))

## like really really shit 
# prior.4 = list(R = list(V = diag(2), nu = 0.002),
#                G = list(G1 = list(V = diag(2), nu = 0.002),
#                         G2 = list(V = diag(2), nu = 0.002)))
# 

#what does diag mean?
#what does at.level() do? when i add it in, priors dont work - think it if you have multiple 
#do i need G1, G2 etc cos on others ive just done G1, G1 

female_LBS_model=MCMCglmm(LBS~ trait+trait:FROH-1,
                          #first fixed term = intercept of the Poisson process 
                          #second fixed term = intercept of the zero-inflation  
                          #remaining terms are trait (Poi & Zi) specific contrasts from the intercept.
                          random=~idh(trait):MumCode+idh(trait):BirthYear, 
                          rcov = ~ idh(trait):units,
                          data = female_LBS_data, 
                          family = "zipoisson", 
                          prior = prior.2,
                          #verbose = FALSE, 
                          pr = TRUE, 
                          pl = TRUE,
                          nitt=150000, burnin=50000, #thin = 500
                          )

plot(female_LBS_model)
summary(female_LBS_model)

#thinning interval=500 for below 

# Iterations = 50001:249501
# Thinning interval  = 500
# Sample size  = 400 
# 
# DIC: 4740.741 
# 
# G-structure:  ~idh(trait):MumCode
# 
# post.mean l-95% CI u-95% CI eff.samp
# traitLBS.MumCode       0.1313  0.08318   0.1891      400
# traitzi_LBS.MumCode    0.6120  0.18874   1.0350      400
# 
# ~idh(trait):BirthYear
# 
# post.mean l-95% CI u-95% CI eff.samp
# traitLBS.BirthYear       0.4009   0.2252   0.6082    400.0
# traitzi_LBS.BirthYear    2.5040   0.9315   4.2734    328.2
# 
# R-structure:  ~idh(trait):units
# 
# post.mean l-95% CI u-95% CI eff.samp
# traitLBS.units       0.1285  0.08134    0.176      400
# traitzi_LBS.units    1.0000  1.00000    1.000        0
# 
# Location effects: LBS ~ trait + trait:FROH - 1 
# 
# post.mean l-95% CI u-95% CI eff.samp  pMCMC   
# traitLBS            1.3630   1.1390   1.6268    400.0 <0.003 **
#   traitzi_LBS        -1.5460  -2.3566  -0.9565    403.7 <0.003 **
#   traitLBS:FROH      -0.7142  -3.5967   2.0237    400.0  0.645   
# traitzi_LBS:FROH   14.1409   7.0072  21.2084    400.0 <0.003 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1