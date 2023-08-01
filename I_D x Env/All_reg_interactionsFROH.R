## All interactions with region and FROH for survival measures and birth wt ###

library(RODBC)
library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)

setwd("H:/")
## birth weight ####
surv_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)%>%
  select(-BirthWt,-MumFROH)

db<-"C:\\Users\\s1881212\\Documents\\Deer_database_2022/RedDeer2.05.accdb" #open connection
con<-odbcConnectAccess2007(db)
birth_wt<-sqlFetch(con, "sys_BirthWt") 
odbcClose(con) #close connection

bw_df=surv_loc_df%>%full_join(birth_wt)%>%na.omit
head(bw_df) ##2378 ids in df

bw_df$Code=as.factor(bw_df$Code)
bw_df$MumCode=as.factor(bw_df$MumCode)
bw_df$BirthYear=as.factor(bw_df$BirthYear)
bw_df$Sex=as.factor(bw_df$Sex)
bw_df$Reg=as.factor(bw_df$Reg)

bw_model_inter=glmmTMB(CaptureWt~ Sex + AgeHrs+ MotherStatus+mum_age+mum_age_sq+Day_seq+
                         (FROH*Reg)+
                          (1|BirthYear)+ (1|MumCode), 
                        family=gaussian(), 
                        data=bw_df, 
                        na.action = na.omit,
)

summary(bw_model_inter)

#AIC w/o interaction
# AIC      BIC   logLik deviance df.resid 
# 6702.8   6812.5  -3332.4   6664.8     2359 
# #with inter
# AIC      BIC   logLik deviance df.resid 
# 6706.8   6845.4  -3329.4   6658.8     2354 


bw_inter=plot(ggpredict(bw_model_inter, terms = c("FROH[all]","Reg","AgeHrs[0]")),show.title=FALSE, colors="metro", line.size=1)+
  labs(x = expression(F["ROH"]), y = "Birth weight (kg)", colour = "Spatial \nregion")+
  theme(text = element_text(size = 15),legend.position = "none") 
bw_inter

ggpredict(bw_model_inter, terms = c("FROH[all]","Reg","AgeHrs[0]"))

bw_trends=emtrends(bw_model_inter, pairwise ~ Reg, var="FROH")
bw_trends
test(bw_trends$emtrends)
7.09 -6.54 
#############################
##### winter survival ######
############################

winter_surv_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/winter_survival_loc.txt", sep = ",", header = TRUE)

winter_surv_df=winter_surv_df%>%
  select(-BirthWt, -MumFROH)%>%
  na.omit()#2403 ids

head(winter_surv_df)
table(winter_surv_df$MotherStatus)

winter_surv_df$Code=as.factor(winter_surv_df$Code)
winter_surv_df$MumCode=as.factor(winter_surv_df$MumCode)
winter_surv_df$BirthYear=as.factor(winter_surv_df$BirthYear)
winter_surv_df$Sex=as.factor(winter_surv_df$Sex)
winter_surv_df$MotherStatus=as.factor(winter_surv_df$MotherStatus)
winter_surv_df$Reg=as.factor(winter_surv_df$Reg)


#fitting simple model of juvenile survival 
winter_surv_inter=glmmTMB(winter_survival~ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+
                            (FROH*Reg)+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial, 
                         data=winter_surv_df, 
                         na.action = na.omit,
)

summary(winter_surv_inter)

#AIC
#w/o
# AIC      BIC   logLik deviance df.resid 
# 2074.0   2172.4  -1020.0   2040.0     2386 
# #with 
# AIC      BIC   logLik deviance df.resid 
# 2079.1   2206.4  -1017.5   2035.1     2381 

ggpredict(winter_surv_inter, terms = c("FROH[all]","Reg","MotherStatus[True yeld]"))
#plot predictions
winter_inter=plot(ggpredict(winter_surv_inter, terms = c("FROH[all]","Reg","MotherStatus[True yeld]")),show.title=FALSE, colors="metro", line.size=1)+
  labs(x = expression(F["ROH"]), y = "Winter survival probability", colour = "Spatial \nregion")+
  theme(text = element_text(size = 15)) 
winter_inter

winter_trends=emtrends(winter_surv_inter, pairwise ~  Reg, var=c("FROH"))
winter_trends
test(winter_trends$emtrends)
##########################
### juvenile survival ####
##########################
## effects of spatial region on juvenile survival (0-2)
surv_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)

surv_loc_df=surv_loc_df%>%
  select(-MumFROH, -BirthWt)%>%
  na.omit()#2621

head(surv_loc_df)
table(surv_loc_df$MotherStatus)

surv_loc_df$Code=as.factor(surv_loc_df$Code)
surv_loc_df$MumCode=as.factor(surv_loc_df$MumCode)
surv_loc_df$BirthYear=as.factor(surv_loc_df$BirthYear)
surv_loc_df$Sex=as.factor(surv_loc_df$Sex)
surv_loc_df$MotherStatus=as.factor(surv_loc_df$MotherStatus)
surv_loc_df$Reg=as.factor(surv_loc_df$Reg)


#fitting simple model of juvenile survival 
suv_model_inter=glmmTMB(juvenile_survival~ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+
                           (FROH*Reg)+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial, 
                         data=surv_loc_df, 
                         na.action = na.omit,
)

summary(suv_model_inter)

#AIC 
#w/o
# AIC      BIC   logLik deviance df.resid 
# 2969.5   3069.3  -1467.8   2935.5     2604
# #with inter
# AIC      BIC   logLik deviance df.resid 
# 2975.0   3104.2  -1465.5   2931.0     2599 

ggpredict(suv_model_inter, terms = c("FROH[all]","Reg","MotherStatus[True yeld]"))
#plot predictions
inter=plot(ggpredict(suv_model_inter, terms = c("FROH[all]","Reg","MotherStatus[True yeld]")),show.title=FALSE, colors="metro", line.size=1)+
  labs(x = expression(F["ROH"]), y = "Juvenile survival probability", colour = "Spatial \nregion")+
  theme(text = element_text(size = 15),legend.position = "none") 
inter

surv_trends=emtrends(suv_model_inter, pairwise ~  Reg, var=c("FROH"))
surv_trends
test(surv_trends$emtrends)


## all plots ###

inters=(bw_inter/winter_inter/inter)+plot_annotation(tag_levels = "A")
inters


ggsave(inters,
       file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/All_FROH_interaction.jpeg",
       width = 6,
       height = 10)

