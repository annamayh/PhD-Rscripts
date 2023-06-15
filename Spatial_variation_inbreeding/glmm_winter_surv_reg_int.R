library(tidyverse)
library(glmmTMB)
library(ggeffects)

setwd("H:/")

winter_surv_reg_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/winter_survival_loc.txt", sep = ",", header = TRUE)

Fgrm=read.table("PhD_4th_yr/2023_ROH_search/Fgrm_052023.ibc", header = TRUE)%>%
  select(IID, Fhat3)%>%rename(Code=IID, Fgrm=Fhat3)
mum_Fgrm=read.table("PhD_4th_yr/2023_ROH_search/Fgrm_052023.ibc", header = TRUE)%>%
  select(IID, Fhat3)%>%rename(MumCode=IID, MumFgrm=Fhat3)
#hist(Fgrm$Fgrm, breaks = 30)

winter_surv_loc_df=winter_surv_reg_df%>%left_join(Fgrm)%>%left_join(mum_Fgrm)%>%
  filter(winter_survival!="NA")%>%
  na.omit()

head(winter_surv_loc_df)
table(winter_surv_loc_df$winter_survival)

winter_surv_loc_df$winter_survival=as.integer(winter_surv_loc_df$winter_survival)
winter_surv_loc_df$Code=as.factor(winter_surv_loc_df$Code)
winter_surv_loc_df$MumCode=as.factor(winter_surv_loc_df$MumCode)
winter_surv_loc_df$BirthYear=as.factor(winter_surv_loc_df$BirthYear)
winter_surv_loc_df$Sex=as.factor(winter_surv_loc_df$Sex)
winter_surv_loc_df$MotherStatus=as.factor(winter_surv_loc_df$MotherStatus)
winter_surv_loc_df$Reg=as.factor(winter_surv_loc_df$Reg)


#fitting simple model of juvenile survival 
wint_suv_model_simple=glmmTMB(winter_survival~ Sex + MotherStatus + mum_age+mum_age_sq+BirthWt+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial, 
                         data=winter_surv_loc_df, 
                         na.action = na.omit,
)

summary(wint_suv_model_simple)


surv_reg_fixed=update(wint_suv_model_simple, ~ . + Reg)
summary(surv_reg_fixed)
## plot showing difference in survival across regions 
plot(ggpredict(surv_reg_fixed, terms = c("Reg")))
ggpredict(surv_reg_fixed, terms = c("Reg"))

surv_reg_fixed_froh=update(surv_reg_fixed, ~ . + FROH)
summary(surv_reg_fixed_froh)
## plot showing difference in survival across regions 
plot(ggpredict(surv_reg_fixed_froh, terms = c("FROH[all]","Reg")))



### survival differs across regions 
surv_FROH4=update(wint_suv_model_simple, ~ . + (Reg*FROH))
R4=plot(ggpredict(surv_FROH4, terms = c("FROH[all]","Reg")),show.title=FALSE, colors="metro", line.size=1)
R4
summary(surv_FROH4)
ggpredict(surv_FROH4, terms = c("FROH[all]","Reg"))
FROH_trend=emtrends(surv_FROH4, pairwise ~ Reg, var="FROH")
test(FROH_trend$emtrends)



surv_Fgrm4=update(wint_suv_model_simple, ~ . + (Reg*Fgrm))
G4=plot(ggpredict(surv_Fgrm4, terms = c("Fgrm[all]","Reg")), show.title=FALSE,colors="metro", line.size=1)
G4
summary(surv_Fgrm4)

winter_survival=R4+G4
winter_survival


ggsave(winter_survival,
       file = "PhD_4th_yr/Spatial_var_inbreeding/glmmTMB/wint_surv_region_interaction.png",
       width = 10,
       height = 6)
