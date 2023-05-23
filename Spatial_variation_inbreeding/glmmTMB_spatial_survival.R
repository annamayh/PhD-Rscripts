library(tidyverse)
library(glmmTMB)
library(ggeffects)

setwd("H:/")

surv_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)

Fgrm=read.table("PhD_4th_yr/2023_ROH_search/Fgrm_052023.ibc", header = TRUE)%>%
  select(IID, Fhat3)%>%rename(Code=IID, Fgrm=Fhat3)
mum_Fgrm=read.table("PhD_4th_yr/2023_ROH_search/Fgrm_052023.ibc", header = TRUE)%>%
  select(IID, Fhat3)%>%rename(MumCode=IID, MumFgrm=Fhat3)
#hist(Fgrm$Fgrm, breaks = 30)

surv_loc_df=surv_loc_df%>%left_join(Fgrm)%>%left_join(mum_Fgrm)%>%
  select(-BirthWt)%>%na.omit()

head(surv_loc_df)
table(surv_loc_df$Reg)


surv_loc_df$Code=as.factor(surv_loc_df$Code)
surv_loc_df$MumCode=as.factor(surv_loc_df$MumCode)
surv_loc_df$BirthYear=as.factor(surv_loc_df$BirthYear)
surv_loc_df$Sex=as.factor(surv_loc_df$Sex)
surv_loc_df$MotherStatus=as.factor(surv_loc_df$MotherStatus)
surv_loc_df$Reg=as.factor(surv_loc_df$Reg)


#fitting simple model of juvenile survival 
suv_model_simple=glmmTMB(juvenile_survival~ Sex + MotherStatus + mum_age+mum_age_sq+
                         (1|BirthYear)+(1|MumCode), 
                       family=binomial, 
                       data=surv_loc_df, 
                       na.action = na.omit,
                       )

summary(suv_model_simple)


surv_reg_fixed=update(suv_model_simple, ~ . + Reg)
summary(surv_reg_fixed)
## plot showing difference in survival across regions 
plot(ggpredict(surv_reg_fixed, terms = c("Reg")))

### survival differs across regions 

surv_FROH=update(suv_model_simple, ~ . + FROH )
summary(surv_FROH)

surv_model_mumFROH=update(surv_FROH, ~ . + MumFROH )
summary(surv_model_mumFROH)
ggpredict(surv_model_mumFROH, terms = c("MumFROH"))


surv_mumFROH_reg=update(surv_model_mumFROH, ~ . + (1|Reg))
summary(surv_mumFROH_reg)
ggpredict(surv_mumFROH_reg, terms = c("MumFROH"))


surv_inter=update(suv_model_simple, ~ . + FROH+(MumFROH*Reg) )
summary(surv_inter)
FROH_pred=plot(ggpredict(surv_inter, terms = c("MumFROH[all]","Reg")))
pred_mumFROH=(ggpredict(surv_inter, terms = c("MumFROH[all]","Reg")))

ggplot(pred_mumFROH, aes(x=x,y=predicted, colour=group))+
  geom_line()

plot(ggpredict(surv_inter, terms = c("FROH","Reg")))

## same with Fgrm

surv_Fgrm=update(suv_model_simple, ~ . + Fgrm)
summary(surv_Fgrm)

surv_model_mumFgrm=update(surv_Fgrm, ~ . + MumFgrm )
summary(surv_model_mumFgrm)

surv_mumFgrm_reg=update(surv_model_mumFgrm, ~ . + (1|Reg))
summary(surv_mumFgrm_reg)

surv_inter_fgrm=update(suv_model_simple, ~ . + Fgrm + (MumFgrm*Reg) )
summary(surv_inter_fgrm)

Fgrm_preed=plot(ggpredict(surv_inter_fgrm, terms = c("MumFgrm[all]","Reg")))

library(patchwork)
FROH_pred+Fgrm_preed





#### compare region as fixed effect with 2 measures of inbreeding 

surv_FROH1=update(suv_model_simple, ~ . + FROH +Reg)
R1=plot(ggpredict(surv_FROH1, terms = c("FROH[all]","Reg")),show.title=FALSE, colors="metro", line.size=1)

surv_Fgrm1=update(suv_model_simple, ~ . + Fgrm +Reg)
G1=plot(ggpredict(surv_Fgrm1, terms = c("Fgrm[all]","Reg")),show.title=FALSE, colors="metro", line.size=1)

A=R1+G1+plot_annotation(title = "Predicted survival with increasing id inbreeding coefficient") 
## shows that for both measures of inbreeding juvenile survival decreases at the same(ish) rate


surv_FROH2=update(suv_model_simple, ~ . + FROH +MumFROH)
R2=plot(ggpredict(surv_FROH2, terms = c("MumFROH[all]")),show.title=FALSE,  line.size=1)
R2

surv_Fgrm2=update(suv_model_simple, ~ . + Fgrm +MumFgrm)
G2=plot(ggpredict(surv_Fgrm2, terms = c("MumFgrm[all]")),show.title=FALSE,  line.size=1)
G2

B=R2+G2+plot_annotation(title = "Predicted survival with Mother inbreeding coefficient (no region fitted)")

surv_FROH3=update(suv_model_simple, ~ . + FROH +(1|Reg)+MumFROH)
R3=plot(ggpredict(surv_FROH3, terms = c("MumFROH[all]")),show.title=FALSE,  line.size=1)
R3
summary(surv_FROH3)


surv_Fgrm3=update(suv_model_simple, ~ . + Fgrm +(1|Reg)+MumFgrm)
G3=plot(ggpredict(surv_Fgrm3, terms = c("MumFgrm[all]")),show.title=FALSE, line.size=1)
G3
summary(surv_Fgrm3)


C=R3+G3+plot_annotation(title = "Predicted survival with region as random effect and Mother inbreeding coefficient")


surv_FROH4=update(suv_model_simple, ~ . + FROH +(Reg*MumFROH))
R4=plot(ggpredict(surv_FROH4, terms = c("MumFROH[all]","Reg")),show.title=FALSE, colors="metro", line.size=1)
R4
summary(surv_FROH4)

surv_Fgrm4=update(suv_model_simple, ~ . + Fgrm +(Reg*MumFgrm))
G4=plot(ggpredict(surv_Fgrm4, terms = c("MumFgrm[all]","Reg")), show.title=FALSE,colors="metro", line.size=1)
G4
summary(surv_Fgrm4)


D=R4+G4+plot_annotation(title = "Predicted survival with interaction between region (fixed) and Mother inbreeding coefficient")

###
# surv_FROH5=update(suv_model_simple, ~ . + FROH +Reg+MumFROH)
# R5=plot(ggpredict(surv_FROH5, terms = c("MumFROH[all]","Reg")),show.title=FALSE, colors="metro", line.size=1)
# R5
# 
# surv_Fgrm5=update(suv_model_simple, ~ . + Fgrm +Reg+MumFgrm)
# G5=plot(ggpredict(surv_Fgrm5, terms = c("MumFgrm[all]","Reg")), show.title=FALSE,colors="metro", line.size=1)
# G5
## when treat regions as a fixed effect only it just eastimates the same gradient of mumFroh for all regions (obviously)


D=R4+G4+plot_annotation(title = "Predicted survival with interaction between region (fixed) and Mother inbreeding coefficient")



A/B/C/D

ggsave(A,
       file = "PhD_4th_yr/Spatial_var_inbreeding/glmmTMB/survival_id_incIBC.png",
       width = 10,
       height = 6)

ggsave(B,
       file = "PhD_4th_yr/Spatial_var_inbreeding/glmmTMB/survival_motherIBC.png",
       width = 10,
       height = 6)

ggsave(C,
       file = "PhD_4th_yr/Spatial_var_inbreeding/glmmTMB/survival_motherIBC_plus_regionRand.png",
       width = 10,
       height = 6)

ggsave(D,
       file = "PhD_4th_yr/Spatial_var_inbreeding/glmmTMB/survival_motherIBC_inter_regionFix.png",
       width = 10,
       height = 6)



#check correlations of mum ibc in IM
surv_loc_df%>%filter(Reg=="IM")%>%
  ggplot(aes(MumFROH, MumFgrm))+
  geom_point()

cor(surv_loc_df$MumFROH, surv_loc_df$MumFgrm)
