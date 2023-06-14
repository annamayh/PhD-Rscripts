library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)

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


#### compare region as fixed effect with 2 measures of inbreeding 

surv_FROH1=update(suv_model_simple, ~ . + FROH +Reg)
R1=plot(ggpredict(surv_FROH1, terms = c("FROH[all]","Reg")),show.title=FALSE, colors="metro", line.size=1)
R1
summary(surv_FROH1)

surv_Fgrm1=update(suv_model_simple, ~ . + Fgrm +Reg)
G1=plot(ggpredict(surv_Fgrm1, terms = c("Fgrm[all]","Reg")),show.title=FALSE, colors="metro", line.size=1)
G1
summary(surv_Fgrm1)

A=R1+G1+plot_annotation(title = "Predicted survival with increasing id inbreeding coefficient") 
## shows that for both measures of inbreeding juvenile survival decreases at the same(ish) rate
 
## could also assume an interaction between region and focal inbreeding 
surv_FROH1_inter=update(suv_model_simple, ~ . + FROH *Reg)
R1_I=plot(ggpredict(surv_FROH1_inter, terms = c("FROH[all]","Reg")),show.title=FALSE, colors="metro", line.size=1)
R1_I
summary(surv_FROH1_inter)
FROH_trend=emtrends(surv_FROH1_inter, pairwise ~ Reg, var="FROH")
test(FROH_trend$emtrends)

surv_Fgrm1_inter=update(suv_model_simple, ~ . + Fgrm *Reg)
G1_I=plot(ggpredict(surv_Fgrm1_inter, terms = c("Fgrm[all]","Reg")),show.title=FALSE, colors="metro", line.size=1)
G1_I
summary(surv_Fgrm1_inter)
emm1.2 = emmeans(surv_Fgrm1_inter, specs = pairwise ~ Reg)
emm1.2$contrasts

A_I=R1_I+G1_I+plot_annotation(title = "Predicted survival with increasing id inbreeding coefficient (interaction of inbreeding and region)") 
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
R4=plot(ggemmeans(surv_FROH4, terms = c("MumFROH[all]","Reg")),show.title=FALSE, colors="metro", line.size=1)
R4
summary(surv_FROH4)

surv_Fgrm4=update(suv_model_simple, ~ . + Fgrm +(Reg*MumFgrm))
G4=plot(ggemmeans(surv_Fgrm4, terms = c("MumFgrm[all]","Reg")), show.title=FALSE,colors="metro", line.size=1)
G4
summary(surv_Fgrm4)

#pairwise comparisons
emm1.1 = emmeans(surv_FROH4, specs = pairwise ~ MumFROH:Reg)
emm1.1$contrasts


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

ggsave(A_I,
       file = "PhD_4th_yr/Spatial_var_inbreeding/glmmTMB/survival_id_incIBC_inter_region.png",
       width = 10,
       height = 6)




#check correlations of mum ibc in IM
intermeiate=surv_loc_df%>%filter(Reg=="SG")#%>%
  ggplot(aes(MumFROH, MumFgrm))+
  geom_point()
length(unique(check$MumCode))
  
cor(surv_loc_df$MumFROH, surv_loc_df$MumFgrm)

##checking models fitted for each reg seperatly
intermeiate=surv_loc_df%>%filter(Reg=="SG")#%>%

suv_model_simple_intermed=glmmTMB(juvenile_survival~ Sex + MotherStatus + mum_age+mum_age_sq+FROH+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial, 
                         data=intermeiate, 
                         na.action = na.omit,
)

summary(suv_model_simple_intermed)
##FROH not sig for SG but highly sig for IM

plot(ggpredict(suv_model_simple_intermed, terms = c("FROH[all]")),show.title=FALSE,  line.size=1)

### estimtes by hand to double check ####

new_FROH <- rnorm(nrow(surv_loc_df), mean(surv_loc_df$FROH), sd(surv_loc_df$FROH))
new_mumFROH <- rnorm(nrow(surv_loc_df), mean(surv_loc_df$MumFROH), sd(surv_loc_df$MumFROH))
new_

head(surv_loc_df)

?rnorm
