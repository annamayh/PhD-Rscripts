library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)

setwd("H:/")

## effects of spatial region on juvenile survival (0-2)
surv_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)

surv_loc_df=surv_loc_df%>%
  select(-MumFROH)%>%
  na.omit()

head(surv_loc_df)
table(surv_loc_df$Reg)

surv_loc_df$Code=as.factor(surv_loc_df$Code)
surv_loc_df$MumCode=as.factor(surv_loc_df$MumCode)
surv_loc_df$BirthYear=as.factor(surv_loc_df$BirthYear)
surv_loc_df$Sex=as.factor(surv_loc_df$Sex)
surv_loc_df$MotherStatus=as.factor(surv_loc_df$MotherStatus)
surv_loc_df$Reg=as.factor(surv_loc_df$Reg)


#fitting simple model of juvenile survival 
suv_model_simple=glmmTMB(juvenile_survival~ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial, 
                         data=surv_loc_df, 
                         na.action = na.omit,
)


surv_reg_fixed=update(suv_model_simple, ~ . + Reg) ##just region as fixed effect
summary(surv_reg_fixed)
reg_surv=plot(ggpredict(surv_reg_fixed, terms = c("Reg")),colors="Red",  dot.size = 3, show.title = FALSE)+
  labs(x = "Spatial region", y = "Juvenile survival probability")+
  theme(text = element_text(size = 18), legend.position = "none") 
reg_surv
reg_pred=ggpredict(surv_reg_fixed, terms = c("Reg"))%>%mutate(group="Juvenile survival")
reg_pred

surv_froh=update(suv_model_simple, ~ . + FROH + Reg) #just FROH as fixed effect - we already know this tho
summary(surv_froh)

## now fitting an interaction between region and FROH
surv_froh_inter=update(suv_model_simple, ~ . + (Reg*FROH)) #interaction between region and froh
summary(surv_froh_inter)
#plot predictions
inter=plot(ggpredict(surv_froh_inter, terms = c("FROH[all]","Reg")),show.title=FALSE, colors="metro", line.size=1)+
  labs(x = expression(F["ROH"]), y = "Juvenile survival probability", colour = "Spatial region")+
  theme(text = element_text(size = 15)) 
inter

ggpredict(surv_froh_inter, terms = c("FROH[all]","Reg"))

### same process to look at winter survival 

winter_surv_reg_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/first_year_survival_loc.txt", sep = ",", header = TRUE)

winter_surv_loc_df=winter_surv_reg_df%>%
  filter(first_year_survival!="NA")%>%
  select(-MumFROH)%>%
  na.omit()

head(winter_surv_loc_df)
table(winter_surv_loc_df$winter_survival)

winter_surv_loc_df$first_year_survival=as.integer(winter_surv_loc_df$first_year_survival)
winter_surv_loc_df$Code=as.factor(winter_surv_loc_df$Code)
winter_surv_loc_df$MumCode=as.factor(winter_surv_loc_df$MumCode)
winter_surv_loc_df$BirthYear=as.factor(winter_surv_loc_df$BirthYear)
winter_surv_loc_df$Sex=as.factor(winter_surv_loc_df$Sex)
winter_surv_loc_df$MotherStatus=as.factor(winter_surv_loc_df$MotherStatus)
winter_surv_loc_df$Reg=as.factor(winter_surv_loc_df$Reg)


#fitting simple model of juvenile survival 
wint_suv_model_simple=glmmTMB(first_year_survival~ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+ ##fit birth weight as cov as 
                                (1|BirthYear)+(1|MumCode), 
                              family=binomial, 
                              data=winter_surv_loc_df, 
                              na.action = na.omit,
)

summary(wint_suv_model_simple)


plot(ggpredict(wint_suv_model_simple, terms = c("Day_seq[all]")))
ggpredict(wint_suv_model_simple, terms = c("Day_seq[17, 47, 77, 108, 139]"))
plot(ggpredict(wint_suv_model_simple, terms = c("BirthWt[all]")))
ggpredict(wint_suv_model_simple, terms = c("BirthWt[all]"))


wint_surv_reg_fixed=update(wint_suv_model_simple, ~ . + Reg) ##just region as fixed effect
summary(wint_surv_reg_fixed)
wint_reg=plot(ggpredict(wint_surv_reg_fixed, terms = c("Reg")),colors="Red",  dot.size = 3, show.title = FALSE)+
  labs(x = "Spatial region", y = "First year survival probability")+
  theme(text = element_text(size = 18), legend.position = "none") 
wint_reg
wint_reg_pred=ggpredict(wint_surv_reg_fixed, terms = c("Reg"))%>%mutate(group="First year survival")
wint_reg_pred

wint_surv_froh=update(wint_suv_model_simple, ~ . + FROH) #FROH as fixed effect - we already know this tho
summary(wint_surv_froh)
wint_surv_froh=update(wint_suv_model_simple, ~ . + FROH+Reg) #FROH as fixed effect - we already know this tho
summary(wint_surv_froh)

## now fitting an interaction between region and FROH
wint_surv_froh_inter=update(wint_suv_model_simple, ~ . + (Reg*FROH)) #interaction between region and froh
summary(wint_surv_froh_inter)
plot(ggpredict(wint_surv_froh_inter, terms = c("FROH[all]","Reg","Sex")))
##plot predictions
wint_inter=plot(ggpredict(wint_surv_froh_inter, terms = c("FROH[all]","Reg")),show.title=FALSE, colors="metro", line.size=1)+
  labs(x = expression(F["ROH"]), y = "First year survival probability", colour = "Spatial region")+
  theme(text = element_text(size = 15),legend.position = "none")
wint_inter
FROH_trend=emtrends(wint_surv_froh_inter, pairwise ~ Reg, var="FROH") ##all pairwise comparisons 
FROH_trend
test(FROH_trend$emtrends)



#plot wth winter and juve survival 
both_inter=wint_inter+inter+plot_annotation(tag_levels = 'A')
both_inter
  
ggsave(both_inter,
       file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/froh_inter_wint_juve_surv.png",
       width = 10,
       height = 5)


## plot w/ differences bteween regions

both_df_reg=rbind(reg_pred,wint_reg_pred)

both_reg=ggplot(both_df_reg, aes(x=x, y=predicted, color=x))+
  geom_point()+
  geom_pointrange(aes(ymin=conf.low, ymax=conf.high))+
  theme_bw()+
  facet_wrap(~group)+
  scale_color_manual(values = c("#a20025","#00aba9","#60a917","chocolate1", "#647687","#f0a30a" ))+
  labs(x="Spatial region", y="Predicted survival probability")+
  theme(text = element_text(size = 18),legend.position = "none")


both_reg

ggsave(both_reg,
       file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/region_surv.png",
       width = 8,
       height = 5)




SG_look=winter_surv_loc_df%>%filter(Reg=="SG")
table(SG_look$first_year_survival)



