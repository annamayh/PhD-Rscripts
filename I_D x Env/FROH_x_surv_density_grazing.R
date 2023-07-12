library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)

setwd("H:/")

surv_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)%>%
  select(-BirthWt, -MumFROH)%>%
  na.omit()


density=read.csv("PhD_4th_yr/Spatial_var_inbreeding/glmmTMB/ValuesforAnna.csv", header=T, stringsAsFactors = F)%>%
  rename(Code=Name, BirthYear=Year)


loc_den_graze=surv_loc_df%>%
  inner_join(density)%>%
  unique()

hist(loc_den_graze$GrazeType)

head(loc_den_graze)

loc_den_graze$Code=as.factor(loc_den_graze$Code)
loc_den_graze$MumCode=as.factor(loc_den_graze$MumCode)
loc_den_graze$BirthYear=as.factor(loc_den_graze$BirthYear)
loc_den_graze$Sex=as.factor(loc_den_graze$Sex)
loc_den_graze$MotherStatus=as.factor(loc_den_graze$MotherStatus)
loc_den_graze$Reg=as.factor(loc_den_graze$Reg)



suv_model_simple=glmmTMB(juvenile_survival~ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial, 
                         data=loc_den_graze, 
                         na.action = na.omit,
)


surv_den=update(suv_model_simple, ~ . + AnnualDensity) ##just region as fixed effect
summary(surv_den)
den_pred=ggpredict(surv_den, terms = c("AnnualDensity[all]"))
den_pred
plot(den_pred)


surv_graze=update(suv_model_simple, ~ . + GrazeType) ##just region as fixed effect
summary(surv_graze)
graze_pred=ggpredict(surv_graze, terms = c("GrazeType[all]"))
graze_pred
plot(graze_pred)


surv_graze_den=update(suv_model_simple, ~ . + GrazeType+AnnualDensity) ##just region as fixed effect
summary(surv_graze_den)
graze_pred=ggpredict(surv_graze_den, terms = c("GrazeType[all]"))
graze_pred
plot(graze_pred)




surv_den_inter=update(suv_model_simple, ~ . -FROH+ FROH:AnnualDensity+ FROH:GrazeType) ##just region as fixed effect
summary(surv_den_inter)
den_pred_inter=ggpredict(surv_den_inter, 
                         terms = c("FROH[all]","AnnualDensity[0, 0.00025,0.0005, 0.00075, 0.001, 0.00125,0.0015, 0.00175, 0.002, 0.00225, 0.0025,0.00275, 0.003]"))
den_pred_inter

library(RColorBrewer)
cols <- (colorRampPalette(brewer.pal(9, "YlGnBu"))(20))


denisty_pred_inter_plot=ggplot(den_pred_inter, aes(x=x, y=predicted, color=group))+
  geom_line(linewidth=1)+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, color=NULL, fill=group),alpha=0.4, show.legend = F)+
  theme_bw()+
  scale_fill_manual(values=cols)+
  scale_color_manual(values = cols)+
  labs(x=expression(F["ROH"]), y="Juvenile survival probability",color="Population \ndensity", tag = "A")+
  theme(text = element_text(size = 18))
denisty_pred_inter_plot


denisty_pred_inter_plot=ggplot(den_pred_inter, aes(x=x, y=group))+
  geom_density_2d()+
  theme_bw()+
  #labs(x=expression(F["ROH"]), y="Juvenile survival probability",color="Population \ndensity", tag = "A")+
  theme(text = element_text(size = 18))
denisty_pred_inter_plot




cols_5 <- brewer.pal(5, "YlGnBu")
cols_5


graze_pred_inter=ggpredict(surv_den_inter, terms = c("FROH[all]","GrazeType[0,0.25,0.5,0.75,1]"))
graze_pred_inter

graze_pred_inter_plot=ggplot(graze_pred_inter, aes(x=x, y=predicted, color=group))+
  geom_line(linewidth=1.5)+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, color=NULL, fill=group),alpha=0.15, show.legend = F)+
  theme_bw()+
  scale_fill_manual(values=cols_5)+
  scale_color_manual(values = cols_5)+
  labs(x=expression(F["ROH"]), y="Juvenile survival probability",color="Graze \nType", tag="B")+
  theme(text = element_text(size = 18))
graze_pred_inter_plot


den_graze=denisty_pred_inter_plot+graze_pred_inter_plot
den_graze


ggsave(den_graze,
       file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/surv_FROH_inter_graze_population_2.png",
       width = 13,
       height = 7)



## check ass of denisty and grazing with region

deni_reg=glmmTMB(AnnualDensity~Reg+GrazeType,
                         family=gaussian, 
                         data=loc_den_graze, 
                         na.action = na.omit)
summary(deni_reg)
reg_den_pred=ggpredict(deni_reg, terms = c("Reg[all]"))
reg_den_pred
reg_den_plot=plot(reg_den_pred)

graze_reg=glmmTMB(GrazeType~Reg,
                 family=gaussian, 
                 data=loc_den_graze, 
                 na.action = na.omit)
summary(graze_reg)
reg_graze_pred=ggpredict(graze_reg, terms = c("Reg[all]"))
reg_graze_pred
reg_graz_plot=plot(reg_graze_pred)

reg_den_graze=reg_graz_plot+reg_den_plot
reg_den_graze

ggsave(reg_den_graze,
       file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/den_graze_region.png",
       width = 8,
       height = 5)

