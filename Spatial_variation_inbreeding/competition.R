library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)
library(colorspace)


setwd("H:/")

surv_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/winter_survival_loc.txt", sep = ",", header = TRUE)%>%
  filter(Sex!=3)%>%
  na.omit()

density=read.csv("PhD_4th_yr/Spatial_var_inbreeding/glmmTMB/ValuesforAnna.csv", header=T, stringsAsFactors = F)%>%
  rename(Code=Name, BirthYear=Year)

loc_den_graze=surv_loc_df%>%
  inner_join(density)%>%
  unique()%>%
  na.omit()%>%
 # mutate(AnnualDensity=AnnualDensity*100)%>%
  mutate(Density_per_graze=GrazeType/AnnualDensity)


loc_den_graze$Code=as.factor(loc_den_graze$Code)
loc_den_graze$Sex=as.factor(loc_den_graze$Sex)
loc_den_graze$Reg=as.factor(loc_den_graze$Reg)
loc_den_graze$MumCode=as.factor(loc_den_graze$MumCode)
loc_den_graze$BirthYear=as.factor(loc_den_graze$BirthYear)
loc_den_graze$MotherStatus=as.factor(loc_den_graze$MotherStatus)

graze_den_simple=glmmTMB(Density_per_graze~ Sex + Reg, 
                     family="gaussian", 
                     data=loc_den_graze, 
                     na.action = na.omit,
)



summary(graze_den_simple)

graze_den_pred=ggpredict(graze_den_simple, terms = c("Reg"))%>%
  mutate(x = fct_relevel(x, "SI", "IM", "LA", "NG","MG",  "SG"))


graz_pred_plot=ggplot(graze_den_pred,aes(x=x, y=predicted, colour=x, ymin=conf.low, ymax=conf.high, group=group))+
  geom_pointrange(linewidth=1)+
  theme_bw()+
  scale_color_manual(values = c("#f0a30a" ,"#a20025","#00aba9","chocolate1", "#60a917","#647687"))+
  labs(x="Spatial region", y="Competition")+
  theme(text = element_text(size = 18),legend.position = "none")#+
# geom_boxplot(data=loc_den_graze, aes(Reg, GrazeType, colour=Reg, alpha=0.1),
#              inherit.aes = F, position=position_nudge(x=0.1), width=0.1)

graz_pred_plot







suv_model_simple=glmmTMB(winter_survival~ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH+Day_seq+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial, 
                         data=loc_den_graze, 
                         na.action = na.omit,
)

summary(suv_model_simple)



surv_graze_den_inter=update(suv_model_simple, ~ . + -FROH +FROH*Density_per_graze) ##just region as fixed effect
summary(surv_graze_den_inter)



library(viridis)
cols_5 <- viridis(8)
#cols_5


graze_pred_inter=ggpredict(surv_graze_den_inter, terms = c("FROH[all]","Density_per_graze[0,1,2,3,4,5]"))
graze_pred_inter

graze_pred_inter_plot=ggplot(graze_pred_inter, aes(x=x, y=predicted, color=group))+
  geom_line(linewidth=1)+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, color=NULL, fill=group),alpha=0.15, show.legend = F)+
  theme_bw()+
  scale_fill_manual(values=cols_5)+
  scale_color_manual(values = cols_5)+
  labs(x=expression(F["ROH"]), y="Juvenile survival probability",color="Competition")+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=rel(0.8)),
        legend.text=element_text(size=rel(0.5)))
graze_pred_inter_plot



