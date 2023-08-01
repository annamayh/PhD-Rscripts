
library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)

setwd("H:/")

## effects of spatial region on juvenile survival (0-2)
surv_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)
head(surv_loc_df)

surv_loc_df$BirthYear=as.numeric(surv_loc_df$BirthYear)


FROH_reg=surv_loc_df%>%
  select(Code, BirthYear, Sex, MumCode, FROH, Reg, juvenile_survival)%>%
  filter(Sex!=3)%>%
  na.omit()%>%
  mutate(year_cont=BirthYear-min(BirthYear))#



FROH_reg$Code=as.factor(FROH_reg$Code)
FROH_reg$MumCode=as.factor(FROH_reg$MumCode)
FROH_reg$BirthYear=as.factor(FROH_reg$BirthYear)
FROH_reg$Sex=as.factor(FROH_reg$Sex)
FROH_reg$Reg=as.factor(FROH_reg$Reg)



FROH_model_simple=glmmTMB(FROH~ year_cont+ Reg+
                           (1|MumCode), 
                         family=gaussian(), 
                         data=FROH_reg, 
                         na.action = na.omit,
)

summary(FROH_model_simple)
plot(ggpredict(FROH_model_simple, terms = c("Reg")))

ggpredict(FROH_model_simple, terms = c("Reg"))

FROH_reg_pred=ggpredict(FROH_model_simple, terms = c("Reg"))
reg=FROH_reg_pred%>%
  mutate(x = fct_relevel(x, "SI", "IM", "LA", "NG","MG",  "SG"))%>%
    ggplot( aes(x=x, y=predicted, color=x, ymin=conf.low, ymax=conf.high))+
  geom_pointrange(linewidth=1)+
  theme_bw()+
  scale_color_manual(values = c("#f0a30a" ,"#a20025","#00aba9","chocolate1","#60a917", "#647687"))+
  labs(x="Spatial region", y=expression(paste("F"[ROH])))+
  theme(text = element_text(size = 18),legend.position = "none")
reg


##inla_froh_gg loaded from INLA_FROH.R 
library(patchwork)
inla_and_glmm=(reg|inla_froh_gg)+plot_annotation(tag_levels = "A")
inla_and_glmm


ggsave(inla_and_glmm,
       file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/FROH.jpeg",
       width = 11,
       height = 7,
       dpi=1000)












#################################################################################
################################################################################
#### compare to only ids that survived ################
################################################################
FROH_reg_surv=FROH_reg%>%filter(juvenile_survival==1)
FROH_reg_surv$Code=as.factor(FROH_reg_surv$Code)
FROH_reg_surv$MumCode=as.factor(FROH_reg_surv$MumCode)
FROH_reg_surv$BirthYear=as.factor(FROH_reg_surv$BirthYear)
FROH_reg_surv$Sex=as.factor(FROH_reg_surv$Sex)
FROH_reg_surv$Reg=as.factor(FROH_reg_surv$Reg)


FROH_model_surv=glmmTMB(FROH~ Sex + year_cont+Reg+
                            (1|MumCode), 
                          family=gaussian(), 
                          data=FROH_reg_surv, 
                          na.action = na.omit,
)

summary(FROH_model_surv)
plot(ggpredict(FROH_model_surv, terms = c("Reg")))

FROH_surv_pred=ggpredict(FROH_model_surv, terms = c("Reg"))%>%mutate(group="2")



both=rbind(FROH_reg_pred, FROH_surv_pred)


both_reg=ggplot(both, aes(x=x, y=predicted, color=x, ymin=conf.low, ymax=conf.high))+
  geom_pointrange(position=position_dodge2(width = 0.25))+
  theme_bw()+
  scale_color_manual(values = c("#a20025","#00aba9","#60a917","chocolate1", "#647687","#f0a30a" ))+
  labs(x="Spatial region", y="FROH (all Vs surviving to age 2)", tag="A")+
  theme(text = element_text(size = 15),legend.position = "none")

both_reg



ggsave(both_reg,
       file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/region_FROH_allVsurvived.png",
       width = 8,
       height = 6)
