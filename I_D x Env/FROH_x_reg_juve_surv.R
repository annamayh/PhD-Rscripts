library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)
library(collapse)

setwd("H:/")

## effects of spatial region on juvenile survival (0-2)
surv_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)

surv_loc_df=surv_loc_df%>%
  select(-MumFROH)%>%
  na.omit()

head(surv_loc_df)
table(surv_loc_df$MotherStatus)

surv_loc_df$Code=as.factor(surv_loc_df$Code)
surv_loc_df$MumCode=as.factor(surv_loc_df$MumCode)
surv_loc_df$BirthYear=as.factor(surv_loc_df$BirthYear)
surv_loc_df$Sex=as.factor(surv_loc_df$Sex)
surv_loc_df$MotherStatus=as.factor(surv_loc_df$MotherStatus)
surv_loc_df$Reg=as.factor(surv_loc_df$Reg)


#fitting simple model of juvenile survival 
suv_model_simple=glmmTMB(juvenile_survival~ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH+BirthWt+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial, 
                         data=surv_loc_df, 
                         na.action = na.omit,
)

summary(suv_model_simple)

surv_reg_fixed=update(suv_model_simple, ~ . + Reg) ##just region as fixed effect
summary(surv_reg_fixed)

emmeans(surv_reg_fixed, pairwise ~  Reg)



ggpredict(surv_reg_fixed, terms = c("Sex"))


reg_pred_f=ggpredict(surv_reg_fixed, terms = c("Reg", "Sex[1]"))%>%
  mutate(x = fct_relevel(x, "SI", "IM", "LA", "NG","MG",  "SG"))

reg_pred_m=ggpredict(surv_reg_fixed, terms = c("Reg", "Sex[2]"))%>%
  mutate(x = fct_relevel(x, "SI", "IM", "LA", "NG","MG",  "SG"))

reg_pred=rbind(reg_pred_f,reg_pred_m)%>%
arrange(x)
  #%>%mutate(group="Juvenile survival")
reg_pred

surv_reg=reg_pred%>%
  ggplot(aes(x=x, y=predicted, color=x, ymin=conf.low, ymax=conf.high, group=group))+
  geom_pointrange(linewidth=1, position = position_dodge(width=0.2))+
  theme_bw()+
  scale_color_manual(values = c("#f0a30a" ,"#a20025","#00aba9","chocolate1", "#60a917","#647687"))+
  labs(x="Spatial region", y="Predicted juvenile survival probability", tag="C")+
  theme(text = element_text(size = 18),legend.position = "none")
surv_reg


# INLA plot taken from INLA_juvenile_surv.R code
# inla_and_reg=surv_reg+inla_surv_plot
# inla_and_reg
# 
# ggsave(inla_and_reg,
#        file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/juve_surv_both.jpeg",
#        width = 11,
#        height = 7, 
#        dpi=300)
# 
# 





#### INTERACTION BETWEEN REGION AND FROH ####################

# 
# # 
# surv_froh=update(suv_model_simple, ~ . + FROH) #just FROH as fixed effect - we already know this tho
# summary(surv_froh)

## now fitting an interaction between region and FROH
surv_froh_inter=update(suv_model_simple, ~ . -FROH + (Reg*FROH)) #interaction between region and froh
summary(surv_froh_inter)
ggpredict(surv_froh_inter, terms = c("FROH[all]","Reg"))
#plot predictions
inter=plot(ggpredict(surv_froh_inter, terms = c("FROH[all]","Reg", "Sex")),show.title=FALSE, colors="metro", line.size=1)+
  labs(x = expression(F["ROH"]), y = "Juvenile survival probability", colour = "Spatial \nregion")+
  theme(text = element_text(size = 15),legend.position = "none") 
inter


FROH_trend=emtrends(surv_froh_inter, pairwise ~  Reg|Sex, var=c("FROH"))
test(FROH_trend$emtrends)

(FROH_trend$contrasts)




inters=inter/bw_inter & theme(legend.position = "right")
inters=inters+ plot_layout(guides = "collect")+plot_annotation(tag_levels = 'A')

ggsave(inters,
       file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/surv_FROH_interaction.jpeg",
       width = 6,
       height = 9)


