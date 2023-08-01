## juvenile survival vs vegetation stuff #####
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)

setwd("H:/")
surv_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)%>%
  select(-MumFROH, -BirthWt)%>%
  na.omit()

density=read.csv("PhD_4th_yr/Spatial_var_inbreeding/glmmTMB/ValuesforAnna.csv", header=T, stringsAsFactors = F)%>%
  rename(Code=Name, BirthYear=Year)

loc_den_graze=surv_loc_df%>%
  inner_join(density)%>%
  unique()%>%
  na.omit()
  #mutate(AnnualDensity=AnnualDensity*100)%>%
  #mutate(Density_per_graze=GrazeType/AnnualDensity)


hist(loc_den_graze$GrazeType)
head(loc_den_graze)

loc_den_graze$Code=as.factor(loc_den_graze$Code)
loc_den_graze$MumCode=as.factor(loc_den_graze$MumCode)
loc_den_graze$BirthYear=as.factor(loc_den_graze$BirthYear)
loc_den_graze$Sex=as.factor(loc_den_graze$Sex)
loc_den_graze$MotherStatus=as.factor(loc_den_graze$MotherStatus)
loc_den_graze$Reg=as.factor(loc_den_graze$Reg)

cor(loc_den_graze$GrazeType, loc_den_graze$AnnualDensity)

surv_graze_inter=glmmTMB(juvenile_survival~ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+
                           (FROH*GrazeType)+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial, 
                         data=loc_den_graze, 
                         na.action = na.omit,
)

summary(surv_graze_inter)

##AICs
#no inter 
# AIC      BIC   logLik deviance df.resid 
# 2484.4   2558.9  -1229.2   2458.4     2269 
#w /inter 
# AIC      BIC   logLik deviance df.resid 
# 2486.3   2566.6  -1229.2   2458.3     226

graze_pred_inter=ggpredict(surv_graze_inter, terms = c("FROH[all]","GrazeType[0,0.2,0.4,0.6,0.8,1]"))
graze_pred_inter

library(viridis)
cols_5 <- viridis(7)
graze_pred_inter_plot=ggplot(graze_pred_inter, aes(x=x, y=predicted, color=group))+
  geom_line(linewidth=1)+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, color=NULL, fill=group),alpha=0.15, show.legend = F)+
  theme_bw()+
  scale_fill_manual(values=cols_5)+
  scale_color_manual(values = cols_5)+
  labs(x=expression(F["ROH"]), y="Juvenile survival probability",color="High quality \ngrazing")+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=12),
        legend.text=element_text(size=rel(0.5)),
        axis.title.y = element_text(size=13))
graze_pred_inter_plot

### pop density ###
surv_pop_inter=glmmTMB(juvenile_survival~ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+
                           (FROH*AnnualDensity)+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial, 
                         data=loc_den_graze, 
                         na.action = na.omit,
)

summary(surv_pop_inter)

## no inter
# AIC      BIC   logLik deviance df.resid 
# 2438.0   2512.5  -1206.0   2412.0     2269 
# 
# AIC      BIC   logLik deviance df.resid 
# 2440.0   2520.2  -1206.0   2412.0     2268 

pop_pred_inter=ggpredict(surv_pop_inter, terms = c("FROH[all]","AnnualDensity[0,0.001,0.002,0.003,0.004, 0.005]"))
pop_pred_inter

pop_pred_inter_plot=ggplot(pop_pred_inter, aes(x=x, y=predicted, color=group))+
  geom_line(linewidth=1)+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, color=NULL, fill=group),alpha=0.15, show.legend = F)+
  theme_bw()+
  scale_fill_manual(values=cols_5)+
  scale_color_manual(values = cols_5)+
  labs(x=expression(F["ROH"]), y="Juvenile survival probability",color="Population \ndensity")+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=12),
        legend.text=element_text(size=rel(0.5)),
        axis.title.y = element_text(size=13))
pop_pred_inter_plot



### annual den vs grazing ####


## first model den on grazing 

mod=lm(AnnualDensity~GrazeType, data=loc_den_graze)
summary(mod)

plot(AnnualDensity~GrazeType, data=loc_den_graze)

standard_res=rstandard(mod)

data_with_res=cbind(loc_den_graze,standard_res)%>%
  mutate(Reg = fct_relevel(Reg, "SI", "IM", "LA", "NG","MG",  "SG"))%>%
  mutate(graz_v_den=AnnualDensity/GrazeType)


#check how competition differs between regions 
comp_reg=ggplot(data_with_res,aes(x=Reg, y=standard_res, colour=Reg))+
  geom_boxplot()+
  theme_bw()+
  scale_color_manual(values = c("#f0a30a" ,"#a20025","#00aba9","chocolate1", "#60a917","#647687"))+
  labs(x="Spatial region", y="'Competition'")+
  theme(text = element_text(size = 18),legend.position = "none")

comp_reg

##run model of residuals on survival
surv_pop_graze_inter=glmmTMB(juvenile_survival~ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+
                         (FROH*standard_res)+
                         (1|BirthYear)+(1|MumCode), 
                       family=binomial, 
                       data=data_with_res, 
                       na.action = na.omit,
)

summary(surv_pop_graze_inter)

#no inter
# AIC      BIC   logLik deviance df.resid 
# 2469.9   2544.5  -1222.0   2443.9     2269 
# with inter

# AIC      BIC   logLik deviance df.resid 
# 2471.6   2551.9  -1221.8   2443.6     2268 

pop_graze_pred_inter=ggpredict(surv_pop_graze_inter, terms = c("FROH[all]","standard_res[all]"))
pop_graze_pred_inter

graze_pop_pred_inter_plot=ggplot(pop_graze_pred_inter, aes(x=x, y=predicted, color=group))+
  geom_line(linewidth=1)+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, color=NULL, fill=group),alpha=0.15, show.legend = F)+
  theme_bw()+
  scale_fill_manual(values=cols_5)+
  scale_color_manual(values = cols_5)+
  labs(x=expression(F["ROH"]), y="Juvenile survival probability",color="'Competition'")+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=12),
        legend.text=element_text(size=rel(0.5)),
        axis.title.y = element_text(size=13))
graze_pop_pred_inter_plot


pop_graz=(graze_pred_inter_plot/pop_pred_inter_plot/graze_pop_pred_inter_plot)+
  plot_annotation(tag_levels = "A")
pop_graz

ggsave(pop_graz,
       file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/juve_graze_pop_inters.jpeg",
       width = 6,
       height = 10)

