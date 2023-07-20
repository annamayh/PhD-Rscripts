library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)

setwd("H:/")

surv_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)%>%
  select(-MumFROH)%>%
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

cor(loc_den_graze$GrazeType, loc_den_graze$AnnualDensity)

suv_model_simple=glmmTMB(juvenile_survival~ Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH+
                           (1|BirthYear)+(1|MumCode), 
                         family=binomial, 
                         data=loc_den_graze, 
                         na.action = na.omit,
)

summary(suv_model_simple)


surv_graze_den=update(suv_model_simple, ~ . + GrazeType*AnnualDensity) ##just region as fixed effect
summary(surv_graze_den)
graze_pred=ggpredict(surv_graze_den, terms = c("GrazeType[all]"))
graze_pred
plot(graze_pred)
den_pred=ggpredict(surv_graze_den, terms = c("AnnualDensity[all]","GrazeType[0,0.25,0.5,0.75,1]"))
den_pred
plot(den_pred)



surv_den_inter=update(suv_model_simple, ~ . -FROH+ (FROH*AnnualDensity*GrazeType)) ##just region as fixed effect
summary(surv_den_inter)
den_pred_inter=ggpredict(surv_den_inter, 
                         terms = c("FROH[all]","AnnualDensity[0, 0.00025,0.0005, 0.00075, 0.001, 0.00125,0.0015, 0.00175, 0.002, 0.00225, 0.0025,0.00275, 0.003]"))
den_pred_inter

library(RColorBrewer)
library(viridis)
#cols <- (colorRampPalette(viridis(9, "Viridis"))(20))
cols=viridis(13)

denisty_pred_inter_plot=ggplot(den_pred_inter, aes(x=x, y=predicted, color=group))+
  geom_line(linewidth=1)+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, color=NULL, fill=group),alpha=0.1, show.legend = F)+
  theme_bw()+
  scale_fill_manual(values=cols)+
  scale_color_manual(values = cols)+
  labs(x=expression(F["ROH"]), y="Juvenile survival probability",color="Population \ndensity")+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=rel(0.8)),
        legend.text=element_text(size=rel(0.5)))
denisty_pred_inter_plot

# 
# denisty_pred_inter_plot=ggplot(den_pred_inter, aes(x=x, y=group))+
#   geom_density_2d()+
#   theme_bw()+
#   #labs(x=expression(F["ROH"]), y="Juvenile survival probability",color="Population \ndensity", tag = "A")+
#   theme(text = element_text(size = 18))
# denisty_pred_inter_plot
# 



cols_5 <- viridis(5)
#cols_5


graze_pred_inter=ggpredict(surv_den_inter, terms = c("FROH[all]","GrazeType[0,0.25,0.5,0.75,1]"))
graze_pred_inter

graze_pred_inter_plot=ggplot(graze_pred_inter, aes(x=x, y=predicted, color=group))+
  geom_line(linewidth=1)+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, color=NULL, fill=group),alpha=0.15, show.legend = F)+
  theme_bw()+
  scale_fill_manual(values=cols_5)+
  scale_color_manual(values = cols_5)+
  labs(x=expression(F["ROH"]), y="Juvenile survival probability",color="High quality \ngrazing")+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=rel(0.8)),
        legend.text=element_text(size=rel(0.5)))
graze_pred_inter_plot


den_graze=denisty_pred_inter_plot+graze_pred_inter_plot
den_graze


ggsave(den_graze,
       file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/surv_FROH_inter_graze_population_2.png",
       width = 13,
       height = 7)




#### BIRTH WEIGHT #####
surv_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)%>%
  select(-BirthWt)


db<-"C:\\Users\\s1881212\\Documents\\Deer_database_2022/RedDeer2.05.accdb" #open connection
con<-odbcConnectAccess2007(db)
birth_wt<-sqlFetch(con, "sys_BirthWt") 

odbcClose(con) #close connection


bw_df=surv_loc_df%>%full_join(birth_wt)%>%inner_join(density)%>%
  unique()%>%na.omit()

head(bw_df)

bw_df$Code=as.factor(bw_df$Code)
bw_df$MumCode=as.factor(bw_df$MumCode)
bw_df$BirthYear=as.factor(bw_df$BirthYear)
bw_df$Sex=as.factor(bw_df$Sex)
bw_df$Reg=as.factor(bw_df$Reg)

bw_model_simple=glmmTMB(CaptureWt~ Sex + AgeHrs+ MotherStatus+FROH+
                          (1|BirthYear)+ (1|MumCode), 
                        family=gaussian(), 
                        data=bw_df, 
                        na.action = na.omit,
)

summary(bw_model_simple)



bw_graze_den=update(bw_model_simple, ~ . + GrazeType+AnnualDensity) ##just region as fixed effect
summary(bw_graze_den)
graze_pred=ggpredict(bw_graze_den, terms = c("GrazeType[all]"))
graze_pred
plot(graze_pred)
den_pred=ggpredict(bw_graze_den, terms = c("AnnualDensity[all]"))
den_pred
plot(den_pred)



bw_den_inter=update(bw_model_simple, ~ . -FROH+ FROH*AnnualDensity+ FROH*GrazeType) ##just region as fixed effect
summary(bw_den_inter)
den_pred_inter_bw=ggpredict(bw_den_inter, 
                         terms = c("FROH[all]","AnnualDensity[0, 0.00025,0.0005, 0.00075, 0.001, 0.00125,0.0015, 0.00175, 0.002, 0.00225, 0.0025,0.00275, 0.003]"))
den_pred_inter_bw

denisty_pred_inter_plot_bw=ggplot(den_pred_inter_bw, aes(x=x, y=predicted, color=group))+
  geom_line(linewidth=1)+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, color=NULL, fill=group),alpha=0.2, show.legend = F)+
  theme_bw()+
  scale_fill_manual(values=cols)+
  scale_color_manual(values = cols)+
  labs(x=expression(F["ROH"]), y="Birth weight (kg)",color="Population \ndensity \n(ids / sq km)")+
  theme(text = element_text(size = 18), 
        legend.title=element_text(size=rel(0.8)),
        legend.text=element_text(size=rel(0.5)))
denisty_pred_inter_plot_bw


graze_pred_inter_bw=ggpredict(bw_den_inter, terms = c("FROH[all]","GrazeType[0,0.25,0.5,0.75,1]"))
graze_pred_inter_bw

graze_pred_inter_plot_bw=ggplot(graze_pred_inter_bw, aes(x=x, y=predicted, color=group))+
  geom_line(linewidth=1)+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, color=NULL, fill=group),alpha=0.15, show.legend = F)+
  theme_bw()+
  scale_fill_manual(values=cols_5)+
  scale_color_manual(values = cols_5)+
  labs(x=expression(F["ROH"]), y="Birth weight (kg)",color="High quality \ngrazing")+
  theme(text = element_text(size = 18), 
        legend.title=element_text(size=rel(0.8)),
        legend.text=element_text(size=rel(0.5)))
graze_pred_inter_plot_bw


den_graze_bw=denisty_pred_inter_plot_bw+graze_pred_inter_plot_bw
den_graze_bw


graz_den_both=den_graze/den_graze_bw
graz_den_both=graz_den_both+plot_annotation(tag_levels = 'A')
graz_den_both

ggsave(graz_den_both,
       file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/graze_den_inter_both.png",
       width = 10,
       height = 9)



den=denisty_pred_inter_plot+denisty_pred_inter_plot_bw
den
ggsave(den,
       file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/den_only_inter.png",
       width = 10,
       height = 6)


graz=graze_pred_inter_plot+graze_pred_inter_plot_bw
graz
ggsave(graz,
       file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/SUPP_graze_only_inter.png",
       width = 10,
       height = 6)



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

