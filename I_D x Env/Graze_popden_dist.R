library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)
library(INLA)
library(ggregplot)
library(inlabru)
library(colorspace)


setwd("H:/")

surv_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)%>%
  select(Code, MumCode, BirthYear,Sex, N, E, Reg)%>%
  filter(Sex!=3)%>%
  na.omit()

density=read.csv("PhD_4th_yr/Spatial_var_inbreeding/glmmTMB/ValuesforAnna.csv", header=T, stringsAsFactors = F)%>%
  rename(MumCode=Name, BirthYear=Year)

loc_den_graze=surv_loc_df%>%
  inner_join(density)%>%
  unique()%>%
  na.omit()#%>%
  #mutate(AnnualDensity=AnnualDensity*1000) #rememnber is in e-03 now

loc_den_graze$Code=as.factor(loc_den_graze$Code)
loc_den_graze$Sex=as.factor(loc_den_graze$Sex)
loc_den_graze$Reg=as.factor(loc_den_graze$Reg)



cor(loc_den_graze$AnnualDensity, loc_den_graze$GrazeType)

cor.test(formula=~GrazeType+AnnualDensity, data=loc_den_graze)

cor=ggplot(loc_den_graze, aes(y=AnnualDensity, x=GrazeType))+
  geom_point()+
  geom_smooth(method=lm)+
  theme_classic()+
  #labs(x="Density", y="Graze quality", subtitle = "r=0.58")+
  theme(text = element_text(size = 18))
        
cor 

graze_simple=glmmTMB(GrazeType~ Sex + Reg, 
                         family="gaussian", 
                         data=loc_den_graze, 
                         na.action = na.omit,
)


summary(graze_simple)

graze_pred=ggpredict(graze_simple, terms = c("Reg"))%>%
  mutate(x = fct_relevel(x, "SI", "IM", "LA", "NG","MG",  "SG"))
  

loc_den_graze_plot=loc_den_graze%>%rename(x=Reg, predicted=GrazeType)

graz_pred_plot=ggplot(graze_pred,aes(x=x, y=predicted, colour=x, ymin=conf.low, ymax=conf.high, group=group))+
  geom_pointrange(linewidth=1)+
  theme_bw()+
  scale_color_manual(values = c("#f0a30a" ,"#a20025","#00aba9","chocolate1", "#60a917","#647687"))+
  labs(x="Spatial region", y="Grazing quality")+
  theme(text = element_text(size = 18),legend.position = "none")+
  geom_boxplot(data=loc_den_graze, aes(Reg, GrazeType, colour=Reg, alpha=0.1),
               inherit.aes = F, position=position_nudge(x=0.2), width=0.1)


graz_pred_plot


### graze unsing inla



IM1  <- inla(GrazeType~ Sex + Reg, 

             family = "gaussian",
             data = loc_den_graze,
             control.compute = list(dic=TRUE)) 


summary(IM1)


## add spatial field

setwd("H:/")

rum_outline=read.csv("PhD_4th_yr/Spatial_var_inbreeding/INLA/RumBoundary.csv")%>%
  rename(E = Easting, N = Northing) 


N=nrow(rum_outline)
rum_line_rev=rum_outline[N:1, c("E","N")]

Mesh=inla.mesh.2d(loc.domain= rum_outline, 
                  max.edge=2, #probs use 1 for actual model
                  boundary=
                    inla.mesh.segment(rum_line_rev))

plot(Mesh, asp=1)

Locations=cbind(loc_den_graze$E, loc_den_graze$N)#locations of ids

A=inla.spde.make.A(Mesh, loc=Locations)

#define SPDE 
spde=inla.spde2.matern(Mesh, alpha = 2) # would need to adjust for time series data

#define spatial field
w.index=inla.spde.make.index(name = 'w', n.spde=spde$n.spde, n.group = 1, n.repl = 1)

#make model matrix
N <- nrow(loc_den_graze)
X0=data.frame(Intercept = rep(1, N),
              Sex = loc_den_graze$Sex, 
              Reg= loc_den_graze$Reg)

x=as.data.frame(X0)


## now fitting spde and birth year as random effects 
# have to re-do the stack
stackfit2=inla.stack(
  data=list(y=loc_den_graze$GrazeType), 
  A = list(1,  A), 
  effects=list(
    X=x,
    w=w.index)
)

#define SPDE 


IM_spde_graze  <- inla(y~ -1 + Intercept+Sex + Reg+ f(w, model=spde), 
                 family = "gaussian",
                 data=inla.stack.data(stackfit2), 
                 control.compute = list(dic=TRUE),
                 control.predictor = list(
                   A=inla.stack.A(stackfit2))) 


summary(IM_spde_graze)

inla_graze_plot=ggField(IM_spde_graze, Mesh)+
  labs(fill = "Grazing Quality")+
  theme_bw()+
  scale_fill_discrete_sequential(palette = "Greens", rev=T, 
                                 labels=c("Low", " "," "," "," "," "," "," ","High"))+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=rel(0.8)))


inla_graze_plot

SpatialList <- list(IM1,  IM_spde_graze)
sapply(SpatialList, function(f) f$dic$dic)
INLADICFig(SpatialList)+theme_classic()











######################
## Annual Denisty ####
######################



den_simple=glmmTMB(AnnualDensity~ Sex + Reg, 
                     family="gaussian", 
                     data=loc_den_graze, 
                     na.action = na.omit,
)


summary(den_simple)

den_pred=ggpredict(den_simple, terms = c("Reg"))%>%
  mutate(x = fct_relevel(x, "SI", "IM", "LA", "NG","MG",  "SG"))



den_pred_plot=den_pred%>%
  ggplot(aes(x=x, y=predicted, colour=x, ymin=conf.low, ymax=conf.high, group=group))+
  geom_pointrange(linewidth=1)+
  theme_bw()+
  scale_color_manual(values = c("#f0a30a" ,"#a20025","#00aba9","chocolate1", "#60a917","#647687"))+
  labs(x="Spatial region", y="Derived population denisty metric")+
  theme(text = element_text(size = 18),legend.position = "none")+
  geom_boxplot(data=loc_den_graze, aes(Reg, AnnualDensity, colour=Reg, alpha=0.1),
               inherit.aes = F, position=position_nudge(x=0.2), width=0.1)


den_pred_plot


graz_pred_plot+den_pred_plot


##unsing INLA ##

loc_den_graze_2=loc_den_graze%>%mutate(AnnualDensity=AnnualDensity*1000) #rememnber is in e-03 now


IM1_den  <- inla(AnnualDensity~ Sex + Reg, 
             family = "gaussian",
             data = loc_den_graze_2,
             control.compute = list(dic=TRUE)) 


summary(IM1_den)


stackfit_den=inla.stack(
  data=list(y=loc_den_graze_2$AnnualDensity), 
  A = list(1,  A), 
  effects=list(
    X=x,
    w=w.index)
)

#define SPDE 


IM_spde_den  <- inla(y~ -1 + Intercept+Sex + Reg+ f(w, model=spde), 
                       family = "gaussian",
                       data=inla.stack.data(stackfit_den), 
                       control.compute = list(dic=TRUE),
                       control.predictor = list(
                         A=inla.stack.A(stackfit_den))) 


summary(IM_spde_den)

inla_den_plot=ggField(IM_spde_den, Mesh)+
  labs(fill = "Derived \npopulation \ndensity")+
  theme_bw()+
  scale_fill_discrete_sequential(labels=c("Low", " "," "," "," "," "," "," ","High"), 
                                 palette = "Purples", rev=T)+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=rel(0.8)))


inla_den_plot

SpatialList <- list(IM1_den,  IM_spde_den)
sapply(SpatialList, function(f) f$dic$dic)
INLADICFig(SpatialList)+theme_classic()





(graz_pred_plot+inla_graze_plot)/(den_pred_plot+inla_den_plot)


#save(graz_pred_plot,inla_graze_plot,den_pred_plot,inla_den_plot, file="PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/Graze_density_plots.RData")


#load(file="PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/Graze_density_plots.RData")

graz_den=((graz_pred_plot+inla_graze_plot)/(den_pred_plot+inla_den_plot))+plot_annotation(tag_levels = "A")
graz_den

ggsave(graz_den,
       file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/graz_den_all.png",
       width = 10,
       height = 10, 
       dpi=2000)


