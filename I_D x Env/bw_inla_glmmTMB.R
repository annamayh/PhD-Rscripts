
library(RODBC)
library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(emmeans)

surv_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)%>%
  select(-BirthWt)


db<-"C:\\Users\\s1881212\\Documents\\Deer_database_2022/RedDeer2.05.accdb" #open connection
con<-odbcConnectAccess2007(db)
birth_wt<-sqlFetch(con, "sys_BirthWt") 

odbcClose(con) #close connection


bw_df=surv_loc_df%>%full_join(birth_wt)%>%na.omit
head(bw_df)

bw_df$Code=as.factor(bw_df$Code)
bw_df$MumCode=as.factor(bw_df$MumCode)
bw_df$BirthYear=as.factor(bw_df$BirthYear)
bw_df$Sex=as.factor(bw_df$Sex)
bw_df$Reg=as.factor(bw_df$Reg)



bw_model_simple=glmmTMB(CaptureWt~ Sex + AgeHrs+ MotherStatus+
                            (1|BirthYear)+ (1|MumCode), 
                          family=gaussian(), 
                          data=bw_df, 
                          na.action = na.omit,
)

summary(bw_model_simple)


bw_reg_fixed=update(bw_model_simple, ~ . + Reg+FROH) ##just region as fixed effect
summary(bw_reg_fixed)

bw_pred=ggpredict(bw_reg_fixed, terms = c("Reg"))
bw_pred

bw_reg=bw_pred%>%
  mutate(x = fct_relevel(x, "SI", "IM", "LA", "NG","MG",  "SG"))%>%
  ggplot(aes(x=x, y=predicted, color=x, ymin=conf.low, ymax=conf.high))+
  geom_pointrange(linewidth=1)+
  theme_bw()+
  scale_color_manual(values = c("#f0a30a" ,"#a20025","#00aba9","chocolate1","#60a917", "#647687"))+
  labs(x="Spatial region", y="Birth Weight (kg)", tag="E")+
  theme(text = element_text(size = 18),legend.position = "none")
bw_reg


bw_reg_inter=update(bw_model_simple, ~ . + Reg:FROH) ##just region as fixed effect
summary(bw_reg_inter)

bw_inter=plot(ggpredict(bw_reg_inter, terms = c("FROH[all]","Reg")),show.title=FALSE, colors="metro", line.size=1)+
  labs(x = expression(F["ROH"]), y = "Birth weight (kg)", colour = "Spatial region")+
  theme(text = element_text(size = 15)) 
bw_inter

FROH_trend=emtrends(bw_reg_inter, pairwise ~ Reg, var="FROH")
test(FROH_trend$emtrends)






### now same using INLA


IM1  <- inla(CaptureWt~ Sex + AgeHrs+ MotherStatus+
               f(BirthYear, model = 'iid')+f(MumCode, model = 'iid'), 
               
             family = "gaussian",
             data = bw_df,
             control.compute = list(dic=TRUE)) 


summary(IM1)

IM2  <- inla(CaptureWt~ Sex + AgeHrs+ MotherStatus+FROH+
               f(BirthYear, model = 'iid')+f(MumCode, model = 'iid'), 
             
             family = "gaussian",
             data = bw_df,
             control.compute = list(dic=TRUE)) 


summary(IM2)

List <- list(IM1, IM2)
sapply(List, function(f) f$dic$dic)
INLADICFig(List)

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

Locations=cbind(bw_df$E, bw_df$N)#locations of ids

A=inla.spde.make.A(Mesh, loc=Locations)

#define SPDE 
spde=inla.spde2.matern(Mesh, alpha = 2) # would need to adjust for time series data

#define spatial field
w.index=inla.spde.make.index(name = 'w', n.spde=spde$n.spde, n.group = 1, n.repl = 1)

#make model matrix
N <- nrow(bw_df)
X0=data.frame(Intercept = rep(1, N),
              Sex = bw_df$Sex, 
              MotherStatus=bw_df$MotherStatus, 
              AgeHrs=bw_df$AgeHrs, 
              FROH=bw_df$FROH)

x=as.data.frame(X0)


## now fitting spde and birth year as random effects 
# have to re-do the stack
stackfit2=inla.stack(
  data=list(y=bw_df$CaptureWt), 
  A = list(1, 1, 1, A), 
  effects=list(
    X=x,
    BirthYear = bw_df$BirthYear,
    MumCode = bw_df$MumCode, # insert vectors of any random effects
    w=w.index)
)

#define SPDE 


IM_spde  <- inla(y~ -1 + Intercept+Sex + MotherStatus + AgeHrs+FROH+
                   f(BirthYear, model = 'iid')+f(MumCode, model = 'iid')+ f(w, model=spde), 
                 family = "gaussian",
                 data=inla.stack.data(stackfit2), 
                 control.compute = list(dic=TRUE),
                 control.predictor = list(
                   A=inla.stack.A(stackfit2))) 


summary(IM_spde)


List <- list(IM1, IM2, IM_spde)
sapply(List, function(f) f$dic$dic)
INLADICFig(List)



inla_bw_plot=ggField(IM_spde, Mesh)+
  labs(tag="F", fill = "Birth weight (kg) \n(as deviation \nfrom mean)")+
  theme_bw()+
  scale_fill_discrete_sequential(palette = "Greens", rev=FALSE)+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=rel(0.8)))


inla_bw_plot





bw=bw_reg+inla_bw_plot
bw



# ggsave(bw,
#        file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/bw_both.png",
#        width = 11,
#        height = 7, 
#        dpi=1000)


big=inla_and_glmm/inla_and_reg/bw
big

ggsave(big,
       file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/all.png",
       width = 11,
       height = 18, 
       dpi=1000)

