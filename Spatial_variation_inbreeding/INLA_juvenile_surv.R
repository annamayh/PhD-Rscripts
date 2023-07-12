library(INLA)
library(inlabru)
library(tidyverse)
library(ggregplot)

setwd("H:/")

surv_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)
  

surv_loc_df_study_area=surv_loc_df%>%
  filter(!E>1385)%>%
  filter(!N<7997.5) %>%#removing ids with no known region or ~10 ids with outside the limits of study area
  dplyr::select(-MumFROH)%>%
  na.omit()

table(surv_loc_df_study_area$MotherStatus)
surv_loc_df_study_area$MotherStatus[surv_loc_df_study_area$MotherStatus == "Naïve"] <- "Naive"


surv_loc_df_study_area$Sex=as.factor(surv_loc_df_study_area$Sex)#
surv_loc_df_study_area$Code=as.factor(surv_loc_df_study_area$Code)#
surv_loc_df_study_area$BirthYear=as.factor(surv_loc_df_study_area$BirthYear)#
surv_loc_df_study_area$MumCode=as.factor(surv_loc_df_study_area$MumCode)#
surv_loc_df_study_area$MotherStatus=as.factor(surv_loc_df_study_area$MotherStatus)#

head(surv_loc_df_study_area)

table(surv_loc_df$Reg)
table(surv_loc_df_study_area$Reg)



###########################################################################################
## model selection using INLA
############################################################################

IM1.1  <- inla(juvenile_survival~Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+BirthWt, 
               family = "binomial",
               data = surv_loc_df_study_area,
               control.compute = list(dic=TRUE)) 


summary(IM1.1)



IM1  <- inla(juvenile_survival~Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+BirthWt+
               f(BirthYear, model = 'iid')+f(MumCode, model = 'iid'), 
             family = "binomial",
             data = surv_loc_df_study_area,
             control.compute = list(dic=TRUE)) 


summary(IM1)




IM2  <- inla(juvenile_survival~1+Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH+BirthWt+
  f(BirthYear, model = 'iid')+f(MumCode, model = 'iid'), 
             family = "binomial",
             data = surv_loc_df_study_area,
             control.compute = list(dic=TRUE)) 


summary(IM2)



modelsel <- list(IM1.1,IM1, IM2)
sapply(modelsel, function(f) f$dic$dic)
INLADICFig(modelsel, ModelNames = c("IM1.1","IM1" ,"IM2"))




##### now add in spatial variation #####
rum_outline=read.csv("PhD_4th_yr/Spatial_var_inbreeding/INLA/RumBoundary.csv")%>%
  rename(E = Easting, N = Northing) 

N=nrow(rum_outline)
rum_line_rev=rum_outline[N:1, c("E","N")]

Mesh=inla.mesh.2d(loc.domain= rum_outline, 
                  max.edge=2, #probs use 1 for actual model
                  boundary=
                    inla.mesh.segment(rum_line_rev))

plot(Mesh, asp=1)
### plotting visual of mesh 
Rum_mesh=ggplot()+
  gg(Mesh)+
  geom_point(aes(surv_loc_df_study_area$E, surv_loc_df_study_area$N, colour=surv_loc_df_study_area$Reg), alpha=0.5)+
  scale_colour_manual(values = c("#a20025","#00aba9","#60a917","chocolate1", "#647687","#f0a30a" ))+
  theme_classic()+
  labs(x="Easting", y="Northing", colour="Region")+
  theme(text = element_text(size = 18))
Rum_mesh

# ggsave(Rum_mesh,
#        file = "PhD_4th_yr/Spatial_var_inbreeding/Chapter_wrting/plots/Rum_mesh.png",
#        width = 7,
#        height = 8, 
#        dpi=300)
# 


Locations=cbind(surv_loc_df_study_area$E, surv_loc_df_study_area$N)#locations of ids
#make A mtrix
A=inla.spde.make.A(Mesh, loc=Locations)
dim(A)

#define SPDE 
spde=inla.spde2.matern(Mesh, alpha = 2) # would need to adjust for time series data

#define spatial field
w.index=inla.spde.make.index(name = 'w', n.spde=spde$n.spde, n.group = 1, n.repl = 1)

#make model matrix
N <- nrow(surv_loc_df_study_area)
X0=data.frame(Intercept = rep(1, N),
              FROH = surv_loc_df_study_area$FROH, 
             Sex = surv_loc_df_study_area$Sex, 
             MotherStatus = surv_loc_df_study_area$MotherStatus, 
             mum_age =surv_loc_df_study_area$mum_age,
             mum_age_sq = surv_loc_df_study_area$mum_age_sq, 
             Day_seq=surv_loc_df_study_area$Day_seq, 
             BirthWt=surv_loc_df_study_area$BirthWt)

x=as.data.frame(X0)

#make stack
stackfit=inla.stack(
  tag="Fit",
  data=list(y=surv_loc_df_study_area$juvenile_survival), 
  A = list(A, 1, 1, 1), #for the 6 fixed and random effects I have 
  effects=list(
    w=w.index,
    X=x,
    BirthYear = surv_loc_df_study_area$BirthYear, 
    MumCode = surv_loc_df_study_area$MumCode)
)


IM_spde  <- inla(y~ -1 + Intercept+Sex + MotherStatus + mum_age+mum_age_sq+Day_seq+FROH+BirthWt+
                   f(BirthYear, model = 'iid')+f(MumCode, model = 'iid')+ f(w, model=spde), 
                 family = "binomial",
                 data=inla.stack.data(stackfit), 
                 control.compute = list(dic=TRUE),
                 control.predictor = list(
                   A=inla.stack.A(stackfit))) 


summary(IM_spde)

SpatialList <- list(IM1, IM1.1, IM2, IM_spde)
sapply(SpatialList, function(f) f$dic$dic)
INLADICFig(SpatialList, ModelNames = c("IM1", "IM1.1", "IM2" ,"SPDE_1"))+theme_classic()


inla_surv_plot=ggField(IM_spde, Mesh)+
  labs(tag="D", fill = "Juvenile survival \n(as untransformed \ndevaition from mean)")+
  theme_bw()+
  scale_fill_discrete_sequential(palette = "Oranges", rev=FALSE)+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=rel(0.8)))


inla_surv_plot







##############################################################################################
#### now try to spliot up the FROH values into groups 
#####################################################################################################

# 
# 
# surv_loc_df_study_area_gr=surv_loc_df_study_area%>%
#   mutate(Groups=case_when(FROH<=0.05~1,
#                          FROH>0.05&FROH<=0.1~2,
#                          FROH>0.1&FROH<=0.15~3, 
#                          FROH>0.15~4
#                          ))
# 
# 
# NGroups <- length(unique(surv_loc_df_study_area_gr$Groups)) 
# 
# 
# A2=inla.spde.make.A(Mesh, loc=Locations,
#                    group = as.numeric(as.factor(surv_loc_df_study_area_gr$Groups)), 
#                    n.group = NGroups)
# dim(A2)
# 
# #define spatial field
# w.index=inla.spde.make.index(name = 'w', n.spde=spde$n.spde, n.group = NGroups)
# 
# #make model matrix
# 
# N <- nrow(surv_loc_df_study_area)
# X0=data.frame(Intercept = rep(1, N),
#               Sex = surv_loc_df_study_area$Sex, 
#               MotherStatus = surv_loc_df_study_area$MotherStatus, 
#               mum_age =surv_loc_df_study_area$mum_age,
#               mum_age_sq = surv_loc_df_study_area$mum_age_sq)
# 
# x=as.data.frame(X0)
# 
# 
# 
# #make stack
# stackfit2=inla.stack(
#   tag="Fit",
#   data=list(y=surv_loc_df_study_area$juvenile_survival), 
#   A = list(A2, 1, 1, 1), #for the 6 fixed and random effects I have 
#   effects=list(
#     w=w.index,
#     X=x,
#     BirthYear = surv_loc_df_study_area$BirthYear, 
#     MumCode = surv_loc_df_study_area$MumCode)
# )
# 
# inla.setOption(num.threads = 8) 
# 
# IM_spde_gr  <- inla(y~ -1 + Intercept+Sex + MotherStatus + mum_age+mum_age_sq+
#                    f(BirthYear, model = 'iid')+f(MumCode, model = 'iid')+ 
#                    f(w, model=spde, group = w.group,control.group = list(model = 'iid')), 
#                  family = "binomial",
#                  data=inla.stack.data(stackfit2), 
#                  control.compute = list(dic=TRUE),
#                  control.predictor = list(
#                    A=inla.stack.A(stackfit2))) 
# 
# 
# summary(IM_spde_gr)
# 
# IM_spde_gr$summary.random$w
# 
# Labels = c("Low FROH", "Med-low FROH", "Med-high FROH", "High FROH")
# names(Labels) <- c(1:NGroups)
# 
# ggField(IM_spde_gr, Mesh, Groups = NGroups) + # Notice the groups argument, using the number of unique months.
#   scale_fill_brewer(palette = "Reds") + 
#   facet_wrap( ~ Groups, labeller = labeller(Groups = Labels), ncol = 3) # Doing this manually changes the facet labels
# 
