

library(INLA)
library(inlabru)
library(tidyverse)
library(ggregplot)
library(colorspace)
library(RColorBrewer)

setwd("H:/")

spatial_FROH=read.table("PhD_4th_yr/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)%>%
  filter(!E>1385)%>%
  filter(!N<7997.5)%>% #remove 38 records outside the syudy areas
  mutate(year_cont=BirthYear-min(BirthYear))%>% #get cntinuous vairable of birth year
  select(Code, BirthYear, Sex, MumCode, FROH, Reg, juvenile_survival)%>%
  na.omit()
###########################################################################################
## model selection using INLA
############################################################################

IM1  <- inla(FROH~year_cont, 
             family = "gaussian",
             data = spatial_FROH,
             control.compute = list(dic=TRUE)) 


summary(IM1)

spatial_FROH$MumCode=as.factor(spatial_FROH$MumCode)


IM2  <- inla(FROH~year_cont+ f(MumCode, model = 'iid'), 
             family = "gaussian",
             data = spatial_FROH,
             control.compute = list(dic=TRUE)) 


summary(IM2)
List <- list(IM1, IM2)
sapply(List, function(f) f$dic$dic)
INLADICFig(List)



##############################################################################################
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


#################################################################################################
## setting up INLA model ########################################################################
######################################################################################################

#set weighting using A matrix

Locations=cbind(spatial_FROH$E, spatial_FROH$N)#locations of ids

A=inla.spde.make.A(Mesh, loc=Locations)

#define SPDE 
spde=inla.spde2.matern(Mesh, alpha = 2) # would need to adjust for time series data

#define spatial field
w.index=inla.spde.make.index(name = 'w', n.spde=spde$n.spde, n.group = 1, n.repl = 1)

#make model matrix
N <- nrow(spatial_FROH)
X0=data.frame(Intercept = rep(1, N),
              year_cont = spatial_FROH$year_cont)

x=as.data.frame(X0)


## now fitting spde and birth year as random effects 
# have to re-do the stack
stackfit2=inla.stack(
  data=list(y=spatial_FROH$FROH), 
  A = list(1, 1, A), 
  effects=list(
    X=x,
    MumCode = spatial_FROH$MumCode, # insert vectors of any random effects
    w=w.index)
)


Fspat2=as.formula(paste0("y ~ -1 + Intercept + year_cont + f(MumCode, model = 'iid') + f(w, model=spde) "))

#run model
IM_sp2=inla(y ~ -1 + Intercept + year_cont + 
              f(MumCode, model = 'iid') + f(w, model=spde),
           family = "gaussian", 
           data=inla.stack.data(stackfit2), 
           control.compute = list(dic=TRUE),
           control.predictor = list(
             A=inla.stack.A(stackfit2))#, verbose = TRUE
           
)

summary(IM_sp2)




## check whether DIC decreases with SPDE 
SpatialList <- list(IM1, IM2, IM_sp2)
sapply(SpatialList, function(f) f$dic$dic)
INLADICFig(SpatialList, ModelNames = c("IM1","IM2" ,"SPDE_1"))

library(viridisLite)
coul <- viridis(100)


inla_froh_gg=ggField(IM_sp2, Mesh, Fill="Continuous")+
  labs(fill = "FROH\n(as deviation \nfrom mean)")+
  theme_bw()+
  scale_fill_continuous_sequential(palette = "BluYl")+
  theme(text = element_text(size = 18),
        legend.title=element_text(size=rel(0.8)))

inla_froh_gg



