

library(INLA)
library(inlabru)
library(tidyverse)
library(ggregplot)

setwd("H:/")

rum_outline=read.csv("PhD_4th_yr/Spatial_var_inbreeding/INLA/RumBoundary.csv")


N=nrow(rum_outline)
rum_line_rev=rum_outline[N:1, c("Easting","Northing")]

mesh=inla.mesh.2d(loc.domain= rum_outline, 
                  max.edge=2, #probs use 1 for actual model
                  boundary=
                  inla.mesh.segment(rum_line_rev))

plot(mesh, asp=1)

spatial_FROH=read.table("PhD_4th_yr/Spatial_var_inbreeding/INLA/Location_plus_FROH.txt", sep = ",", header = TRUE)%>%
  filter(!E>1385)%>%
  filter(!N<7997.5) #remove 38 records outside the syudy areas


LocToReg6 <- function(E, N) {
  ifelse(N < 8019, "SG", #south glen
         ifelse(E < 1361, "LA", # laundry greens
                ifelse(E < 1366,
                       ifelse(N > 8033, "NG", "MG"), # north glen, mid  glen 
                       ifelse(E < 1373 , "IM", "SI")))) }


spatial_FROH$Reg <- with(spatial_FROH, LocToReg6(E, N))

### plotting visual of mesh 
(Rum_mesh=ggplot()+
  gg(mesh)+
  geom_point(aes(spatial_FROH$E, spatial_FROH$N, colour=spatial_FROH$Reg), alpha=0.4)+
  theme_classic()+
  labs(x="Easting", y="Northing"))



###########################################################################################
## model selection using INLA
############################################################################

resp_var="FROH" #defining response variable

covar=c("BirthYear") #define covariates for model

## specify basic formula with birth year as random effect
f1 <- as.formula(paste0(resp_var, " ~ ", 
                        paste(covar, collapse = " + "))) 
IM1  <- inla(f1, 
               family = "gaussian",
               data = spatial_FROH) 


summary(IM1)


#################################################################################################
## setting up INLA model ########################################################################
######################################################################################################

#set weighting using A matrix

Locations=cbind(spatial_FROH$E, spatial_FROH$N)#locations of ids


A=inla.spde.make.A(mesh, loc=Locations)


#define SPDE 

spde=inla.spde2.matern(mesh, alpha = 2) # would need to adjust for time series data

#define spatial field

w.index=inla.spde.make.index(name = 'w', n.spde=spde$n.spde, n.group = 1, n.repl = 1)

#make model matrix
Xm=model.matrix(~-1 + BirthYear, data = spatial_FROH)
X=data.frame(BirthYear=Xm[,1])


N <- nrow(spatial_FROH)

#make stack
stackfit=inla.stack(
  data=list(y=spatial_FROH$FROH), 
  A = list(1, 1, A), 
  effects=list(
    Intercept=rep(1, N), 
    X=X,
    w=w.index)
  )



# define model formula

Fspat=as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "), "+ f(w, model=spde)"))

#run model
IM_sp=inla(Fspat,
  family = "gaussian", 
  data=inla.stack.data(stackfit), 
  control.compute = list(dic=TRUE),
  control.predictor = list(
    A=inla.stack.A(stackfit))#, verbose = TRUE
  
)

summary(IM_sp)


SpatialList <- list(IM1, IM_sp)

INLADICFig(SpatialList, ModelNames = c("Base","SPDE"))


ggField(IM_sp, mesh)

IM_sp$


IM_sp$Spatial$Model %>%
  ggField(IM_sp$Spatial$Mesh, Points = IM_sp$Data[,c("X", "Y")],
          PointAlpha = 0.1) +
  scale_fill_discrete_sequential(palette = AlberPalettes[[1]])
