

library(INLA)
library(inlabru)
library(tidyverse)
library(ggregplot)

setwd("H:/")

spatial_FROH=read.table("PhD_4th_yr/Spatial_var_inbreeding/INLA/Location_plus_FROH.txt", sep = ",", header = TRUE)%>%
  filter(!E>1385)%>%
  filter(!N<7997.5)%>% #remove 38 records outside the syudy areas
  rename(X = E, Y = N) %>% #rename easting and northing to X and Y
  mutate(year_cont=BirthYear-min(BirthYear))#%>% #get cntinuous vairable of birth year
 # mutate(FROH=FROH*10)

###########################################################################################
## model selection using INLA
############################################################################

IM1  <- inla(FROH~year_cont, 
             family = "gaussian",
             data = spatial_FROH,
             control.compute = list(dic=TRUE)) 


summary(IM1)


spatial_FROH$BirthYear=as.factor(spatial_FROH$BirthYear)#treating birth yr as a factor

IM2  <- inla(FROH~year_cont+f(BirthYear, model = 'iid'), 
             family = "gaussian",
             data = spatial_FROH,
             control.compute = list(dic=TRUE)) 


summary(IM2)
## looks like year as a random variable explains some variation but not year as a continuous: makes sense 


##############################################################################################
setwd("H:/")

rum_outline=read.csv("PhD_4th_yr/Spatial_var_inbreeding/INLA/RumBoundary.csv")%>%
  rename(X = Easting, Y = Northing) 
  


N=nrow(rum_outline)
rum_line_rev=rum_outline[N:1, c("X","Y")]

Mesh=inla.mesh.2d(loc.domain= rum_outline, 
                  max.edge=2, #probs use 1 for actual model
                  boundary=
                  inla.mesh.segment(rum_line_rev))

plot(Mesh, asp=1)

LocToReg6 <- function(E, N) {
  ifelse(N < 8019, "SG", #south glen
         ifelse(E < 1361, "LA", # laundry greens
                ifelse(E < 1366,
                       ifelse(N > 8033, "NG", "MG"), # north glen, mid  glen 
                       ifelse(E < 1373 , "IM", "SI")))) }


spatial_FROH$Reg <- with(spatial_FROH, LocToReg6(X, Y))

### plotting visual of mesh 
Rum_mesh=ggplot()+
  gg(Mesh)+
  geom_point(aes(spatial_FROH$X, spatial_FROH$Y, colour=spatial_FROH$Reg), alpha=0.4)+
  theme_classic()+
  labs(x="Easting", y="Northing")
Rum_mesh

#################################################################################################
## setting up INLA model ########################################################################
######################################################################################################

#set weighting using A matrix

Locations=cbind(spatial_FROH$X, spatial_FROH$Y)#locations of ids


A=inla.spde.make.A(Mesh, loc=Locations)


#define SPDE 

spde=inla.spde2.matern(Mesh, alpha = 2) # would need to adjust for time series data

#define spatial field

w.index=inla.spde.make.index(name = 'w', n.spde=spde$n.spde, n.group = 1, n.repl = 1)

#make model matrix
Xm=model.matrix(~-1 + year_cont, data = spatial_FROH)
X=data.frame(year_cont=Xm[,1])
head(X)

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

Fspat=as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "), 
                        " + f(w, model=spde)"))

#run model
IM_sp=inla(Fspat,
  family = "gaussian", 
  data=inla.stack.data(stackfit), 
  control.compute = list(dic=TRUE),
  control.predictor = list(
    A=inla.stack.A(stackfit))#, verbose = TRUE
  
)

summary(IM_sp)


## now fitting spde and birth year as random effects 
# have to re-do the stack
stackfit2=inla.stack(
  data=list(y=spatial_FROH$FROH), 
  A = list(1, 1, A), 
  effects=list(
    Intercept=rep(1, N), 
    BirthYear = spatial_FROH$BirthYear, # insert vectors of any random effects
    w=w.index)
)


Fspat2=as.formula(paste0("y ~ -1 + Intercept + f(BirthYear, model = 'iid') + f(w, model=spde) "))

#run model
IM_sp2=inla(Fspat2,
           family = "gaussian", 
           data=inla.stack.data(stackfit2), 
           control.compute = list(dic=TRUE),
           control.predictor = list(
             A=inla.stack.A(stackfit))#, verbose = TRUE
           
)

summary(IM_sp2)




## check whether DIC decreases with SPDE 
SpatialList <- list(IM1, IM2, IM_sp, IM_sp2)
sapply(SpatialList, function(f) f$dic$dic)
INLADICFig(SpatialList, ModelNames = c("IM1","IM2" ,"SPDE_1","SPDE_2"))


ggField(IM_sp2, Mesh, Fill = "Continuous")

