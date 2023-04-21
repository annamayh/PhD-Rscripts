
library(ggregplot)
library("INLA")
library(tidyverse)

## much of this is modified from greg Alberys github tutorial: https://ourcodingclub.github.io/tutorials/inla/ 
#and 
##https://datascienceplus.com/spatial-regression-in-r-part-2-inla/

setwd("H:/")

FROH_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/INLA/Location_plus_FROH.txt", sep=",", header=T)%>%
  mutate(year_cont=BirthYear-min(BirthYear)) %>%
  mutate(N_div=N/100)%>%
  mutate(E_div=E/100)

head(FROH_loc_df)


phen <- c("Code", "N", "E") # Base columns with spatial information we'll need

resp <- "FROH" # Response variable

covar <- c("year_cont") # 


FROH_loc_df$MumCode <- as.factor(FROH_loc_df$MumCode)
FROH_loc_df$BirthYear <- as.factor(FROH_loc_df$BirthYear)

# Setting up a custom theme
THEME <- theme(axis.text.x = element_text(size = 12,colour = "black"),
               axis.text.y = element_text(size = 12, colour = "black"),
               axis.title.x = element_text(vjust = -0.35),
               axis.title.y = element_text(vjust = 1.2)) + theme_bw()
#viewing sample regions on Rum coloured by region
(samp_locations <- ggplot(FROH_loc_df, aes(E, N)) + 
    geom_jitter(aes(colour = factor(Reg))) + coord_fixed() + 
    THEME )



# Run a simple model
IM0.1  <- inla(FROH ~ year_cont , 
               family = "gaussian", # Specify the family. Can be a wide range (see r-inla.org).
               data = FROH_loc_df) # Specify the data

summary(IM0.1)



### CREATING MESH OF LOCATIONS #####

Locations = cbind(FROH_loc_df$E_div, FROH_loc_df$N_div) # using the sampling locations 
loc_samp=Locations[sample(nrow(Locations), size=100),]

#using smaller trinagles increases precision but also increse computer power
#from book
mesh0 <- inla.mesh.2d(loc=Locations, max.edge=20)#when only one max-edge specified does not have an outer character domain 
#from greg
Mesh <- inla.mesh.2d(Locations, max.edge = c(20, 50))
plot(Mesh)

Mesh2 <- inla.mesh.2d(loc_samp, max.edge = c(3, 5))
plot(Mesh2)


plot(mesh0)
plot(Mesh)


### MAKING THE PROJECTION MATRIX ###

Amat <- inla.spde.make.A(Mesh, loc = Locations)

### setting priors of SPDE

spde <- inla.spde2.pcmatern(Mesh, 
                            prior.range = c(500, 0.5), # c(500, 0.5) == the probability that the range of the spatial effect is below 500 meters is 0.5.
                            prior.sigma = c(2, 0.05)) #

#For the variation parameter the prior is set in a similar fashion:
##So in the lines above we said: the probability that the variation in the spatial effect is larger than 2 is 0.05.
## probably need to try out different priors






