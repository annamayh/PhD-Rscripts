## brms post ###
library(dplyr)
library(purrr)
library(tidybayes)
library(ggplot2)

setwd("M:/")

load("PhD_3rdYR/Model outputs/m2_birthwt.RData")


summary(m2)

get_variables(m2)

plot(m2)
