library(tidyverse)
library(data.table)

setwd("H:/")

load(file="PhD_4th_yr/Inbreeding_depression_models/LBS/male_LBS_model_output_1.5mit.RData")

summary(male_LBS_model)

FROH_sum_sol_pois=summary(male_LBS_model)$solutions[3,1]
FROH_sum_sol_zi=summary(male_LBS_model)$solutions[4,1]



sols_pois<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(matches("FROH_"))%>% ## taking out sols with FROH included
  dplyr::mutate(across(1:33, ~.x + FROH_sum_sol_pois)) ## adding FROHsum to chrFROH values


names <- sols_pois %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_pois,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(sols_pois,2,quantile,probs = c(0.95)) #gets upper confidence interval for all solutions 
CI_lower<-apply(sols_pois,2,quantile,probs = c(0.05)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH

FROH_sols$CHR<-as.factor(FROH_sols$CHR)


ggplot(data=FROH_sols, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  geom_hline(yintercept=FROH_sum_sol_pois, lty=1,colour="red") + ## red line is the average effect of all chromosomes 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="solution + CI ") +
  theme_bw()+  # use a white background                 
  theme(legend.position = "none")+
  scale_color_manual("grey50")

