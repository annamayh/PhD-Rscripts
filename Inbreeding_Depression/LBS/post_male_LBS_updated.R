

## jarrods function for converting male LBS to total LBS 



library(tidyverse)
library(data.table)
library(MCMCglmm)
library(patchwork)


setwd("H:/")

load(file="PhD_4th_yr/Inbreeding_depression_models/LBS/male_LBS_model_output_2.RData")

summary(male_LBS_model)
#plot(male_LBS_model)

FROH_sum_sol_zi=summary(male_LBS_model)$solutions[4,1]
FROH_sum_zi_upr=summary(male_LBS_model)$solutions[4,2]
FROH_sum_zi_lwr=summary(male_LBS_model)$solutions[4,3]


sols_zi<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(matches("2):FROH_"))%>% ## taking out sols with FROH included
  dplyr::mutate(across(1:33, ~.x + FROH_sum_sol_zi)) ## adding FROHsum to chrFROH values

names <- sols_zi %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_zi,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(sols_zi,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
CI_lower<-apply(sols_zi,2,quantile,probs = c(0.025)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols_zi<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH


male_LBS_zi=
  ggplot(data=FROH_sols_zi, aes(x=as.factor(CHR), y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="Posterior mean estimate + CI ") +
  theme_classic()+  # use a white background                 
  theme(legend.position = "none")+
  annotate("segment", x = 0, xend = 34, y = FROH_sum_sol_zi, yend = FROH_sum_sol_zi, colour = "mediumseagreen", alpha=0.6, linewidth=1)+
  ##add in estimate of FROHsum
  annotate("pointrange", x = 34, y = FROH_sum_sol_zi, ymin = FROH_sum_zi_upr, ymax = FROH_sum_zi_lwr,
           colour = "mediumseagreen", linewidth = 1, alpha=0.5, size=0.2)+
  expand_limits(x = 36)

male_LBS_zi


## same plot for poisson process

FROH_sum_sol_pois=summary(male_LBS_model)$solutions[3,1]
FROH_sum_pois_upr=summary(male_LBS_model)$solutions[3,2]
FROH_sum_pois_lwr=summary(male_LBS_model)$solutions[3,3]

sols_pois<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(matches("1):FROH_"))%>% ## taking out sols with FROH included
  dplyr::mutate(across(1:33, ~.x + FROH_sum_sol_pois)) ## adding FROHsum to chrFROH values

names <- sols_pois %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_pois,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(sols_pois,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
CI_lower<-apply(sols_pois,2,quantile,probs = c(0.025)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols_pois<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33)%>% ##filtering all random variables for those including FROH
  add_column(model = "poisson") 

male_LBS_pois=
  ggplot(data=FROH_sols_pois, aes(x=as.factor(CHR), y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="Posterior mean estimate + CI ") +
  theme_classic()+  # use a white background                 
  theme(legend.position = "none")+
  annotate("segment", x = 0, xend = 34, y = FROH_sum_sol_pois, yend = FROH_sum_sol_pois, colour = "chocolate3", alpha=0.6, linewidth=1)+
  ##add in estimate of FROHsum
  annotate("pointrange", x = 34, y = FROH_sum_sol_pois, ymin = FROH_sum_pois_upr, ymax = FROH_sum_pois_lwr,
           colour = "chocolate3", linewidth = 1, alpha=0.5, size=0.2)+
  annotate("segment", x = 35, xend = 35, y = FROH_sum_sol_pois, yend = 0, colour = "chocolate3", alpha=0.5, linewidth=1)+
  annotate("segment", x = 35, xend = 34.5, y = 0,  yend=0, colour = "chocolate3", alpha=0.5, linewidth=1)+
  annotate("segment", x = 35, xend = 34.5, y = FROH_sum_sol_pois, yend = FROH_sum_sol_pois, colour = "chocolate3", alpha=0.5, linewidth=1)+

  geom_text(aes(x=35.5, y=-0.38, label="***"), colour="chocolate3", alpha=0.5)+
  expand_limits(x = 36)

male_LBS_pois


male_LBS_zi+male_LBS_pois




###################################################################################
### transforming the data ########################################################
##################################################################################


#### function written by jarrod to tranform the data #### 
# Mu takes two values: the untrsanformed equation output for the poisson part followed by the untrsanformed 
# equation output for the Bernoulli part. V is a 2x2 covariance matrix describing the variance of the 
# random effects you want to average over. Because you fitted idh structures V is a diagonal matrix. 
# You can sum the variances over random terms (mum, year AND residuals) to get V.

normal.hupoisson<-function(mu, V){
  
  int.foo<-function(x, mu, V){(1-plogis(x[2]))*exp(x[1])/(1-exp(-exp(x[1])))*mvtnorm::dmvnorm(x, mu, V)}

  cubature::adaptIntegrate(int.foo, qnorm(0.0001, mu,sqrt(diag(V))), qnorm(0.9999, mu,sqrt(diag(V))), mu=mu, V=V)[[1]]
  
}


##V<-rIW(diag(2), 10) example


summary_table=summary(male_LBS_model)$solutions
summary_table

# G-structure:  ~idh(trait):MumCode
# 
# post.mean  l-95% CI u-95% CI eff.samp
# traitLBS.MumCode       0.2932 0.0002927   0.7501    272.2
# traitzi_LBS.MumCode    0.1578 0.0002135   0.7392    532.4
# 
# ~idh(trait):BirthYear
# 
# post.mean  l-95% CI u-95% CI eff.samp
# traitLBS.BirthYear      0.04717 0.0002359    0.202    59.95
# traitzi_LBS.BirthYear   3.19705 1.2145816    5.672  1056.13

# R-structure:  ~idh(trait):units
# 
# post.mean  l-95% CI u-95% CI eff.samp
# traitLBS.units       0.5488 0.0004595    1.007    227.6
# traitzi_LBS.units    1.0000 1.0000000    1.000      0.0
# 

#mu<-rnorm(2)


rand_var_zi=0.1578+3.19705+1.0000 ##variance estimatd for all radoms effects in zi

rand_var_pois=0.2932+0.04717+0.5488 ##variance estimatd for all radoms effects in pois

v=diag(c(rand_var_pois,rand_var_zi),2) ##assuming no covariance 


zi_mean_chr<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(matches("2):FROH_"))%>%## taking out sols with FROH included
  apply(2, mean)

pois_mean_chr<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(matches("1):FROH_"))%>%## taking out sols with FROH included
  apply(2, mean)

#as.numeric(zi_mean_chr[1])

ibc=0.1

chr=1

ibc_qua=c(0,0.05,0.1,0.15,0.2,0.25,0.3)

LBS_male_full=list()

for (c in 1:length(ibc_qua)){
  
  ibc=ibc_qua[c]
  
  LBS_male=list()
  
  for (chr in 1:33){
  
    zi_inter_mean=summary_table[2,1]
    zi_FROHsum_mean=summary_table[4,1]
    
    untrans_LBS_zi=(zi_inter_mean+(ibc*zi_FROHsum_mean*33)+(ibc*as.numeric(zi_mean_chr[chr])))#have to * 33 as FROHsum is sum of all 33 chrs 
    
    pois_inter_mean=summary_table[1,1]
    pois_FROHsum_mean=summary_table[3,1]
    
    untrans_LBS_pois=(pois_inter_mean+(ibc*pois_FROHsum_mean)+(ibc*as.numeric(pois_mean_chr[chr])))
    mu=(c(untrans_LBS_pois,untrans_LBS_zi))

    LBS_trans=normal.hupoisson(mu,v)
    LBS_trans_mat=matrix(c(LBS_trans, paste0(chr), paste0(ibc)), nrow = 1)
    
    LBS_male[[chr]]=LBS_trans_mat
      }

  LBS_male_all_list=do.call(rbind.data.frame,LBS_male)
  LBS_male_full[[c]]=LBS_male_all_list
  
}

LBS_male_all=do.call(rbind.data.frame,LBS_male_full)

names(LBS_male_all)=(c("LBS", "CHR", "ibc"))

LBS_male_all <- sapply(LBS_male_all, as.numeric)%>%as.data.frame()


LBS_pred_independent=ggplot(data=LBS_male_all,aes(x=ibc,y=LBS, group=as.factor(CHR), color=as.factor(CHR)))+
  geom_line(linewidth=1)+
  theme_bw()+
  labs(x="Inbreeding Coefficient", y="Predicted LBS", colour="Chromosome", title="Predicted male LBS assuming chromosome independence")

LBS_pred_independent

##need to add in real data too
