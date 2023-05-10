

## jarrods function for converting male LBS to total LBS 



library(tidyverse)
library(data.table)
library(MCMCglmm)
library(patchwork)


setwd("H:/")

load(file="PhD_4th_yr/Inbreeding_depression_models/LBS/male_LBS_model_output_prior.6.RData")

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
  labs(x="Chromosome", y="Posterior mean estimate + CI ",title = "Zero inflation", tag = "B") +
  theme_classic()+  # use a white background                 
  theme(legend.position = "none")+
  annotate("segment", x = 0, xend = 34, y = FROH_sum_sol_zi, yend = FROH_sum_sol_zi, colour = "mediumseagreen", alpha=0.6, linewidth=1)+
  ##add in estimate of FROHsum
  annotate("pointrange", x = 34, y = FROH_sum_sol_zi, ymin = FROH_sum_zi_upr, ymax = FROH_sum_zi_lwr,
           colour = "mediumseagreen", linewidth = 1, alpha=0.5, size=0.2)+
  expand_limits(x = 35)

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

FROH_sols_pois<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33)#%>% ##filtering all random variables for those including FROH
  #add_column(model = "poisson") 

male_LBS_pois=
  ggplot(data=FROH_sols_pois, aes(x=as.factor(CHR), y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="Posterior mean estimate + CI ", title = "Poisson process") +
  theme_classic()+  # use a white background                 
  theme(legend.position = "none")+
  annotate("segment", x = 0, xend = 34, y = FROH_sum_sol_pois, yend = FROH_sum_sol_pois, colour = "chocolate3", alpha=0.6, linewidth=1)+
  ##add in estimate of FROHsum
  annotate("pointrange", x = 34, y = FROH_sum_sol_pois, ymin = FROH_sum_pois_upr, ymax = FROH_sum_pois_lwr,
           colour = "chocolate3", linewidth = 1, alpha=0.5, size=0.2)+

  geom_text(aes(x=34.5, y=-0.68, label="***"), colour="chocolate3", alpha=0.5)+
  expand_limits(x = 35)

male_LBS_pois


forest=male_LBS_zi+male_LBS_pois

forest


###################################################################################
### transforming the data ########################################################
##################################################################################


#### function written by jarrod to tranform the data #### 
# Mu takes two values: the untrsanformed equation output for the poisson part followed by the untrsanformed 
# equation output for the Bernoulli part. V is a 2x2 covariance matrix describing the variance of the 
# random effects you want to average over. Because you fitted idh structures V is a diagonal matrix. 
# You can sum the variances over random terms (mum, year AND residuals) to get V.

normal.zipooisson<-function(mu, V){
  
  int.foo<-function(x, mu, V){(1-plogis(x[2]))*exp(x[1])*mvtnorm::dmvnorm(x, mu, V)}
  
  
  cubature::adaptIntegrate(int.foo, qnorm(0.0001, mu,sqrt(diag(V))), qnorm(0.9999, mu,sqrt(diag(V))), mu=mu, V=V)[[1]]
  
}

##V<-rIW(diag(2), 10) example


summary_table=summary(male_LBS_model)$solutions
summary_table

summary(male_LBS_model)

# G-structure:  ~idh(trait):MumCode
# 
# post.mean  l-95% CI u-95% CI eff.samp
# traitLBS.MumCode       0.2887 8.535e-06   0.6842    580.6
# traitzi_LBS.MumCode    0.4543 9.402e-07   1.5940    446.0
# 
# ~idh(trait):BirthYear
# 
# post.mean  l-95% CI u-95% CI eff.samp
# traitLBS.BirthYear      0.05659 1.745e-07   0.2178    559.9
# traitzi_LBS.BirthYear   3.61229 1.366e+00   6.3368    639.9

# R-structure:  ~idh(trait):units
# 
# post.mean l-95% CI u-95% CI eff.samp
# traitLBS.units       0.6026   0.2432    1.073    689.9
# traitzi_LBS.units    1.0000   1.0000    1.000      0.0
# 

#mu<-rnorm(2)


rand_var_zi=0.4543+3.61229+1.0000 ##variance estimatd for radoms effects in zi: mum code, birth year and residual

rand_var_pois=0.2887+0.05659+0.6026 ##variance estimatd for all radoms effects in pois

v=diag(c(rand_var_pois,rand_var_zi),2) ##assuming no covariance 


zi_mean_chr<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(matches("2):FROH_"))%>%## taking out sols with FROH included
  apply(2, mean)

pois_mean_chr<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(matches("1):FROH_"))%>%## taking out sols with FROH included
  apply(2, mean)

#as.numeric(zi_mean_chr[1])

# ibc=0.1
# 
# chr=1

zi_inter_mean=summary_table[2,1]
zi_FROHsum_mean=summary_table[4,1]


pois_inter_mean=summary_table[1,1]
pois_FROHsum_mean=summary_table[3,1]

ibc_qua=c(0,0.05,0.1,0.15,0.2,0.25,0.3)

LBS_male_full=list()

for (c in 1:length(ibc_qua)){
  
  ibc=ibc_qua[c]
  
  LBS_male=list()
  
  for (chr in 1:33){

    untrans_LBS_zi=(zi_inter_mean+(ibc*(as.numeric(FROH_sols_zi[chr,1])))) 
    untrans_LBS_pois=(pois_inter_mean+(ibc*(as.numeric(FROH_sols_pois[chr,1]))))
    mu=(c(untrans_LBS_pois,untrans_LBS_zi))#  zi and pois effects

    LBS_trans=normal.zipooisson(mu,v) # function written by jarrod taking in to account variance 
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
  labs(x="Inbreeding Coefficient", y="Predicted LBS", colour="Chromosome", title="Predicted male LBS assuming chromosome independence",tag = "C")

LBS_pred_independent

##need to add in real data too



###############################################################
## now just assuming same inbreeding on all chrs using FROsum ####
#################################################################
# mean
zi_inter_mean=summary_table[2,1]
zi_FROHsum_mean=summary_table[4,1]
pois_inter_mean=summary_table[1,1]
pois_FROHsum_mean=summary_table[3,1]

ibc=0.125
mean_untrans_LBS_zi_sum=(zi_inter_mean+(ibc*zi_FROHsum_mean*33))#have to * 33 as FROHsum is sum of all 33 chrs
mean_untrans_LBS_pois_sum=(pois_inter_mean+(ibc*pois_FROHsum_mean*33))
mean_mu=(c(mean_untrans_LBS_pois_sum,mean_untrans_LBS_zi_sum))

LBS_trans=normal.zipooisson(mean_mu,v)
LBS_trans

#F=0 -F=0.125
10.29534-0.2761225 #difference between estimates
10.01922/10.29534 ##total decrese in LBS
# ##==97%

###########################################################################################################
# for all iterations to get 95% CIs 
######################################################################################

sol_pois<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(matches("tLBS:FROHsum"))## taking out sols with FROH included
sol_zi<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(matches("zi_LBS:FROHsum"))## taking out sols with FROH included
sol_inter_pois<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(("traitLBS"))## taking out sols of intercept
sol_inter_zi<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(("traitzi_LBS"))## taking out sols of intercept



LBS_transformed=list()

for (c in 1:length(ibc_qua)){
    
    ibc=ibc_qua[c]
    LBS_sum=list()
    
        for (row in 1:nrow(sol_pois)){
        
            untrans_LBS_zi_sum=(sol_inter_zi[row,]+(ibc*sol_zi[row,]*33))#have to * 33 as FROHsum is sum of all 33 chrs 
            untrans_LBS_pois_sum=(sol_inter_pois[row,]+(ibc*sol_pois[row,]*33))
            mu=(c(untrans_LBS_pois_sum,untrans_LBS_zi_sum))
            
            LBS_trans=normal.zipooisson(mu,v)
            LBS_sum[[row]]=LBS_trans
              }  
    
        LBS_sum_unlist=do.call(rbind.data.frame,LBS_sum)
        
        mean<-apply(LBS_sum_unlist,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
        CI_upper<-apply(LBS_sum_unlist,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
        CI_lower<-apply(LBS_sum_unlist,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
        
        LBS_trans_mat=matrix(c(mean,CI_upper,CI_lower, paste0(ibc)),nrow = 1, ncol = 4)
        LBS_transformed[[c]]=LBS_trans_mat
    
}

LBSsum_df=do.call(rbind.data.frame,LBS_transformed)
colnames(LBSsum_df) <- c("LBS_mean","CI_upr","CI_lwr","ibc")
LBSsum_df <- sapply(LBSsum_df, as.numeric)%>%as.data.frame()
    
    






LBS_pred_sum=ggplot(data=LBSsum_df,aes(x=ibc,y=LBS_mean))+
  geom_line(linewidth=1)+
  theme_bw()+
  labs(x="Inbreeding Coefficient", y="Predicted LBS", title="Predicted male LBS assuming equal inbreeding", tag = "A")+
  geom_ribbon(aes(ymin = CI_lwr, ymax = CI_upr), alpha = 0.3)


LBS_pred_sum

# male_LBS_df%>%
#   filter(FROH_sum_div==0)%>%
  mean(male_LBS_df$LBS)


library(patchwork)

LBS_all3=(LBS_pred_sum|forest|LBS_pred_independent)
LBS_all3

ggsave(LBS_all3,
       file = "PhD_4th_yr/Inbreeding_depression_models/LBS/Chapter_plots/male_LBS_all3.png",
       width = 16,
       height = 6)





#     
#     LBSsum=list()
# 
# for (c in 1:length(ibc_qua)){
#   
#   #ibc=0
#   ibc=ibc_qua[c]
#   untrans_LBS_zi_sum=(zi_inter_mean+(ibc*zi_FROHsum_mean*33))#have to * 33 as FROHsum is sum of all 33 chrs 
#   untrans_LBS_pois_sum=(pois_inter_mean+(ibc*pois_FROHsum_mean*33))
# 
#   mean_mu=(c(untrans_LBS_pois_sum,untrans_LBS_zi_sum))
# 
#   LBS_sum_mean=normal.hupoisson(mean_mu,v)
#   LBS_mat=matrix(c(LBS_sum_mean, paste0(ibc)),nrow=1,ncol=2)
#   
#   LBSsum[[c]]=LBS_mat
#   
# }
#   
# LBSsum_df=do.call(rbind.data.frame,LBSsum)
# colnames(LBSsum_df) <- c("LBS","ibc")
# LBSsum_df <- sapply(LBSsum_df, as.numeric)%>%as.data.frame()