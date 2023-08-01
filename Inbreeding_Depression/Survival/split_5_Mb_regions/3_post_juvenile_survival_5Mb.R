
library(tidyverse)
library(data.table)

setwd("H:/")

load(file="PhD_4th_yr/Inbreeding_depression_models/survival/split_regions/juvenile_surv_5Mb_output.RData")


summary_table

FROH_sum_sol=summary_table[7,1]
FROH_sum_upr=summary_table[7,2]
FROH_sum_lwr=summary_table[7,3]





sols_full<-as.data.frame(sols)%>%dplyr::select(matches("FROH_mat"))%>% ## taking out sols with FROH included
 dplyr::mutate(across(1:487, ~.x + FROH_sum_sol)) ## adding FROHsum to chrFROH values

names <- sols_full %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_full,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(sols_full,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
CI_lower<-apply(sols_full,2,quantile,probs = c(0.025)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols<-Random_table%>%add_column(window = 1:487) ##filtering all random variables for those including FROH

FROH_sols$window<-as.factor(FROH_sols$window)


#FROH_sols$overlap_zero<-ifelse(FROH_sols$CHR %in% c("12"),"Yes","No")


surv_forest=ggplot(data=FROH_sols, aes(x=window, y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  geom_hline(yintercept=FROH_sum_sol, lty=1, colour="red") +  # black line is 0
  #coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="Posterior estimate + CI", 
       title="Deviation of chromosomal inbreeding effects \nfrom combined effect of all chromosomes", tag = "B") +
  theme_classic()+  # use a white background                 
  theme(legend.position = "none")

surv_forest

# 
# ggsave(surv_forest, 
#        file = "PhD_4th_yr/Inbreeding_depression_models/survival/chapter plots/surv_forest.png", 
#        width = 7, 
#        height = 8)
# 


sum(FROH_sols$solution)

Sex=0.5#1=female 2=male good to just take average 
mum_age=mean(juvenile_surv_df_na_rm$mum_age)
mum_age_sq=mum_age^2
day_seq=33
#Bw=mean(juvenile_surv_df_na_rm$BirthWt)
#ibc=0.125


### working out for average calf
# fixed_means=(as.matrix(summary_table[c(1:2,7:9),1]))
# coeffs=matrix(c(1,Sex,(ibc*33),mum_age,mum_age_sq),ncol=1)
# pred=colSums(fixed_means*coeffs)
# pnorm(pred,0,mean(rowSums(juvenile_surv_model$VCV)))

#females
# 0.6946248 ibc=0
# 0.2864378 ibc=0.25
# 0.4890666 ibc=0.125

##MALES##
# 0.65143 ibc=0
# 0.2471012 ibc=0.25
# 0.441476 ibc = 0.125

## plot if you were totally outbred on all other chromosomes 
  ################################################################################################
  ################################################################################################


    
summary_table
    
inter=as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(matches("Interc"))

sex_ests=as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(matches("Sex"))
sex_ests=sex_ests*Sex

mumage_ests=as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(("mum_age"))
mumage_ests=mumage_ests*mum_age

mumagesq_ests=as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(("mum_age_sq"))
mumagesq_ests=mumagesq_ests*(mum_age^2)

day_seq_ests=as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(("Day_seq"))
day_seq_ests=day_seq_ests*day_seq


# 
# BWests=as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(("BirthWt"))
# BWests=BWests*Bw
# 

all_effects=as.matrix(rowSums(inter%>%cbind(sex_ests)%>%cbind(mumage_ests)%>%cbind(mumagesq_ests)%>%cbind(day_seq_ests)))

#chr_sols=as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(matches("FROH_"))
chr_sols=as.matrix(sols_full)

transform_prediction=function(prediction){
  survival_prediction_table=pnorm(prediction,0,mean(rowSums(juvenile_surv_model$VCV)))##mean=0 sd=random effects
}



ibc_qua=c(0,0.05,0.1,0.15,0.2,0.25,0.3)

#v=3
#chr=1  

pred_chrs_inde=list()

for (chr in 1:33){
  
  chr_col=chr_sols[,chr]#selecting first chr 
  pred_matrix_indep=list()
    
    for(v in 1:length(ibc_qua)){
      ibc=ibc_qua[v]
      
      if (ibc==0){
        #when ibc=0 does work as all ibc are 0   
        all_its=all_effects
        
      }else{
    
        sols_mult_ibc=chr_col*ibc
        all_its=sols_mult_ibc+all_effects
      }
    
      transformed=apply(all_its,2,transform_prediction)
      
      mean_surv=mean(transformed)
      quantU_surv=quantile(transformed, prob=c(0.975)) 
      quantL_surv=quantile(transformed, prob=c(0.025)) 
      pred_mat=matrix(c(mean_surv,quantU_surv,quantL_surv,(paste0(round(ibc, digits = 4))), paste0(chr)), nrow=1, ncol=5)
      
      pred_matrix_indep[[v]]=pred_mat
      }
    
  pred_chrs=do.call(rbind.data.frame,pred_matrix_indep)
  colnames(pred_chrs) <- c("mean_surv_predictin","CI_upr" ,"CI_lwr" ,"ibc", "CHR")
  pred_chrs_inde[[chr]]=pred_chrs
    
}
  

pred_survival_ind_chr=do.call(rbind.data.frame,pred_chrs_inde)
pred_survival_ind_chr <- sapply(pred_survival_ind_chr, as.numeric)%>%as.data.frame()


suvr_pred_independent=ggplot(data=pred_survival_ind_chr,aes(x=ibc,y=mean_surv_predictin, group=as.factor(CHR), color=as.factor(CHR)))+
    geom_line(linewidth=1)+
    theme_bw()+
    labs(x="Inbreeding Coefficient", y="Juvenile survival", 
         colour="Chromosome", title="Predicted survival assuming chromosome independence", tag = "C")+
  theme(text = element_text(size = 18), plot.title = element_text(size=13), 
        legend.title = element_text(size=10), 
        legend.text = element_text(size=9))

suvr_pred_independent




#################################################################################################
## now assumng same inbreeding on all chrs 
##############################################################################################

#now effects of inbreeding 

df_all_eff=all_effects%>%as.data.frame() #df of all other effect of survival created above 

random_effs=as.data.frame(rowSums(juvenile_surv_model$VCV)) ##variation of random effects




survival_transformed=list()


  for(v in 1:length(ibc_qua)){
    ibc=ibc_qua[v]
    
    pred_per_it=list()
    
    for (row in 1:nrow(df_all_eff)){
    
    iter=sum(chr_sols[row,]*ibc)#sum of all chromosomal effects x inbreeding coefficient 
      
    pred=df_all_eff%>%
      mutate(eff_plus_ibc=iter+V1, .keep = "none")
    
      prediction=as.numeric(pred[row,])
      random_effects=as.numeric(random_effs[row,])
      
      survival_prediction=pnorm(prediction,0,random_effects)#
      pred_per_it[[row]]=survival_prediction
      
    }
    
    surv_pred_unlist=do.call(rbind.data.frame,pred_per_it)
    
    mean<-apply(surv_pred_unlist,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
    CI_upper<-apply(surv_pred_unlist,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
    CI_lower<-apply(surv_pred_unlist,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
    
    surv_trans_mat=matrix(c(mean,CI_upper,CI_lower, paste0(ibc)),nrow = 1, ncol = 4)
    survival_transformed[[v]]=surv_trans_mat
  }
  

  pred_all=do.call(rbind.data.frame,survival_transformed)
  colnames(pred_all) <- c("mean_surv_prediction","CI_upr" ,"CI_lwr" ,"ibc")



pred_survival_nonind_chr <- sapply(pred_all, as.numeric)%>%as.data.frame()

juvenile_surv_df_na_rm_plot=juvenile_surv_df_na_rm%>%mutate(juvenile_survival=as.numeric(juvenile_survival))

combined_chr=ggplot(data=pred_survival_nonind_chr,aes(x=ibc,y=mean_surv_prediction))+
  geom_line(linewidth=1)+
  theme_bw()+
  labs(x="Inbreeding coefficient", y="Juvenile survival", 
       title="Predicted survival assuming equal inbreeding \non all chromosomes",tag = "A")+
  geom_ribbon(aes(ymin = CI_lwr, ymax = CI_upr), alpha = 0.3)+
  theme(text = element_text(size = 18), plot.title = element_text(size=13))+
  geom_point(data=juvenile_surv_df_na_rm_plot, aes(x=FROH_sum_div, y=juvenile_survival), inherit.aes = F, alpha=0.2)


combined_chr


# 0.6741455 - 0.2698610 
# 0.4042845/0.6741455

library(patchwork)

surv_3in1=(combined_chr|surv_forest|suvr_pred_independent)#+plot_annotation(tag_levels = 'A')
surv_3in1


ggsave(surv_3in1,
       file = "PhD_4th_yr/Inbreeding_depression_models/survival/chapter plots/survival_all3.png",
       width = 16,
       height = 6)



