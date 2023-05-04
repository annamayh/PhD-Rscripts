
library(tidyverse)
library(data.table)

setwd("H:/")

load(file="PhD_4th_yr/Inbreeding_depression_models/survival/juvenile_survival_model_output_full.RData")

summary(juvenile_surv_model)

summary_table=summary(juvenile_surv_model)$solutions
summary_table

FROH_sum_sol=summary(juvenile_surv_model)$solutions[7,1]
FROH_sum_upr=summary(juvenile_surv_model)$solutions[7,2]
FROH_sum_lwr=summary(juvenile_surv_model)$solutions[7,3]


mean(rowSums(juvenile_surv_model$VCV))


sols_full<-as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(matches("FROH_"))%>% ## taking out sols with FROH included
 dplyr::mutate(across(1:33, ~.x + FROH_sum_sol)) ## adding FROHsum to chrFROH values

names <- sols_full %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_full,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(sols_full,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
CI_lower<-apply(sols_full,2,quantile,probs = c(0.025)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH

FROH_sols$CHR<-as.factor(FROH_sols$CHR)


FROH_sols$overlap_zero<-ifelse(FROH_sols$CHR %in% c("12"),"Yes","No")


surv_forest=ggplot(data=FROH_sols, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper, color = overlap_zero)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  geom_hline(yintercept=FROH_sum_sol, lty=1,colour="red") + ## red line is the average effect of all chromosomes 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="Posterior estimate + CI", title="Deviation of chromosomal inbreeding effects \nfrom combined effect of all chromosomes") +
  theme_classic()+  # use a white background                 
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey50","red"),guide=FALSE)+
  annotate("segment", x = 0, xend = 34, y = FROH_sum_sol, yend = FROH_sum_sol, colour = "red", alpha=0.6)+
  ##add in estimate of FROHsum
  annotate("pointrange", x = 34, y = FROH_sum_sol, ymin = FROH_sum_upr, ymax = FROH_sum_lwr,
           colour = "red", linewidth = 1, alpha=0.5, size=0.2)+
  annotate("segment", x = 35, xend = 35, y = 0, yend = FROH_sum_sol, colour = "red", alpha=0.5, linewidth=1)+
  #two lines between sig line 
  annotate("segment", x = 35, xend = 34.5, y = 0,  yend=0, colour = "red", alpha=0.5, linewidth=1)+
  annotate("segment", x = 35, xend = 34.5, y = FROH_sum_sol, yend = FROH_sum_sol, colour = "red", alpha=0.5, linewidth=1)+
  ##significance of FROHsum
  geom_text(aes(x=35.5, y=-0.14, label="***"), colour="red", alpha=0.5)+
  expand_limits(x = 36)

surv_forest

# 
# ggsave(surv_forest, 
#        file = "PhD_4th_yr/Inbreeding_depression_models/survival/chapter plots/surv_forest.png", 
#        width = 7, 
#        height = 8)
# 


sum(FROH_sols$solution)



## plot if you were totally outbred on all other chromosomes 

  
  ################################################################################################
  ################################################################################################

Sex=1.5#1=female 2=male good to just take average 
mum_age=mean(juvenile_surv_df_na_rm$mum_age)
mum_age_sq=mum_age^2

    
summary_table
    
inter=as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(matches("Interc"))


sex_ests=as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(matches("Sex"))
sex_ests=sex_ests*Sex

mumage_ests=as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(("mum_age"))
mumage_ests=mumage_ests*mum_age

mumagesq_ests=as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(("mum_age_sq"))
mumagesq_ests=mumagesq_ests*(mum_age^2)

all_effects=as.matrix(rowSums(inter%>%cbind(sex_ests)%>%cbind(mumage_ests)%>%cbind(mumagesq_ests)))

#now effects of inbreeding 
FROHsumest=as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(matches("FROHs"))

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
    labs(x="Inbreeding Coefficient", y="Predicted survival", colour="Chromosome", title="Predicted survival assuming chromosome independence")

suvr_pred_independent




#################################################################################################
## now assumng same inbreeding on all chrs 
##############################################################################################
chr_sols=as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(matches("FROH_"))
chr_sols=as.matrix(sols_full)


pred_matrix_all=list()


  for(v in 1:length(ibc_qua)){
    ibc=ibc_qua[v]
    
    if (ibc==0){
      #when ibc=0 does work as all ibc are 0   
      all_its=all_effects
      
    }else{
      
      sols_mult_ibc=chr_sols*ibc
      FROH_sums_mult=ibc*FROHsumest*33
      
      total_ibc_eff=sols_mult_ibc%>%cbind(FROH_sums_mult)
      total_ibc_eff=rowSums(total_ibc_eff)#add all together
      
      all_its=total_ibc_eff+all_effects
    }
    
    transformed=apply(all_its,2,transform_prediction)
    
    mean_surv=mean(transformed)
    quantU_surv=quantile(transformed, prob=c(0.975)) 
    quantL_surv=quantile(transformed, prob=c(0.025)) 
    pred_mat=matrix(c(mean_surv,quantU_surv,quantL_surv,(paste0(round(ibc, digits = 4))), paste0(chr)), nrow=1, ncol=5)
    
    pred_matrix_all[[v]]=pred_mat
  }
  
  pred_all=do.call(rbind.data.frame,pred_matrix_all)
  colnames(pred_all) <- c("mean_surv_predictin","CI_upr" ,"CI_lwr" ,"ibc", "CHR")



pred_survival_nonind_chr <- sapply(pred_all, as.numeric)%>%as.data.frame()

combined_chr=ggplot(data=pred_survival_nonind_chr,aes(x=ibc,y=mean_surv_predictin))+
  geom_line(size=1)+
  theme_bw()+
  labs(x="Inbreeding coefficient", y="Predicted survival", title="Predicted survival assuming equal inbreeding \non all chromosomes")+
  geom_ribbon(aes(ymin = CI_lwr, ymax = CI_upr), alpha = 0.3)

combined_chr


library(patchwork)

surv_3in1=surv_forest+combined_chr+suvr_pred_independent+plot_annotation(tag_levels = 'A')
surv_3in1


ggsave(surv_3in1,
       file = "PhD_4th_yr/Inbreeding_depression_models/survival/chapter plots/survival_all3.png",
       width = 17,
       height = 6)



