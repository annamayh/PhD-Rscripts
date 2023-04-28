
library(tidyverse)
library(data.table)

setwd("H:/")

load(file="PhD_4th_yr/Inbreeding_depression_models/survival/juvenile_survival_model_output_full.RData")

summary(juvenile_surv_model)

summary_table=summary(juvenile_surv_model)$solutions
summary_table

FROH_sum_sol=summary(juvenile_surv_model)$solutions[7,1]

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


(surv_forest=ggplot(data=FROH_sols, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper, color = overlap_zero)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  geom_hline(yintercept=FROH_sum_sol, lty=1,colour="red") + ## red line is the average effect of all chromosomes 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="solution + CI") +
  theme_bw()+  # use a white background                 
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey50","red"),guide=FALSE))



ggsave(surv_forest, 
       file = "PhD_4th_yr/Inbreeding_depression_models/survival/chapter plots/surv_forest.png", 
       width = 7, 
       height = 8)



sum(FROH_sols$solution)



###################################################################################################################
###################################################################################################################


Sex=1#1=female 2=male
mum_age=mean(juvenile_surv_df_na_rm$mum_age)
mum_age_sq=mum_age^2

## plot if you were totally outbred on all other chromosomes 

  
  ################################################################################################
  ################################################################################################

  ibc_qua=c(0,0.05,0.1,0.15,0.2,0.25,0.3)
  
  sols_fixed_df=as.data.frame(juvenile_surv_model$Sol[,1:42])%>%
    dplyr::mutate(across(10:42, ~.x + FROH_sum_sol))

  
  survival_prediction_table=list()
  
  for(v in 1:length(ibc_qua)){
    
   ibc=ibc_qua[v]
    
  surv_pred_all_rows=list()
  

    for (i in 1:nrow(sols_fixed_df)){

          surv_pred_one_row=list()
          
         for (ch in 1:33){
           
            chr_col=9+ch
            row=sols_fixed_df[i,]
            
             surv_pred=row$`(Intercept)`+#intercept
              (Sex *(row$Sex2))+ 
              ((row$`MotherStatusTrue yeld`)) + 
              (mum_age*(row$mum_age))+
              (mum_age_sq*(row$mum_age_sq))+
              
              (ibc*(as.numeric(row[chr_col])))
            
             surv_pred_one_row[[ch]]=pnorm(surv_pred,0,mean(rowSums(juvenile_surv_model$VCV)))
            #surv_pred_list
         }
        
          surv_pred_one_row_bind=do.call(cbind.data.frame,surv_pred_one_row)
          colnames(surv_pred_one_row_bind) <- c(paste0("chr_predict", 1:33))
          surv_pred_all_rows[[i]]=surv_pred_one_row_bind
    }  
  
  surv_pred_all_rows_bind=do.call(rbind.data.frame,surv_pred_all_rows)

  prediction_est<-apply(surv_pred_all_rows_bind,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
  CI_upper<-apply(surv_pred_all_rows_bind,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
  CI_lower<-apply(surv_pred_all_rows_bind,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
  
  survival_prediction_table[[v]]<-tibble(prediction_est)%>%add_column(CI_upper)%>%
   add_column(CI_lower)%>%
   add_column(paste0(ibc))%>%
   tibble::rownames_to_column("CHR")
   
  
  }
    
  
  
survival_prediction_all=do.call(rbind.data.frame,survival_prediction_table)

names(survival_prediction_all)[5]="ibc"

survival_prediction_all$CHR=as.numeric(survival_prediction_all$CHR)

survival_prediction_all=survival_prediction_all%>%arrange(CHR)

#survival_prediction_all$CHR=as.factor(survival_prediction_all$CHR)


(suvr_pred_independent=ggplot(data=survival_prediction_all,aes(x=ibc,y=prediction_est, group=as.factor(CHR), color=CHR))+
    geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper, fill = CHR, color = NULL), alpha = .035) +
    geom_line(linewidth=1)+
    theme_bw()+
    labs(x="Inbreeding Coefficient", y="Predicted survival on each chromosome")+
    theme(legend.position="none")
)

    
    
    ggsave(suvr_pred_independent, 
           file = "PhD_4th_yr/Inbreeding_depression_models/survival/chapter plots/surv_prediction_independent_chr.png", 
           width = 8, 
           height = 6)
    

    survival_prediction_all$CHR=as.factor(survival_prediction_all$CHR)
    
    
    
    (suvr_pred_independent2=ggplot(data=survival_prediction_all,aes(x=ibc,y=prediction_est, group=as.factor(CHR), color=CHR))+
        geom_line()+
        theme_bw()+
        labs(x="Inbreeding Coefficient", y="Predicted survival for an individual outbred on all other chromosomes", col="Chromosome")#+
       # theme(legend.position="none")
    )
    
    ggsave(suvr_pred_independent2, 
           file = "PhD_4th_yr/Inbreeding_depression_models/survival/chapter plots/surv_prediction_independent_chr2.png", 
           width = 6, 
           height = 6)
    
    