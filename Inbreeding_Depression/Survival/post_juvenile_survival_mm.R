
library(tidyverse)
library(data.table)

setwd("H:/")

load(file="PhD_4th_yr/Inbreeding_depression_models/survival/juvenile_survival_model_output_full_exBW.RData")

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


ggplot(data=FROH_sols, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper, color = overlap_zero)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  geom_hline(yintercept=FROH_sum_sol, lty=1,colour="red") + ## red line is the average effect of all chromosomes 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="solution + CI", title = "Chromosome FROH on survival (inc BW, including shot ids that survived to adult)") +
  theme_bw()+  # use a white background                 
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey50","red"),guide=FALSE)



sum(FROH_sols$solution)



###################################################################################################################
###################################################################################################################


Sex=1#1=female 2=male
mum_age=mean(juvenile_surv_df_na_rm$mum_age)
mum_age_sq=mum_age^2
#FROHsum=mean(juvenile_surv_df_na_rm$FROHsum)
#mum_age_sq=mean(juvenile_surv_df_na_rm$mum_age_sq)
#BirthWt=mean(juvenile_surv_df_na_rm$BirthWt)
# 
# FROH_full<-read.table("PhD_4th_yr/2023_ROH_search/2021_calves_ROH_UpdatedSortedMb_032023.hom.indiv", header=T, stringsAsFactors = F)%>%
#   dplyr::select(IID,KB) %>% dplyr::rename(Code=IID)%>%mutate(FROH=KB/2591865)%>%filter(nchar(Code)==5)
# ibc_qua=unname(quantile(FROH_full$FROH, probs = seq(0, 1, 1/20)))#getting the quantiles of FROH values for all ids

ibc_qua=c(0,0.05,0.1,0.15,0.2,0.25,0.3)

pred_ibcs=list()

for(v in 1:length(ibc_qua)){

        ibc=ibc_qua[v]
        #FROHsum=(ibc*33)##choose inbreeding coefficient and * by the number of chr (which is what FROHsum is)
        FROHchrs=ibc


        surv_pred_list=list()
        for (i in 1:33){
          
        
            surv_pred=summary_table[1,1]+#intercept
              (Sex *(summary_table[2,1]))+ 
              ((summary_table[5,1])) + 
              (mum_age*(summary_table[8,1]))+
              (mum_age_sq*(summary_table[9,1]))+
             # (FROHsum*(summary_table[7,1]))+
              #(BirthWt*(summary_table[10,1]))+
              (FROHchrs*(as.numeric(FROH_sols[i,1])))
            
            surv_pred_list[[i]]=pnorm(surv_pred,0,mean(rowSums(juvenile_surv_model$VCV)))
            
        }
        
        surv_pred_chr=do.call(rbind.data.frame,surv_pred_list)%>%
          tibble::rownames_to_column("CHR")%>%
          add_column(paste0(round(ibc, digits = 4)))
        
        colnames(surv_pred_chr) <- c("CHR", "surv_prediction","ibc")
        
        pred_ibcs[[v]]=surv_pred_chr
        
        }
            



pred_ibcs_all=do.call(rbind.data.frame,pred_ibcs)


ggplot(data=pred_ibcs_all,aes(x=ibc,y=surv_prediction, group=CHR, colour=CHR))+
  geom_line(linewidth=1)+
  theme_bw()#+
  ylim(0,1)



  
  
  
  
  ################################################################################################
  ################################################################################################
  
  
  sols_fixed_df=as.data.frame(juvenile_surv_model$Sol[,1:42])
  i=1
  
  row=sols_fixed_df[1,]
  
  
  surv_pred_list=list()

    chr_col=9+i
    surv_pred=row$`(Intercept)`+#intercept
      (Sex *(row$Sex2))+ 
      ((row$`MotherStatusTrue yeld`)) + 
      (mum_age*(row$mum_age))+
      (mum_age_sq*(row$mum_age_sq))+
      
      (FROHchrs*row[chr_col])
    
    surv_pred_list[[i]]=pnorm(surv_pred,0,mean(rowSums(juvenile_surv_model$VCV)))
    
  
