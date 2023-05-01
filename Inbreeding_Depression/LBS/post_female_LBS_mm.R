library(tidyverse)
library(data.table)
library(patchwork)


setwd("H:/")

load(file="PhD_4th_yr/Inbreeding_depression_models/LBS/Female_LBS_model_output_final.RData")

summary(female_LBS_model)
## in female LBS hurdle is significantly affected by inbreeding 
FROH_sum_sol_hu=summary(female_LBS_model)$solutions[4,1]
FROH_sum_hu_upr=summary(female_LBS_model)$solutions[4,2]
FROH_sum_hu_lwr=summary(female_LBS_model)$solutions[4,3]



### now hurdle 

sols_hu<-as.data.frame(female_LBS_model$Sol)%>%dplyr::select(matches("FROH_"))%>% ## taking out sols with FROH included
  dplyr::mutate(across(1:33, ~.x + FROH_sum_sol_hu)) ## adding FROHsum to chrFROH values


names <- sols_hu %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_hu,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(sols_hu,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
CI_lower<-apply(sols_hu,2,quantile,probs = c(0.025)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols_hu<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH

FROH_sols_hu$CHR<-as.factor(FROH_sols_hu$CHR)



female_LBS_f_hu=
  ggplot(data=FROH_sols_hu, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  #geom_vline(xintercept=33.5, lty=1) + ## red line is the average effect of all chromosomes 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="Posterior mean estimate of hurdle + CI ", title="Deviation of chromosomal inbreeding effcts from \ncombined effect of all chromosomes") +
  theme_classic()+  # use a white background                 
  theme(legend.position = "none")+
  scale_color_manual("grey50")+
  ## line of FROHsum estimate
  annotate("segment", x = 0, xend = 34, y = FROH_sum_sol_hu, yend = FROH_sum_sol_hu, colour = "red", alpha=0.6)+
  ##add in estimate of FROHsum
    annotate("pointrange", x = 34, y = FROH_sum_sol_hu, ymin = FROH_sum_hu_upr, ymax = FROH_sum_hu_lwr,
           colour = "red", linewidth = 1, alpha=0.5, size=0.2)+
  annotate("segment", x = 35, xend = 35, y = 0, yend = FROH_sum_sol_hu, colour = "red", alpha=0.5, linewidth=1)+
  #two lines between sig line 
  annotate("segment", x = 35, xend = 34.5, y = 0,  yend=0, colour = "red", alpha=0.5, linewidth=1)+
  annotate("segment", x = 35, xend = 34.5, y = FROH_sum_sol_hu, yend = FROH_sum_sol_hu, colour = "red", alpha=0.5, linewidth=1)+
  ##significance of FROHsum
  geom_text(aes(x=35.5, y=0.25, label="***"), colour="red", alpha=0.5)+
  expand_limits(x = 36)


female_LBS_f_hu





## for hurdle when chr considered seperatly 

ibc_qua=c(0,0.05,0.1,0.15,0.2,0.25,0.3)

pred_ibcs=list()



for(v in 1:length(ibc_qua)){
  
  ibc=ibc_qua[v]
    
  
  LBS_pred_list=list()
  for (i in 1:33){
    


        summary_table=summary(female_LBS_model)$solutions
        
        
        LBS_pred_hu=summary_table[2,1]+#intercept
          (ibc*(as.numeric(FROH_sols_hu[i,1])))
        
        
        LBS_pred_list[[i]]=(1-plogis(LBS_pred_hu,0,mean(rowSums(female_LBS_model$VCV))))
        ## 1 - Pr LBS = 0 
        
        }
        
  female_LBS_pred_chr=do.call(rbind.data.frame,LBS_pred_list)%>%
    tibble::rownames_to_column("CHR")%>%
    add_column(paste0(round(ibc, digits = 4)))
  
  colnames(female_LBS_pred_chr) <- c("CHR", "female_LBS_prediction_hurdle","ibc")
  
  pred_ibcs[[v]]=female_LBS_pred_chr
}




pred_ibcs_all=do.call(rbind.data.frame,pred_ibcs)

pred_ibcs_all$CHR=as.numeric(pred_ibcs_all$CHR)


independ_chr=pred_ibcs_all%>%
  arrange(CHR)%>%
  ggplot(aes(x=ibc,y=female_LBS_prediction_hurdle, group=as.factor(CHR), colour=as.factor(CHR)))+
  geom_line(linewidth=1)+
  theme_bw()+
  labs(x="Inbreeding coefficient", y="Predicted probability of LBS>0", title = "Predicted probability of LBS>0 assuming chromosome \nindependence", colour="Chromosome")

independ_chr


## for all chromosomes combined 
## need to account for the fact that FROHsum is also estimted in the models 


sols_hu_chrs=as.data.frame(female_LBS_model$Sol)%>%dplyr::select(("traithu_LBS"))

FROHsum_hu=as.data.frame(female_LBS_model$Sol)%>%dplyr::select(("traithu_LBS:FROHsum"))

FROHchr_hu=as.data.frame(female_LBS_model$Sol)%>%dplyr::select(matches("FROH_"))

transform_prediction=function(prediction){
LBS_pred_hu_all=(1-plogis(prediction,0,mean(rowSums(female_LBS_model$VCV))))
}


pred_nonind_hu=list()

for(v in 1:length(ibc_qua)){
  
 ibc=ibc_qua[v]
 
 #multiply all chr by ibc
 chr_multt_ibc=ibc*FROHchr_hu
 #multiply FROHsums by ibc and number of chrs (cos reasons)
 FROH_sums_mult=ibc*FROHsum_hu*33
 #add ibc effects together
 total_ibc_eff=chr_multt_ibc%>%cbind(FROH_sums_mult)
 total_ibc_eff=rowSums(total_ibc_eff)#add all together
 # add on to effect of hurdle 
 LBS_pred_hu=sols_hu_chrs+total_ibc_eff
 #transform every single iteration
 transformed<-apply(LBS_pred_hu,2,transform_prediction) #gets upper confidence interval for all solutions 
 
 mean=mean(transformed)
 quantU=quantile(transformed, prob=c(0.975)) 
 quantL=quantile(transformed, prob=c(0.025)) 
 
 pred_mat=matrix(c(mean,quantU,quantL,(paste0(round(ibc, digits = 4)))), nrow=1, ncol=4)
 
 
 pred_nonind_hu[[v]]=pred_mat  


}







pred_hu=do.call(rbind.data.frame,pred_nonind_hu)
colnames(pred_hu) <- c("female_LBS_prediction_hurdle","CI_upr" ,"CI_lwr" ,"ibc")
pred_hu <- sapply(pred_hu, as.numeric)%>%as.data.frame()


combined_chr=ggplot(data=pred_hu,aes(x=ibc,y=female_LBS_prediction_hurdle))+
  geom_line(size=1)+
  theme_bw()+
  labs(x="Inbreeding coefficient", y="Predicted probability of LBS>0", title="Predicted probability of LBS>0 assuming equal \ninbreeding on all chromosomes")+
  geom_ribbon(aes(ymin = CI_lwr, ymax = CI_upr), alpha = 0.3)

combined_chr


library(patchwork)

hurdle_3in1=female_LBS_f_hu+combined_chr+independ_chr+plot_annotation(tag_levels = 'A')
hurdle_3in1



# ggsave(hurdle_3in1,
#        file = "PhD_4th_yr/Inbreeding_depression_models/LBS/Chapter_plots/female_LBS_hu_3in1.png",
#        width = 17,
#        height = 6)
# 
# 


################################################################################################################
###################################################################################################################
########### FOR POISSON PROCESS #######################################################################
#####################################################################################################


FROH_sum_sol_poi=summary(female_LBS_model)$solutions[3,1]
FROH_sum_poi_upr=summary(female_LBS_model)$solutions[3,2]
FROH_sum_poi_lwr=summary(female_LBS_model)$solutions[3,3]



sols_pois<-as.data.frame(female_LBS_model$Sol)%>%dplyr::select(matches("FROH_"))%>% ## taking out sols with FROH included
  dplyr::mutate(across(1:33, ~.x + FROH_sum_sol_pois)) ## adding FROHsum to chrFROH values


names <- sols_pois %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_pois,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(sols_pois,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
CI_lower<-apply(sols_pois,2,quantile,probs = c(0.025)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols_pois<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH

FROH_sols_pois$CHR<-as.factor(FROH_sols_pois$CHR)


female_LBS_f_pois=
  ggplot(data=FROH_sols_pois, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  #geom_vline(xintercept=33.5, lty=1) + ## red line is the average effect of all chromosomes 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="Posterior mean estimate for truncated poisson process + CI ") +
  theme_classic()+  # use a white background                 
  theme(legend.position = "none")+
  scale_color_manual("grey50")+
  annotate("segment", x = 0, xend = 34, y = FROH_sum_sol_pois, yend = FROH_sum_sol_pois, colour = "red", alpha=0.6)+
  annotate("pointrange", x = 34, y = FROH_sum_sol_pois, ymin = FROH_sum_poi_upr, ymax = FROH_sum_poi_lwr,
           colour = "red", linewidth = 1, alpha=0.5, size=0.2)+
  expand_limits(x = 36)


female_LBS_f_pois



