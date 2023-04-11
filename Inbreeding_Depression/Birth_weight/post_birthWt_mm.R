
library(tidyverse)

load(file="PhD_4th_yr/Inbreeding_depression_models/birth_weight/birth_wt_model_output_incmumage_50kit_50bin.RData")


summary_table=summary(birth_wt_model)$solutions
summary_table

FROH_sum_sol=summary(birth_wt_model)$solutions[10,1]



sols_full<-as.data.frame(birth_wt_model$Sol)%>%dplyr::select(matches("FROH"))%>% ## taking out sols with FROH included
  dplyr::mutate(across(2:34, ~.x + FROHsum)) ## adding FROHsum to chrFROH values

names <- sols_full %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_full,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(sols_full,2,quantile,probs = c(0.95)) #gets upper confidence interval for all solutions 
CI_lower<-apply(sols_full,2,quantile,probs = c(0.05)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH

FROH_sols$CHR<-as.factor(FROH_sols$CHR)



ggplot(data=FROH_sols, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange(colour="grey50") + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  geom_hline(yintercept=FROH_sum_sol, lty=1,colour="red") + ## red line is the average effect of all chromosomes 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="solution + CI", title = "Chromosome FROH on BirthWt (inc mum age)") +
  theme_bw()+  # use a white background                 
  theme(legend.position = "none")


###############################################################################################################
#########################################################################################
##########################################################################################################


# model Sex + AgeHrs + MotherStatus + mum_age+mum_age_sq+FROHsum+(mm chr FROH)+mumcode(random)+birthyr(random)

#read in table to see what averages are
birth_wt_df=read.table("PhD_4th_yr/Inbreeding_depression_models/birth_weight/birth_wt_df2023.txt", sep=",", header=T)
birth_wt_df_na_rm=birth_wt_df%>%
  dplyr::select(-mum_birthyear)%>%na.omit()

##check mean birth weight
mean(birth_wt_df_na_rm$CaptureWt) ##this is what we are predicting 
#7.097226
#mean(birth_wt_df_na_rm$BirthWt)
#6.483708

FROH_full<-read.table("PhD_4th_yr/2023_ROH_search/2021_calves_ROH_UpdatedSortedMb_032023.hom.indiv", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,KB) %>% dplyr::rename(Code=IID)%>%mutate(FROH=KB/2591865)%>%filter(nchar(Code)==5)

## finding means coefficients from df used 
Sex=1#1=female 2=male
AgeHrs=median(birth_wt_df_na_rm$AgeHrs)#median is 24 hrs 
mum_age=mean(birth_wt_df_na_rm$mum_age) #mean mum age in dataset
mum_age_sq=mean(birth_wt_df_na_rm$mum_age_sq)


ibc_qua=unname(quantile(FROH_full$FROH, probs = seq(0, 1, 1/20)))#getting the quantiles of FROH values for all ids
##also choose quantiles 
ibc_qua=c(0,0.05,0.1,0.15,0.2,0.25,0.3)
# ibc_qua=seq(0.05,0.25, 0.01)
# 




pred_ibcs=list()

for(v in 1:length(ibc_qua)){

    ibc=ibc_qua[v] #picking inbreeding coeff from quantiles vector
    FROHsum=(ibc*33)##choose inbreeding coefficient and * by the number of chr (which is what FROHsum is)
    FROHchrs=ibc
    
    bw_pred_list=list()
    
    for (i in 1:33){
      bw_pred=summary_table[1,1]+#intercept
        (Sex *(summary_table[2,1]))+ 
        (AgeHrs*(summary_table[3,1])) + #median age in hrs * coefficient predicted from model
        #(MotherStatus(summary_table[5,1])) + 
        (mum_age*(summary_table[8,1]))+
        (mum_age_sq*(summary_table[9,1]))+
        (FROHsum*(summary_table[10,1]))+
        (FROHchrs*(as.numeric(FROH_sols[i,1])))
      bw_pred_list[i]=bw_pred
    }
    
    bw_pred_chr=do.call(rbind.data.frame,bw_pred_list)%>%
      tibble::rownames_to_column("CHR")%>%
      add_column(paste0(round(ibc, digits = 4)))
     
    colnames(bw_pred_chr) <- c("CHR", "bw_prediction","ibc")
    
    pred_ibcs[[v]]=bw_pred_chr
    
}


#pred_ibcs[[1]]


pred_ibcs_all=do.call(rbind.data.frame,pred_ibcs)


ggplot(data=pred_ibcs_all,aes(x=ibc,y=bw_prediction, group=CHR, colour=CHR))+
  geom_line()+
  theme_bw()
