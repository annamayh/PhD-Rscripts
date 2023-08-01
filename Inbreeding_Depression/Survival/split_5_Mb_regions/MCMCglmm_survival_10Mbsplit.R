library(tidyverse)
library("MCMCglmm")
library(data.table)

setwd("H:/")

juvenile_surv_df_10mb=read.table("PhD_4th_yr/Inbreeding_depression_models/survival/10MBsplit_juvenile_survival_df.txt", sep=",", header=T)

#head(juvenile_surv_df)
juvenile_surv_df_na_rm=juvenile_surv_df_10mb%>%
  select(-mum_birthyear)%>% #removing birthwt from df cos not going to fit it this time and dont need mum birth yr
  #filter(!Sex=="3")%>%
  na.omit() #remove all NAs
## change some variables to factors 
juvenile_surv_df_na_rm$BirthYear=as.factor(juvenile_surv_df_na_rm$BirthYear)
juvenile_surv_df_na_rm$Code=as.factor(juvenile_surv_df_na_rm$Code)
juvenile_surv_df_na_rm$MumCode=as.factor(juvenile_surv_df_na_rm$MumCode)
juvenile_surv_df_na_rm$MotherStatus=as.factor(juvenile_surv_df_na_rm$MotherStatus)
juvenile_surv_df_na_rm$Sex=as.factor(juvenile_surv_df_na_rm$Sex)

k<-100
prior<-list(R=list(V=1,fix=1),
            G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=k), ## multimemberhsip part
                   G1=list(V=1,nu=1,alpha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,alpha.mu=0,alpha.V=k)))





juvenile_surv_model<-MCMCglmm(juvenile_survival~1 + Sex + MotherStatus + FROH_split_sum + mum_age+mum_age_sq, #need to fit sum chrFROH  as continuous covariate,
                              random= ~ idv(Froh_c_1_w_1+    Froh_c_1_w_11+Froh_c_1_w_21+Froh_c_1_w_31+Froh_c_1_w_41+Froh_c_1_w_51+
                      Froh_c_1_w_61+  Froh_c_1_w_71+Froh_c_1_w_81+Froh_c_2_w_1+Froh_c_2_w_11+Froh_c_2_w_21+ 
                      Froh_c_2_w_31+  Froh_c_2_w_41+Froh_c_3_w_1+Froh_c_3_w_11+Froh_c_3_w_21+Froh_c_3_w_31+  
                      Froh_c_3_w_41+  Froh_c_3_w_51+Froh_c_4_w_1+Froh_c_4_w_11+Froh_c_4_w_21+Froh_c_4_w_31+  
                      Froh_c_4_w_41+  Froh_c_4_w_51+Froh_c_4_w_61+Froh_c_5_w_1+Froh_c_5_w_11+Froh_c_5_w_21+  
                      Froh_c_5_w_31+  Froh_c_5_w_41+Froh_c_5_w_51+Froh_c_5_w_61+Froh_c_5_w_71+Froh_c_5_w_81+  
                      Froh_c_5_w_91+  Froh_c_5_w_101+Froh_c_5_w_111+Froh_c_5_w_121+Froh_c_5_w_131+Froh_c_6_w_1+   
                      Froh_c_6_w_11+  Froh_c_6_w_21+Froh_c_6_w_31+Froh_c_6_w_41+Froh_c_7_w_1+Froh_c_7_w_11+  
                      Froh_c_7_w_21+  Froh_c_7_w_31+Froh_c_8_w_1+Froh_c_8_w_11+Froh_c_8_w_21+Froh_c_8_w_31+  
                      Froh_c_8_w_41+  Froh_c_8_w_51+Froh_c_9_w_1+Froh_c_9_w_11+Froh_c_9_w_21+Froh_c_9_w_31+  
                      Froh_c_9_w_41+  Froh_c_9_w_51+Froh_c_9_w_61+Froh_c_9_w_71+Froh_c_9_w_81+Froh_c_9_w_91+  
                      Froh_c_9_w_101+ Froh_c_10_w_1+Froh_c_10_w_11+Froh_c_10_w_21+Froh_c_10_w_31+Froh_c_10_w_41+ 
                      Froh_c_10_w_51+ Froh_c_11_w_1+Froh_c_11_w_11+Froh_c_11_w_21+Froh_c_11_w_31+Froh_c_11_w_41+ 
                      Froh_c_11_w_51+ Froh_c_11_w_61+Froh_c_11_w_71+Froh_c_11_w_81+Froh_c_11_w_91+Froh_c_11_w_101+
                      Froh_c_12_w_1+  Froh_c_12_w_11+Froh_c_12_w_21+Froh_c_12_w_31+Froh_c_12_w_41+Froh_c_12_w_51+ 
                      Froh_c_12_w_61+ Froh_c_12_w_71+Froh_c_12_w_81+Froh_c_12_w_91+Froh_c_13_w_1+Froh_c_13_w_11+ 
                      Froh_c_13_w_21+ Froh_c_13_w_31+Froh_c_13_w_41+Froh_c_13_w_51+Froh_c_13_w_61+Froh_c_13_w_71+ 
                      Froh_c_14_w_1+  Froh_c_14_w_11+Froh_c_14_w_21+Froh_c_14_w_31+Froh_c_14_w_41+Froh_c_14_w_51+ 
                      Froh_c_14_w_61+ Froh_c_14_w_71+Froh_c_14_w_81+Froh_c_15_w_1+Froh_c_15_w_11+Froh_c_15_w_21+ 
                      Froh_c_15_w_31+ Froh_c_15_w_41+Froh_c_15_w_51+Froh_c_15_w_61+Froh_c_15_w_71+Froh_c_15_w_81+ 
                      Froh_c_15_w_91+ Froh_c_16_w_1+Froh_c_16_w_11+Froh_c_16_w_21+Froh_c_16_w_31+Froh_c_17_w_1+  
                      Froh_c_17_w_11+ Froh_c_17_w_21+Froh_c_17_w_31+Froh_c_17_w_41+Froh_c_17_w_51+Froh_c_17_w_61+ 
                      Froh_c_18_w_1+  Froh_c_18_w_11+Froh_c_18_w_21+Froh_c_18_w_31+Froh_c_18_w_41+Froh_c_18_w_51+ 
                      Froh_c_18_w_61+ Froh_c_18_w_71+Froh_c_18_w_81+Froh_c_18_w_91+Froh_c_18_w_101+Froh_c_18_w_111+
                      Froh_c_18_w_121+Froh_c_19_w_1+Froh_c_19_w_11+Froh_c_19_w_21+Froh_c_19_w_31+Froh_c_19_w_41+ 
                      Froh_c_19_w_51+ Froh_c_19_w_61+Froh_c_19_w_71+Froh_c_19_w_81+Froh_c_19_w_91+Froh_c_20_w_1+  
                      Froh_c_20_w_11+ Froh_c_20_w_21+Froh_c_20_w_31+Froh_c_20_w_41+Froh_c_20_w_51+Froh_c_20_w_61+ 
                      Froh_c_20_w_71+ Froh_c_20_w_81+Froh_c_20_w_91+Froh_c_20_w_101+Froh_c_20_w_111+Froh_c_20_w_121+
                      Froh_c_20_w_131+Froh_c_21_w_1+Froh_c_21_w_11+Froh_c_21_w_21+Froh_c_21_w_31+Froh_c_21_w_41+ 
                      Froh_c_21_w_51+Froh_c_21_w_61+Froh_c_21_w_71+Froh_c_21_w_81+Froh_c_22_w_1+Froh_c_22_w_11+ 
                      Froh_c_22_w_21+Froh_c_22_w_31+Froh_c_22_w_41+Froh_c_22_w_51+Froh_c_23_w_1+Froh_c_23_w_11+ 
                      Froh_c_23_w_21+Froh_c_23_w_31+Froh_c_23_w_41+Froh_c_23_w_51+Froh_c_23_w_61+Froh_c_23_w_71+ 
                      Froh_c_23_w_81+Froh_c_24_w_1+Froh_c_24_w_11+Froh_c_24_w_21+Froh_c_24_w_31+Froh_c_24_w_41+ 
                      Froh_c_24_w_51+Froh_c_24_w_61+Froh_c_24_w_71+Froh_c_25_w_1+Froh_c_25_w_11+Froh_c_25_w_21+ 
                      Froh_c_25_w_31+Froh_c_25_w_41+Froh_c_25_w_51+Froh_c_25_w_61+Froh_c_25_w_71+Froh_c_26_w_1+  
                      Froh_c_26_w_11+Froh_c_26_w_21+Froh_c_26_w_31+Froh_c_26_w_41+Froh_c_27_w_1+Froh_c_27_w_11+ 
                      Froh_c_27_w_21+Froh_c_27_w_31+Froh_c_27_w_41+Froh_c_27_w_51+Froh_c_27_w_61+Froh_c_28_w_1+  
                      Froh_c_28_w_11+Froh_c_28_w_21+Froh_c_28_w_31+Froh_c_28_w_41+Froh_c_28_w_51+Froh_c_28_w_61+ 
                      Froh_c_29_w_1+ Froh_c_29_w_11+Froh_c_29_w_21+Froh_c_29_w_31+Froh_c_29_w_41+Froh_c_29_w_51+ 
                      Froh_c_29_w_61+Froh_c_30_w_1+Froh_c_30_w_11+Froh_c_30_w_21+Froh_c_30_w_31+Froh_c_30_w_41+ 
                      Froh_c_30_w_51+Froh_c_30_w_61+Froh_c_30_w_71+Froh_c_31_w_1+Froh_c_31_w_11+Froh_c_31_w_21+ 
                      Froh_c_31_w_31+Froh_c_31_w_41+Froh_c_31_w_51+Froh_c_32_w_1+Froh_c_32_w_11+Froh_c_32_w_21+ 
                      Froh_c_32_w_31+Froh_c_32_w_41+Froh_c_33_w_1+Froh_c_33_w_11+Froh_c_33_w_21+Froh_c_33_w_31+ 
                      Froh_c_33_w_41+Froh_c_33_w_51+Froh_c_33_w_61+Froh_c_33_w_71)+BirthYear +MumCode,
                              family="threshold",
                              data=juvenile_surv_df_na_rm,
                              prior = prior,
                              pr=TRUE,#saves posterior dist for random effects i.e. what we want
                              nitt=200000,burnin=50000
)##





plot(juvenile_surv_model)







save(juvenile_surv_model, file="PhD_4th_yr/Inbreeding_depression_models/survival/10Mb_split_200k.RData")





summary_table=summary(juvenile_surv_model)$solutions
summary_table

FROH_sum_sol=summary(juvenile_surv_model)$solutions[8,1]

mean(rowSums(juvenile_surv_model$VCV))


sols_full<-as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(matches("FROH_"))%>% ## taking out sols with FROH included
dplyr::mutate(across(2:256, ~.x + FROH_split_sum)) ## adding FROHsum to chrFROH values

names <- sols_full %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_full,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(sols_full,2,quantile,probs = c(0.95)) #gets upper confidence interval for all solutions 
CI_lower<-apply(sols_full,2,quantile,probs = c(0.05)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols<-Random_table%>%filter(model_variable %like% "Froh_c")%>% add_column(CHR = 1:255) ##filtering all random variables for those including FROH

FROH_sols$CHR<-as.factor(FROH_sols$CHR)




ggplot(data=FROH_sols, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  geom_hline(yintercept=FROH_sum_sol, lty=1,colour="red") + ## red line is the average effect of all chromosomes 
 # coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="solution + CI", title = "Chromosome FROH on survival (Exl shot ids that didnt make it to adulthood)") +
  theme_bw()+  # use a white background                 
  theme(legend.position = "none")+
  scale_color_manual(values=c("grey50","red"),guide=FALSE)





