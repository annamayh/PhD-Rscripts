library(tidyverse)
library("MCMCglmm")
library(data.table)

setwd("H:/")

juvenile_surv_df_5mb=read.table("PhD_4th_yr/Inbreeding_depression_models/survival/split_regions/5split_juvenile_survival_df.txt", sep=",", header=T)

#head(juvenile_surv_df)
juvenile_surv_df_na_rm=juvenile_surv_df_5mb%>%
  select(-BirthWt)%>% #removing birthwt from df cos not going to fit it this time and dont need mum birth yr
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


FROH_mat=as.matrix(juvenile_surv_df_na_rm[9:506])

juvenile_surv_df_na_rm$FROH_mat=FROH_mat

juvenile_surv_model<-MCMCglmm(juvenile_survival~1 + Sex + MotherStatus + FROH_split_sum + mum_age+mum_age_sq, #need to fit sum chrFROH  as continuous covariate,
                              random= ~ idv(FROH_mat)+BirthYear +MumCode,
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





