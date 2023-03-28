library(tidyverse)
library("MCMCglmm")
library(data.table)

setwd("H:/")

juvenile_surv_df=read.table("PhD_4th_yr/Inbreeding_depression_models/survival/juvenile_survival_df.txt", sep=",", header=T)

#head(juvenile_surv_df)
juvenile_surv_df_na_rm=juvenile_surv_df%>%
  select(-mum_birthyear)%>% #removing birthwt from df cos not going to fit it this time and dont need mum birth yr
  na.omit() #remove all NAs
## change some variables to factors 
juvenile_surv_df_na_rm$BirthYear=as.factor(juvenile_surv_df_na_rm$BirthYear)
juvenile_surv_df_na_rm$Code=as.factor(juvenile_surv_df_na_rm$Code)
juvenile_surv_df_na_rm$MumCode=as.factor(juvenile_surv_df_na_rm$MumCode)
juvenile_surv_df_na_rm$MotherStatus=as.factor(juvenile_surv_df_na_rm$MotherStatus)
juvenile_surv_df_na_rm$Sex=as.factor(juvenile_surv_df_na_rm$Sex)



k<-100
prior<-list(R=list(V=1,fix=1),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k), ## multimemberhsip part
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))



juvenile_surv_model<-MCMCglmm(juvenile_survival~1 + Sex + MotherStatus + FROHsum + mum_age+mum_age_sq +BirthWt, #need to fit sum chrFROH  as continuous covariate,
                 random= ~ 
                   idv(FROH_chr1+FROH_chr2+FROH_chr3+FROH_chr4+FROH_chr5+FROH_chr6+FROH_chr7+FROH_chr8+
                         FROH_chr9+FROH_chr10+FROH_chr11+FROH_chr12+FROH_chr13+FROH_chr14+FROH_chr15+FROH_chr16+
                         FROH_chr17+FROH_chr18+FROH_chr19+FROH_chr20+FROH_chr21+FROH_chr22+FROH_chr23+FROH_chr24+
                         FROH_chr25+FROH_chr26+FROH_chr27+FROH_chr28+FROH_chr29+FROH_chr30+FROH_chr31+FROH_chr32+
                         FROH_chr33)+BirthYear +MumCode,
                 family="threshold",
                 data=juvenile_surv_df_na_rm,
                 prior = prior,
                 pr=TRUE,#saves posterior dist for random effects i.e. what we want
                 nitt=200000,burnin=50000
                 )##



# 200k finished in <20 mins 


plot(juvenile_surv_model)

#trace plots look good .. mum code a bit dodge?

summary(juvenile_surv_model)


save(juvenile_surv_model, file="PhD_4th_yr/Inbreeding_depression_models/survival/juvenile_survival_model_output_200kit_50kbu_incBW.RData")

