
library(MCMCglmm)
library(MasterBayes)
library(data.table)

setwd("H:/")

female_LBS_df=read.table("PhD_4th_yr/Inbreeding_depression_models/LBS/female_LBS_df.txt", sep=",", header=T)

female_LBS_df$Code=as.factor(female_LBS_df$Code)
female_LBS_df$BirthYear=as.factor(female_LBS_df$BirthYear)
female_LBS_df$MumCode=as.factor(female_LBS_df$MumCode)


k<-10000
prior.2 = list(R = list(V = diag(2), nu = 2,fix = 2),
               G = list(G1 = list(V = diag(2), nu = 2),
                        G2 = list(V = diag(2), nu = 2),
                        G3 = list(V=1,nu=1,alpha.mu=0,alpha.V=k)))#G3 referring to chrFROH
#nu=0

female_LBS_model=MCMCglmm(LBS~ trait+trait:FROHsum-1,
                          
                          random=~idh(trait):MumCode+idh(trait):BirthYear + 
                            idv(FROH_chr1+ FROH_chr2+FROH_chr3+FROH_chr4+FROH_chr5+FROH_chr6+FROH_chr7+FROH_chr8+
                                  FROH_chr9+FROH_chr10+FROH_chr11+FROH_chr12+FROH_chr13+FROH_chr14+FROH_chr15+FROH_chr16+
                                  FROH_chr17+FROH_chr18+FROH_chr19+FROH_chr20+FROH_chr21+FROH_chr22+FROH_chr23+FROH_chr24+
                                  FROH_chr25+FROH_chr26+FROH_chr27+FROH_chr28+FROH_chr29+FROH_chr30+FROH_chr31+FROH_chr32+
                                  FROH_chr33), 
                          rcov = ~ idh(trait):units,
                          data = female_LBS_df, 
                          family = "zipoisson", 
                          prior = prior.2,
                          #verbose = FALSE, 
                          pr = TRUE, 
                          pl = TRUE,
                          nitt=1000000, burnin=100000, thin = 500
)




plot(female_LBS_model)


summary(female_LBS_model)



save(female_LBS_model, file="PhD_4th_yr/Inbreeding_depression_models/LBS/Female_LBS_model_output_1mit.RData")
