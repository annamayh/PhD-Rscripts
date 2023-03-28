library(RODBC)
library(tidyverse)
library(MCMCglmm)
library(MasterBayes)
library(data.table)


db<-"C:\\Users\\s1881212\\Documents\\Deer_database_08_2021\\RedDeer1.99.accdb" #open connection
con<-odbcConnectAccess2007(db)

female_LBS<-sqlFetch(con, "sys_LBS+LRS") #built in query - Calculate LBS and LRS for dead hinds

life<-sqlFetch(con, "tbllife") %>% 
  dplyr::select(Code, BirthYear, MumCode)

odbcClose(con)



setwd("H:/")

FROH<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,CHR,KB) %>% dplyr::rename(Code=IID)

KB_perLG<-FROH%>%dplyr::group_by(Code, CHR)%>%
  dplyr::summarise(KB_chr=sum(KB))%>%
  ungroup %>% 
  complete(Code, CHR, fill = list(KB_chr = 0)) #completes for all chromosomes 

deermap <- read.csv("PhD_3rdYR/Data_files/Genome_assembly_mCerEla1.1.csv", header = T, stringsAsFactors = F)%>%
  dplyr::filter(!CHR %in% c("All","All_auto","X","unplaced"))

froh_per_chr<-plyr::join(KB_perLG,deermap)%>% mutate(chr_froh=KB_chr/length_Kb)%>%
  dplyr::select(-length_Mb,-length)%>% reshape2::dcast(Code~CHR) %>% mutate(FROHsum = rowSums(.[2:34]))%>%mutate(FROH_sum_div=FROHsum/33)

colnames(froh_per_chr) <- c("Code", paste0("FROH_chr", 1:33),"FROHsum","FROH_sum_div")









female_LBS_df=female_LBS%>%
  plyr::join(froh_per_chr)%>%plyr::join(life)%>%
  na.omit()




k<-10000
prior.2 = list(R = list(V = diag(2), nu = 2,fix = 2),
               G = list(G1 = list(V = diag(2), nu = 2),
                        G2 = list(V = diag(2), nu = 2),
                        G3 = list(V=1,nu=1,aplha.mu=0,alpha.V=k)))#G3 referring to chrFROH


female_LBS_df$Code=as.factor(female_LBS_df$Code)
female_LBS_df$BirthYear=as.factor(female_LBS_df$BirthYear)
female_LBS_df$MumCode=as.factor(female_LBS_df$MumCode)




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

# R-structure:  ~idh(trait):units
# 
# post.mean l-95% CI u-95% CI eff.samp
# traitLBS.units       0.1256  0.08355   0.1766     1589
# traitzi_LBS.units    1.0000  1.00000   1.0000        0
# 
# Location effects: LBS ~ trait + trait:FROHsum - 1 
# 
# post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# traitLBS              1.36751  1.13945  1.61879     1969 <6e-04 ***
#   traitzi_LBS          -1.52672 -2.16207 -0.78001     1800 <6e-04 ***
#   traitLBS:FROHsum     -0.02409 -0.11718  0.06927     1800  0.609    
# traitzi_LBS:FROHsum   0.42917  0.22179  0.68191     1800 <6e-04 ***
# 



save(female_LBS_model, file = "PhD_4th_yr/Inbreeding_depression_models/LBS/OUTPUT_Female_LBS_mm_ZIPOIS_1mit.RData")

load(file = "PhD_4th_yr/Inbreeding_depression_models/LBS/OUTPUT_Female_LBS_mm_ZIPOIS_1mit.RData")


summary(female_LBS_model)

sols_full<-as.data.frame(female_LBS_model$Sol)%>%dplyr::select(matches("FROH"))%>% ## taking out sols with FROH included
  rename(LBS_zi="traitzi_LBS:FROHsum")%>%
  dplyr::mutate(across(3:35, ~.x + LBS_zi)) ## adding FROHsum to chrFROH values

names <- apply(sols_full,2,mean) %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_full,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(sols_full,2,quantile,probs = c(0.95)) #gets upper confidence interval for all solutions 
CI_lower<-apply(sols_full,2,quantile,probs = c(0.05)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH




## Plot chromosomal inbreeding effects on male LBS ##
#Red line is average chromosomal effect
#dashed line is 0 effect 

LBS_froh_all=sum(FROH_sols$solution) 
sum(FROH_sols$solution) 


FROH_sols$CHR<-as.factor(FROH_sols$CHR)

ggplot(data=FROH_sols, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange(colour="grey50") + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  geom_hline(yintercept=0.42917, lty=1,colour="red") + ## red line is the average effect of all chromosomes 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="solution + CI", title = "Effect of chromosome FROH on female LBS poisson") +
  theme_bw()+  # use a white background                 
  theme(legend.position = "none")








###################################################################################################


female_LBS_model_hurdle=MCMCglmm(LBS~ trait+trait:FROHsum-1,
                          
                          random=~idh(trait):MumCode+idh(trait):BirthYear + 
                            idv(FROH_chr1+ FROH_chr2+FROH_chr3+FROH_chr4+FROH_chr5+FROH_chr6+FROH_chr7+FROH_chr8+
                                  FROH_chr9+FROH_chr10+FROH_chr11+FROH_chr12+FROH_chr13+FROH_chr14+FROH_chr15+FROH_chr16+
                                  FROH_chr17+FROH_chr18+FROH_chr19+FROH_chr20+FROH_chr21+FROH_chr22+FROH_chr23+FROH_chr24+
                                  FROH_chr25+FROH_chr26+FROH_chr27+FROH_chr28+FROH_chr29+FROH_chr30+FROH_chr31+FROH_chr32+
                                  FROH_chr33), 
                          rcov = ~ idh(trait):units,
                          data = female_LBS_df, 
                          family = "hupoisson", 
                          prior = prior.2,
                          #verbose = FALSE, 
                          pr = TRUE, 
                          pl = TRUE,
                          nitt=500000, burnin=50000, thin = 500
)


##need to runf for 1mil 

plot(female_LBS_model_hurdle)


summary(female_LBS_model_hurdle)


plogis(0.36699)

# Location effects: LBS ~ trait + trait:FROHsum - 1 
# 
# post.mean l-95% CI u-95% CI eff.samp   pMCMC   
# traitLBS              1.43530  1.20037  1.65246    900.0 < 0.001 **
#   traithu_LBS          -1.11606 -1.79453 -0.53781    965.5 0.00222 **
#   traitLBS:FROHsum     -0.03408 -0.12990  0.05238    900.0 0.53333   
# traithu_LBS:FROHsum   0.36699  0.14284  0.57011    900.0 < 0.001 **
