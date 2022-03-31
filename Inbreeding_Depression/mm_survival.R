library(RODBC)
library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)
library("purrr")
library("MCMCglmm")
library("data.table")
library(tibble)
library(ggplot2)
library(MasterBayes)

###########################################################################
### reading pedigree etc from database ####################################
###########################################################################


db<-"C:\\Users\\s1881212\\Documents\\Deer_database_08_2021\\RedDeer1.99.accdb" #open connection
con<-odbcConnectAccess2007(db)
life<-sqlFetch(con, "tbllife") %>% 
  dplyr::select(Code, BirthDay,BirthYear, BirthMonth, BirthYear,DeathDay,DeathMonth,DeathYear, Sex, MumCode, DeathType)%>%
  mutate(age_at_death=DeathYear-BirthYear)%>%
  filter(DeathType!= c("S","A"))
  

# getting table life with birth yr etc
mum_stat<-sqlFetch(con, "sys_HindStatusAtConception") %>% dplyr::select(-CalfBirthYear)%>%dplyr::rename(MumCode=Mum, Code=Calf)
# also getting the status of the hind when calf born
#remember need the R-friendly button on!!###
ped<-sqlFetch(con, "sys_Pedigree")%>%dplyr::select(Code, MumCode,Sire)


ped$Code <- as.factor(ped$Code)
ped$MumCode <- as.factor(ped$MumCode)
ped$Sire <- as.factor(ped$Sire)

odbcClose(con)


## yearling survival i.e. up to May 1st of second yr
life$yearling_surv<-NA
life$yearling_surv<-if_else(life$age_at_death==0|life$age_at_death==1, 0,
                            if_else(life$age_at_death==2 & life$DeathMonth < 5,0, 1))

#if ids die over winter (Oct 1st - May 1st)                       
life$winter_surv<-NA
life$winter_surv<-if_else((life$age_at_death==0 & life$DeathMonth<5)|#dying in first few months of year less than a yr old
                            (life$age_at_death==0 & life$DeathMonth>=5 & life$BirthMonth<=life$DeathMonth), 0, 1)#dying between Summer of of yr born and end of yr

#ids that die before Oct 1st in their first yr
life$neonatal_survival<-NA
life$neonatal_survival<-if_else(life$age_at_death==0 & life$DeathMonth<10 & life$BirthMonth<=life$DeathMonth, 0, 1)



#set up df

setwd("H:/")

###need to make this with updated map positions
FROH<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,CHR,KB) %>% dplyr::rename(Code=IID)

KB_perLG<-FROH%>%dplyr::group_by(Code, CHR)%>%
  dplyr::summarise(KB_chr=sum(KB))%>%
  ungroup %>% 
  complete(Code, CHR, fill = list(KB_chr = 0)) #completes for all chromosomes 


deermap <- read.table("PhD_3rdYR/Data_files/TableS1_CervusElaphus_Final_Linkage_Map.txt", header = T, stringsAsFactors = F)
deersum <- deermap %>% 
  dplyr::group_by(CEL.LG) %>%
  dplyr::summarise(Max_EstMb = max(Estimated.Mb.Position))%>%
  dplyr::mutate(chr_length_KB=Max_EstMb/1000)%>%
  dplyr::rename(CHR=CEL.LG)


froh_per_chr<-join(KB_perLG,deersum)%>% mutate(chr_froh=KB_chr/chr_length_KB)%>%
  dplyr::select(-Max_EstMb)%>% reshape2::dcast(Code~CHR) %>% mutate(FROHsum = rowSums(.[2:34]))%>%mutate(FROH_sum_div=FROHsum/33)


colnames(froh_per_chr) <- c("Code", paste0("FROH_chr", 1:33),"FROHsum","FROH_sum_div")

FROH_full<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom.indiv", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,KB) %>% dplyr::rename(Code=IID)%>%mutate(FROH=KB/2495700)



surv_df<-life%>%join(mum_stat)%>%join(froh_per_chr)%>%
  join(FROH_full)%>%dplyr::select(-BirthDay,-DeathDay,-DeathMonth,-BirthMonth)%>%na.omit()

#### first run normal model leaving one chr out then try mm ###

# prepping and inverting pedigree
ids_in<-as.matrix(surv_df%>%dplyr::select(Code))
pruned<-prunePed(ped, ids_in)
ped_ordered<-orderPed(pruned)
ped_ordered$Code <- as.factor(ped_ordered$Code)
ped_ordered$MumCode <- as.factor(ped_ordered$MumCode)
ped_ordered$Sire <- as.factor(ped_ordered$Sire)

Ainv<-inverseA(ped_ordered, nodes="ALL")$Ainv


save(surv_data,Ainv, file="PhD_3rdYR/Model_inputs/df_survival_forMCMCglmm.RData")

#####################################################################
##### this will now all be done on ash server #####################
###################################################################

k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=0.25,nu=0.002,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k), ## multimemberhsip part
                   G1=list(V=0.25,nu=0.002,aplha.mu=0,alpha.V=k),
                   G1=list(V=0.25,nu=0.002,aplha.mu=0,alpha.V=k)))

#### RUNNING MULTI-MEMBERSHIP MODEL ####
# changing FROH sum to FROH or FROH sum/33 is same really 
model2<-MCMCglmm(yearling_surv~1 + Sex + MotherStatus + FROHsum, #need to fit sum chrFROH  as continuous covariate,
                random= ~ Code +
                  idv(FROH_chr1+FROH_chr2+FROH_chr3+FROH_chr4+FROH_chr5+FROH_chr6+FROH_chr7+FROH_chr8+
                        FROH_chr9+FROH_chr10+FROH_chr11+FROH_chr12+FROH_chr13+FROH_chr14+FROH_chr15+FROH_chr16+
                        FROH_chr17+FROH_chr18+FROH_chr19+FROH_chr20+FROH_chr21+FROH_chr22+FROH_chr23+FROH_chr24+
                        FROH_chr25+FROH_chr26+FROH_chr27+FROH_chr28+FROH_chr29+FROH_chr30+FROH_chr31+FROH_chr32+
                        FROH_chr33)+BirthYear +MumCode,
                family="categorical",
                data=surv_df,
                ginverse = list(Code=Ainv), # code is recognised as the pedigree 
                prior = prior,
                pr=TRUE,#saves posterior dist for random effects i.e. what we want 
                nitt=1000000,burnin=1500000, thin = 200)

#save(model2, file="PhD_3rdYR/Model outputs/mm_survival_17.03.RData")

summary(model2)
#plot(model$Sol)

plot(model2)


names <- apply(model2$Sol,2,mean) %>% names ## gets names of all random variables 
sols<-apply(model2$Sol,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(model2$Sol,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
CI_lower<-apply(model2$Sol,2,quantile,probs = c(0.025)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH

FROH_sols$CHR<-as.factor(FROH_sols$CHR)

ggplot(data=FROH_sols, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Linkage group", y="solution + CI", title = "Effect size of Linkage group FROH on Birthweight (kg)") +
  theme_bw()+  # use a white background                 
  theme(legend.position = "none")

#### saving plot in wd
ggsave(file="test_gg.png",
       plot = ggplot( data=FROH_sols, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Linkage group", y="solution + CI", title = "Effect size of Linkage group FROH on Birthweight (kg)") +
  theme_bw()+  # use a white background                 
  theme(legend.position = "none"))

  
  
  
  





### saving trace plots 
library(lattice)

xyplot(
  as.mcmc(model2$VCV),  ## convert matrix to `mcmc` object
  layout=c(2,1)            ## customize panel arrangement
)

## help from https://bbolker.github.io/morelia_2018/notes/bayeslab.html
## variance components 
jpeg(file="test2.jpeg")
xyplot(
  as.mcmc(model2$VCV)  ## convert matrix to `mcmc` object
)
dev.off()

#fixed effects plots 
jpeg(file="test3.jpeg")
xyplot(
  as.mcmc((model2$Sol)[,c(1:4)])  ## convert matrix to `mcmc` object
  ## customize panel arrangement
)
dev.off()

jpeg(file="test4.jpeg")
xyplot(
  as.mcmc((model2$Sol)[,c(5:7)])  ## convert matrix to `mcmc` object
  ## customize panel arrangement
)
dev.off()
