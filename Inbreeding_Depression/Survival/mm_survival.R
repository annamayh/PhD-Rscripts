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


db<-"C:\\Users\\s1881212\\Documents\\Deer_database_08_2021\\RedDeer1.99.accdb" #open connection
con<-odbcConnectAccess2007(db)

mum_birthyr<-sqlFetch(con, "tbllife")%>%dplyr::select(Code,BirthYear)%>%dplyr::rename(mum_birthyear=BirthYear,MumCode=Code)
mum_stat<-sqlFetch(con, "sys_HindStatusAtConception") %>% dplyr::select(-CalfBirthYear)%>%dplyr::rename(MumCode=Mum, Code=Calf)
birth_wt<-sqlFetch(con, "sys_BirthWt")%>%dplyr::select(BirthWt,Code)

life<-sqlFetch(con, "tbllife") %>% 
  dplyr::select(Code, BirthDay,BirthYear, BirthMonth, BirthYear,DeathDay,DeathMonth,DeathYear, Sex, MumCode, DeathType)%>%
  filter(DeathType!= c("S","A","D"))%>%
  drop_na(BirthYear,DeathYear,DeathMonth,BirthMonth)%>%
  join(mum_birthyr)%>%mutate(mum_age=BirthYear-mum_birthyear)%>%
  join(mum_stat)%>%join(birth_wt)

ped<-sqlFetch(con, "sys_Pedigree")%>%dplyr::select(Code, MumCode,Sire)


ped$Code <- as.factor(ped$Code)
ped$MumCode <- as.factor(ped$MumCode)
ped$Sire <- as.factor(ped$Sire)

odbcClose(con)

########################
### neonatl survival ###
#######################

#ids that die before Oct 1st in their first yr
neo<-life%>%mutate(mum_age_sq=mum_age^2)

neo$neonatal_survival<-NA
neo$neonatal_survival<-if_else(neo$DeathMonth<10 & neo$DeathYear==neo$BirthYear, 0, 1)


########################
### winter survival ####
########################
winter<-neo%>%filter(neonatal_survival==1) #removing 714 ids that died as neonates
#if ids die over winter (Oct 1st - May 1st)                       
winter$winter_surv<-NA
winter$winter_surv<-if_else((winter$DeathYear==winter$BirthYear|
                               winter$DeathYear==(winter$BirthYear)+1&winter$DeathMonth<5), 0, 1)#dying between Summer of of yr born and end of yr



########################
## yearling survival i.e. up to May 1st of second yr
###########################
year<-winter%>%filter(winter_surv==1) #removing  ids that die over winter

year$yearling_surv<-NA
year$yearling_surv<-if_else(year$DeathYear==(year$BirthYear)+2 & year$DeathMonth<5|
                              year$DeathYear==(year$BirthYear)+1 & year$DeathMonth>5, 0, 1)

##########################
## juvenile survival #####
##########################

juv<-year%>%filter(yearling_surv==1)%>%dplyr::select(-yearling_surv,-winter_surv,-neonatal_survival)## removing 864 ids that died before first yr
juv$juvenile_surv<-1

juv_all<-life%>%join(juv)%>%mutate_at(vars("juvenile_surv"),~replace_na(.,0))%>%mutate(mum_age_sq=mum_age^2)


#set up df

setwd("H:/")

###need to make this with updated map positions
FROH<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,CHR,KB) %>% dplyr::rename(Code=IID)

KB_perLG<-FROH%>%dplyr::group_by(Code, CHR)%>%
  dplyr::summarise(KB_chr=sum(KB))%>%
  ungroup %>% 
  complete(Code, CHR, fill = list(KB_chr = 0)) #completes for all chromosomes 


deermap <- read.csv("PhD_3rdYR/Data_files/Genome_assembly_mCerEla1.1.csv", header = T, stringsAsFactors = F)%>%
  dplyr::filter(!CHR %in% c("All","All_auto","X","unplaced"))

froh_per_chr<-join(KB_perLG,deermap)%>% mutate(chr_froh=KB_chr/length_Kb)%>%
  dplyr::select(-length_Mb,-length)%>% reshape2::dcast(Code~CHR) %>% mutate(FROHsum = rowSums(.[2:34]))%>%mutate(FROH_sum_div=FROHsum/33)


colnames(froh_per_chr) <- c("Code", paste0("FROH_chr", 1:33),"FROHsum","FROH_sum_div")

FROH_full<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom.indiv", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,KB) %>% dplyr::rename(Code=IID)%>%mutate(FROH=KB/2591865)

#####################################################################################################################

##### set up dfs for models ###


#neonatal_survival
neonat_df<-neo%>%join(froh_per_chr)%>%
  join(FROH_full)%>%dplyr::select(-BirthDay,-DeathDay,-DeathMonth,-BirthMonth)%>%na.omit()

## winter survival
winter_df<-winter%>%join(froh_per_chr)%>%
  join(FROH_full)%>%dplyr::select(-BirthDay,-DeathDay,-DeathMonth,-BirthMonth)%>%na.omit()

#yearling survival
year_df<-year%>%join(froh_per_chr)%>%
  join(FROH_full)%>%dplyr::select(-BirthDay,-DeathDay,-DeathMonth,-BirthMonth)%>%na.omit()

#juvenile survival
juve_df<-juv_all%>%join(froh_per_chr)%>%
  join(FROH_full)%>%dplyr::select(-BirthDay,-DeathDay,-DeathMonth,-BirthMonth)%>%na.omit()



# prepping and inverting pedigree
#ids_in<-as.matrix(surv_df%>%dplyr::select(Code))
#pruned<-prunePed(ped, ids_in)
ped_ordered<-orderPed(ped)
ped_ordered$Code <- as.factor(ped_ordered$Code)
ped_ordered$MumCode <- as.factor(ped_ordered$MumCode)
ped_ordered$Sire <- as.factor(ped_ordered$Sire)

Ainv<-inverseA(ped_ordered, nodes="ALL")$Ainv


#save(neonat_df,winter_df,year_df,juve_df,Ainv, file="PhD_3rdYR/Model_inputs/df_ALLsurvival_forMCMCglmm.RData")

setwd("H:/")
load(file="PhD_3rdYR/Model_inputs/df_ALLsurvival_forMCMCglmm.RData")
#####################################################################################################################################
# 
k<-100
prior<-list(R=list(V=1,fix=1),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k), ## multimemberhsip part
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#### RUNNING MULTI-MEMBERSHIP MODEL ####
# changing FROH sum to FROH or FROH sum/33 is same really
model3<-MCMCglmm(juvenile_surv~1 + Sex + MotherStatus + FROHsum + mum_age+mum_age_sq, #need to fit sum chrFROH  as continuous covariate,
                random= ~ 
                  idv(FROH_chr1+FROH_chr2+FROH_chr3+FROH_chr4+FROH_chr5+FROH_chr6+FROH_chr7+FROH_chr8+
                        FROH_chr9+FROH_chr10+FROH_chr11+FROH_chr12+FROH_chr13+FROH_chr14+FROH_chr15+FROH_chr16+
                        FROH_chr17+FROH_chr18+FROH_chr19+FROH_chr20+FROH_chr21+FROH_chr22+FROH_chr23+FROH_chr24+
                        FROH_chr25+FROH_chr26+FROH_chr27+FROH_chr28+FROH_chr29+FROH_chr30+FROH_chr31+FROH_chr32+
                        FROH_chr33)+BirthYear +MumCode,
                family="threshold",
                data=juve_df,
                prior = prior,
                pr=TRUE,#saves posterior dist for random effects i.e. what we want
                nitt=150000,burnin=50000)##

save(model3, file="PhD_3rdYR/Model outputs/Juvenile_survival_0-2/mm_juvenile2.RData")


#maybe try without birthwt?


# 
# setwd("H:/")
load(file="PhD_3rdYR/Model outputs/Juvenile_survival_0-2/mm_juvenile2.RData")
# 
# 
summary(model3)
# #plot(model$Sol)

plot(model3)
#

model<-model3

names <- apply(model2$Sol,2,mean) %>% names ## gets names of all random variables
sols<-apply(model2$Sol,2,mean)#gets mean of all solutions i.e. the effect size of random effects
CI_upper<-apply(model2$Sol,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions
CI_lower<-apply(model2$Sol,2,quantile,probs = c(0.025)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH

FROH_sols$CHR<-as.factor(FROH_sols$CHR)
#
ggplot(data=FROH_sols, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="solution + CI", title = "Effect of Chr FROH on juvenile survival (no ped, no BthWt)") +
  theme_bw()+  # use a white background
  theme(legend.position = "none")

###################################################################################################################
###################################################################################################################


sols_full<-as.data.frame(model$Sol)%>%dplyr::select(matches("FROH"))%>% ## taking out sols with FROH included
  dplyr::mutate(across(2:34, ~.x + FROHsum)) ## adding FROHsum to chrFROH values

names <- apply(sols_full,2,mean) %>% names ## gets names of all random variables, 2 = all down row
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
  geom_hline(yintercept=-0.216671, lty=1,colour="red") + ## red line is the average effect of all chromosomes 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="solution + CI", title = "Effect size of chromosome FROH on survival") +
  theme_bw()+  # use a white background                 
  theme(legend.position = "none")


sum(FROH_sols$solution)

#### saving plot in wd
# ggsave(file="test_gg.png",
#        plot = ggplot( data=FROH_sols, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper)) +
#   geom_pointrange() + #plots lines based on Y and lower and upper CI
#   geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
#   coord_flip() +  # flip coordinates (puts labels on y axis)
#   labs(x="Linkage group", y="solution + CI", title = "Effect size of Linkage group FROH on Birthweight (kg)") +
#   theme_bw()+  # use a white background
#   theme(legend.position = "none"))
#
#
#

  





### saving trace plots 
# library(lattice)
# 
# xyplot(
#   as.mcmc(model2$VCV),  ## convert matrix to `mcmc` object
#   layout=c(2,1)            ## customize panel arrangement
# )
# 
# ## help from https://bbolker.github.io/morelia_2018/notes/bayeslab.html
# ## variance components 
# jpeg(file="test2.jpeg")
# xyplot(
#   as.mcmc(model2$VCV)  ## convert matrix to `mcmc` object
# )
# dev.off()
# 
# #fixed effects plots 
# jpeg(file="test3.jpeg")
# xyplot(
#   as.mcmc((model2$Sol)[,c(1:4)])  ## convert matrix to `mcmc` object
#   ## customize panel arrangement
# )
# dev.off()
# 
# jpeg(file="test4.jpeg")
# xyplot(
#   as.mcmc((model2$Sol)[,c(5:7)])  ## convert matrix to `mcmc` object
#   ## customize panel arrangement
# )
# dev.off()
