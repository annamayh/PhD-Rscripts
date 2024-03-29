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
birth_wt<-sqlFetch(con, "sys_BirthWt") 
# getting trait 
life<-sqlFetch(con, "tbllife") %>% dplyr::select(Code, BirthYear, MumCode, Sex) 
# getting table life with birth yr etc
mum_stat<-sqlFetch(con, "sys_HindStatusAtConception") %>% dplyr::select(-CalfBirthYear)%>%dplyr::rename(MumCode=Mum, Code=Calf)
# also getting the status of the hind when calf born
#remember need the R-friendly button on!!###
ped<-sqlFetch(con, "sys_Pedigree")%>%dplyr::select(Code, MumCode,Sire)
ped$Code <- as.factor(ped$Code)
ped$MumCode <- as.factor(ped$MumCode)
ped$Sire <- as.factor(ped$Sire)

## ordering and inverting pedigree



odbcClose(con) #close connection

################################################
### read in ROH file to work out LG FROH ###
#################################################
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

#FROH per chr each id
froh_per_chr<-join(KB_perLG,deermap)%>% mutate(chr_froh=KB_chr/length_Kb)%>%
  dplyr::select(-length_Mb,-length)%>% reshape2::dcast(Code~CHR) %>% mutate(FROHsum = rowSums(.[2:34]))%>%mutate(FROH_sum_div=FROHsum/33)

colnames(froh_per_chr) <- c("Code", paste0("FROH_chr", 1:33),"FROHsum","FROH_sum_div")

#FROH per chr mother 

mum_froh_chr=froh_per_chr
colnames(mum_froh_chr) <- c("MumCode", paste0("Mum_FROH_chr", 1:33),"Mum_FROHsum","Mum_FROH_sum_div")


# FROH_full<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom.indiv", header=T, stringsAsFactors = F)%>%
#   dplyr::select(IID,KB) %>% dplyr::rename(Code=IID)%>%mutate(FROH=KB/2591865)


Birth_wt_df<-birth_wt%>%join(life)%>%join(mum_stat)%>%join(froh_per_chr)%>%join(mum_froh_chr)%>%na.omit()

#need to set factors for categorical variables#
Birth_wt_df$Sex <- factor(Birth_wt_df$Sex)
Birth_wt_df$Code <- factor(Birth_wt_df$Code)
Birth_wt_df$MotherStatus <- factor(Birth_wt_df$MotherStatus)
Birth_wt_df$MumCode <- factor(Birth_wt_df$MumCode)

### sorting out pedigree

ids_in<-as.matrix(Birth_wt_df%>%dplyr::select(Code))

pruned<-prunePed(ped, ids_in)
ped_ordered<-orderPed(pruned)
ped_ordered$Code <- as.factor(ped_ordered$Code)
ped_ordered$MumCode <- as.factor(ped_ordered$MumCode)
ped_ordered$Sire <- as.factor(ped_ordered$Sire)

Ainv<-inverseA(ped_ordered, nodes="ALL")$Ainv


## set priors
# kirstys priors to see if they help
k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#### RUNNING MULTI-MEMBERSHIP MODEL ####
# changing FROH sum to FROH or FROH sum/33 is same really 
model<-MCMCglmm(CaptureWt~1 + Sex + AgeHrs + MotherStatus + FROHsum +Mum_FROHsum, #need to fit sum chrFROH  as continuous covariate,
                random= ~ Code +
                  idv(FROH_chr1+FROH_chr2+FROH_chr3+FROH_chr4+FROH_chr5+FROH_chr6+FROH_chr7+FROH_chr8+
                        FROH_chr9+FROH_chr10+FROH_chr11+FROH_chr12+FROH_chr13+FROH_chr14+FROH_chr15+FROH_chr16+
                        FROH_chr17+FROH_chr18+FROH_chr19+FROH_chr20+FROH_chr21+FROH_chr22+FROH_chr23+FROH_chr24+
                        FROH_chr25+FROH_chr26+FROH_chr27+FROH_chr28+FROH_chr29+FROH_chr30+FROH_chr31+FROH_chr32+
                        FROH_chr33) +
                  idv(Mum_FROH_chr1+Mum_FROH_chr2+Mum_FROH_chr3+Mum_FROH_chr4+Mum_FROH_chr5+Mum_FROH_chr6+Mum_FROH_chr7+Mum_FROH_chr8+
                        Mum_FROH_chr9+Mum_FROH_chr10+Mum_FROH_chr11+Mum_FROH_chr12+Mum_FROH_chr13+Mum_FROH_chr14+Mum_FROH_chr15+Mum_FROH_chr16+
                        Mum_FROH_chr17+Mum_FROH_chr18+Mum_FROH_chr19+Mum_FROH_chr20+Mum_FROH_chr21+Mum_FROH_chr22+Mum_FROH_chr23+Mum_FROH_chr24+
                        Mum_FROH_chr25+Mum_FROH_chr26+Mum_FROH_chr27+Mum_FROH_chr28+Mum_FROH_chr29+Mum_FROH_chr30+Mum_FROH_chr31+Mum_FROH_chr32+
                        Mum_FROH_chr33)+BirthYear +MumCode,
                data=Birth_wt_df,
                ginverse = list(Code=Ainv), # code is recognised as the pedigree 
                prior = prior,
                pr=TRUE,#saves posterior dist for random effects i.e. what we want 
                nitt=50000,burnin=5000)


#save(model, file="PhD_3rdYR/Model outputs/Birth_weight/mm_birthwt_interpretable.RData")


###############################################################################################
###############################################################################################
######### LOAD MODEL ALREADY RUN ##############################################################
################################################################################################



load("PhD_3rdYR/Model outputs/Birth_weight/mm_birthwt.RData")
#### getting output info ###
# 
summary(model)

xyplot(
  as.mcmc((model2$Sol)[,c(1:4)])  
)

plot(model$Sol)
# 
plot(model)
# 
# names <- apply(model$Sol,2,mean) %>% names ## gets names of all random variables, 2 = all down row
# sols<-apply(model$Sol,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
# CI_upper<-apply(model$Sol,2,quantile,probs = c(0.95)) #gets upper confidence interval for all solutions 
# CI_lower<-apply(model$Sol,2,quantile,probs = c(0.05)) #gets lower CI for all solutions
# 
# Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)
# 
# names(Random_table)[1]<-"solution"
# names(Random_table)[2]<-"model_variable"
# 
# FROH_sols<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH
# 
# FROH_sols$CHR<-as.factor(FROH_sols$CHR)
# 
# ggplot(data=FROH_sols, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper)) +
#   geom_pointrange(colour="grey50") + #plots lines based on Y and lower and upper CI
#   geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
#   coord_flip() +  # flip coordinates (puts labels on y axis)
#   labs(x="Linkage group", y="solution + CI", title = "Effect size of Linkage group FROH on Birthweight (kg)") +
#   theme_bw()+  # use a white background                 
#   theme(legend.position = "none")
# 

#####################################################################################################################
###### getting in terms of other bw model ##########################################################################
###################################################################################################################

sols_full<-as.data.frame(model$Sol)%>%dplyr::select(matches("FROH"))%>% ## taking out sols with FROH included
  dplyr::mutate(across(3:35, ~.x + FROHsum))%>%
  dplyr::mutate(across(36:68, ~.x + Mum_FROHsum))## adding FROHsum to chrFROH values

names <- apply(sols_full,2,mean) %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_full,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(sols_full,2,quantile,probs = c(0.95)) #gets upper confidence interval for all solutions 
CI_lower<-apply(sols_full,2,quantile,probs = c(0.05)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols<-Random_table%>%filter(model_variable %like% "FROH_c" & !model_variable %like% "Mum" )%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH

FROH_sols$CHR<-as.factor(FROH_sols$CHR)

ggplot(data=FROH_sols, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange(colour="grey50") + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  geom_hline(yintercept=-0.14, lty=1,colour="red") + ## red line is the average effect of all chromosomes 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="solution + CI", title = "Effect size of chromosome FROH on Birthweight (kg)") +
  theme_bw()+  # use a white background                 
  theme(legend.position = "none")


sum(FROH_sols$solution) ## givs the total effect of all chromosomes shown to be 4.3kg


## for mum froh


FROH_sols_mum<-Random_table%>%filter(model_variable %like% "Mum_FROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH

FROH_sols_mum$CHR<-as.factor(FROH_sols_mum$CHR)

ggplot(data=FROH_sols_mum, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange(colour="grey50") + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  geom_hline(yintercept=-0.027, lty=1,colour="red") + ## red line is the average effect of all chromosomes 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="solution + CI", title = "Effect size of chromosome FROH on Birthweight (kg)") +
  theme_bw()+  # use a white background                 
  theme(legend.position = "none")



