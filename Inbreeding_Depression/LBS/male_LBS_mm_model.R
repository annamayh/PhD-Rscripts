
## MALE LBS ###
library(RODBC)
library(tidyverse)
library(MCMCglmm)
library(MasterBayes)

# Set up dataframe #

db<-"C:\\Users\\s1881212\\Documents\\Deer_database_08_2021\\RedDeer1.99.accdb" #open connection
con<-odbcConnectAccess2007(db)
ped<-sqlFetch(con, "sys_Pedigree")%>%dplyr::select(Code, BirthYear,MumCode,Sire)
life<-sqlFetch(con, "tbllife") %>% 
  dplyr::select(Code,DeathYear,Sex)
odbcClose(con)

sire_count=ped%>%dplyr::select(Code,Sire)%>%
  group_by(Sire)%>%tally()%>% ##grouping by sire and counting up ids sired by id
  rename(Code=Sire, LBS=n)

dead_males=life%>%
  filter(Sex=="2")%>% #select males
  na.omit() ##removing males that are still alive (i.e. no death recorded)

all_male_LBS=left_join(dead_males,sire_count)
all_male_LBS[is.na(all_male_LBS)] <- 0 #males with no recorded sired ids given 0 LBS


# Read in ROH per chr #

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

froh_per_chr<-plyr::join(KB_perLG,deermap)%>% mutate(chr_froh=KB_chr/length_Kb)%>%
  dplyr::select(-length_Mb,-length)%>% reshape2::dcast(Code~CHR) %>% mutate(FROHsum = rowSums(.[2:34]))%>%mutate(FROH_sum_div=FROHsum/33)

colnames(froh_per_chr) <- c("Code", paste0("FROH_chr", 1:33),"FROHsum","FROH_sum_div")

# Combine dataframes #

male_LBS_df=all_male_LBS%>%
  plyr::join(froh_per_chr)%>%plyr::join(ped)%>%
  dplyr::select(-Sire, -DeathYear, -Sex)%>%
  na.omit()

male_LBS_df$Code=as.factor(male_LBS_df$Code)
male_LBS_df$BirthYear=as.factor(male_LBS_df$BirthYear)
male_LBS_df$MumCode=as.factor(male_LBS_df$MumCode)


## priors 
k<-10000
prior.2 = list(R = list(V = diag(2), nu = 2,fix = 2),
               G = list(G1 = list(V = diag(2), nu = 2),
                        G2 = list(V = diag(2), nu = 2),
                        G3 = list(V=1,nu=1,aplha.mu=0,alpha.V=k)))#G3 referring to chrFROH


#################################################################################
# Run model #
#################################################################################
male_LBS_model=MCMCglmm(LBS~ trait+trait:FROHsum-1,

                        random=~idh(trait):MumCode+idh(trait):BirthYear + 
                          idv(FROH_chr1+ FROH_chr2+FROH_chr3+FROH_chr4+FROH_chr5+FROH_chr6+FROH_chr7+FROH_chr8+
                                FROH_chr9+FROH_chr10+FROH_chr11+FROH_chr12+FROH_chr13+FROH_chr14+FROH_chr15+FROH_chr16+
                                FROH_chr17+FROH_chr18+FROH_chr19+FROH_chr20+FROH_chr21+FROH_chr22+FROH_chr23+FROH_chr24+
                                FROH_chr25+FROH_chr26+FROH_chr27+FROH_chr28+FROH_chr29+FROH_chr30+FROH_chr31+FROH_chr32+
                                FROH_chr33), 
                        rcov = ~ idh(trait):units,
                        data = male_LBS_df, 
                        family = "zipoisson", 
                        prior = prior.2,
                        #verbose = FALSE, 
                        pr = TRUE, 
                        pl = TRUE,
                        nitt=1500000, burnin=150000, thin = 500
)


#first fixed term = intercept of the Poisson process 
#second fixed term = intercept of the zero-inflation  
#remaining terms are trait (Poi & Zi) specific contrasts from the intercept.

## trace plots ##
plot(male_LBS_model)


## model output ##
summary(male_LBS_model)


#############################################################################################################
## saving output #############################################################################################
###########################################################################################################
save(male_LBS_model, file = "PhD_4th_yr/Inbreeding_depression_models/LBS/OUTPUT_Male_LBS_mm_1.5mit.RData")

load(file = "PhD_4th_yr/Inbreeding_depression_models/LBS/OUTPUT_Male_LBS_mm_1.5mit.RData")


sols_full<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(matches("FROH"))%>% ## taking out sols with FROH included
  rename(LBS_poiss="traitLBS:FROHsum")%>%
  dplyr::mutate(across(3:35, ~.x + LBS_poiss)) ## adding FROHsum to chrFROH values

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


FROH_sols$CHR<-as.factor(FROH_sols$CHR)

ggplot(data=FROH_sols, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange(colour="grey50") + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  geom_hline(yintercept=-0.75470641, lty=1,colour="red") + ## red line is the average effect of all chromosomes 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="solution + CI", title = "Effect of chromosome FROH on male LBS poisson") +
  theme_bw()+  # use a white background                 
  theme(legend.position = "none")




### Overall effect of FROH (across all chr) on poisson process of LBS ###
## Same output as the simple LBS model .. and highly significant in that one


LBS_froh_all=sum(FROH_sols$solution) 

## extracting mean estimates from values

traitLBS=mean(male_LBS_model$Sol[,1])


#assuming FROH of 0.25/0.125
1-plogis((0.125*LBS_froh_all)+(traitLBS))
1-plogis((0.25*LBS_froh_all)+(traitLBS))
1-plogis((0.25*LBS_froh_all))
1-plogis((0.125*LBS_froh_all))

plogis(14.2)

plogis(male_LBS_model$Sol[, 2]/sqrt(1 + c2))



plogis((0.125*LBS_froh_all)+(traitLBS))


exp(plogis((0.125*LBS_froh_all)+(traitLBS)))

exp(plogis((0.125*LBS_froh_all)+(traitLBS)))/(1+exp(plogis((0.125*LBS_froh_all)+(traitLBS))))



traitzi_LBS=mean(male_LBS_model$Sol[,2])
traitLBS_FROHsum=mean(male_LBS_model$Sol[,3])
traitzi_LBS_FROHsum=mean(male_LBS_model$Sol[,4])



