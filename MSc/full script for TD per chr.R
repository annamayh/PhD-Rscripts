library(GenABEL)
library(plyr)
library(reshape2)
library(dplyr)
library(magrittr)
##
setwd("M:/project")

load("20180209_SoayPlates1-83_GenABEL.RData")

basedata <- read.table("20190208_SoayBaseData_edit.txt", header = T, stringsAsFactors = F)
head(basedata)
sexdata <- select(basedata, ID, SEX)


for(h in 21:27){
  ## for 1:30 chr (ccmputer couldne come with all at the same time)
  
  soayabel6 <- soayabel[,chromosome(soayabel) == h] #
  
  
  
  nsnps(soayabel6)
  
  pedi <- read.table("20190208_Full_Pedigree.txt", header = T)  ################################ header = T
  
  
  restab <- pedi                 ##################### Names are ID, MOTHER, FATHER
  head(restab) #creating master table 
  restab <- join(restab, sexdata)
  
  genotab <- as.double.gwaa.data(soayabel6)#table with id snp and id geno
  genotab <-melt(genotab)#combining 
  
  names(genotab) <- c("ID", "SNP", "ID_geno")
  summary(genotab)
  
  rtab <- join(restab,genotab) ###########doesnt work for whole genome ==need to subset by chromosome and do a loop :S
  
  head(rtab)
  tail(rtab)
  
  #adding mothers geno
  
  names(genotab) <- c("MOTHER", "SNP", "MUM_geno")
  
  rtab <- join(rtab,genotab)
  
  
  #adding fathers geno
  
  names(genotab) <- c("FATHER", "SNP", "DAD_geno")
  
  rtab <- join(rtab,genotab)
  
  head(rtab)
  
  
  #Determining allele tranmission 
  
  rtab$MUM_allele <- NA #creating columns to populate 
  rtab$DAD_allele <- NA
  
  #MUM allele
  
  rtab$MUM_allele <- ifelse(rtab$MUM_geno == 1 & rtab$ID_geno == 1 & rtab$DAD_geno == 2, "A",   #################### CHECK THESE
                            ifelse(rtab$MUM_geno == 1 & rtab$ID_geno == 1 & rtab$DAD_geno == 0, "B",
                                   ifelse(rtab$ID_geno == 0 & rtab$MUM_geno == 1, "A",
                                          ifelse(rtab$ID_geno == 2 & rtab$MUM_geno == 1, "B", NA))))          ##### DONT PUT NA in brackets
  
  #DAD ALLELE
  rtab$DAD_allele <- ifelse(rtab$DAD_geno == 1 & rtab$ID_geno == 1 & rtab$MUM_geno == 2,"A", 
                            ifelse(rtab$DAD_geno == 1 & rtab$ID_geno == 1 & rtab$MUM_geno == 0, "B",
                                   ifelse(rtab$ID_geno == 0 & rtab$DAD_geno == 1, "A",
                                          ifelse(rtab$ID_geno == 2 & rtab$DAD_geno == 1, "B", NA))))
  
  head(rtab)
  
  
  rtab <- subset(rtab, !is.na(SNP))
  
  
  otab <- select(rtab, -MUM_geno, -DAD_geno) #otab is same as rtab but remove mum and dad geno
  
  otab <- melt(otab, id.vars = c("ID", "MOTHER", "FATHER", "SEX", "SNP", "ID_geno"))
  
  
  
  
  
  
  #### Want to summarise which allele was transmitted by parent SEX
  
  #mum transmissions 
  mumtrans <- rtab %>%
    group_by(SNP, MUM_allele) %>%
    summarize(n()) %>%
    na.omit
  
  mumtrans
  
  #dad transmission 
  dadtrans <- rtab %>%
    group_by(SNP, DAD_allele) %>%
    summarize(n())  %>%  #number of individuals for each snp 
    na.omit
  
  #offspring transmission 
  offspringtrans <- otab %>%
    group_by(SNP, SEX, value) %>%
    summarize(n()) %>%
    na.omit
  
  ###~~ MUM TD ~~~#####
  
  mumtrans <- dcast(mumtrans, SNP ~ MUM_allele) #dcast = oppposite of melt 
  mumtrans[is.na(mumtrans)] <- 0 
  mumtrans$p.value <- NA
  
  for(i in 1:nrow(mumtrans)){
    
    
    x1 <- mumtrans$A[i] #observed
    x2 <- mumtrans$B[i]
    x3 <- (x1 + x2)/2
    x4 <- x3 #expected
    
    x <- matrix(c(x1, x2, x3, x4), nrow = 2)# creating table with expected vs observed to test
    suppressWarnings(fit1 <- fisher.test(x))
    
    mumtrans$p.value[i] <-  fit1$p.value
    
  }
  
  2.245e-6
  
  mumtrans <- arrange(mumtrans, p.value)
  mumtrans
  
  
  
  #####~~~DAD TD~~~~####
  
  dadtrans <- dcast(dadtrans, SNP ~ DAD_allele) #dcast = oppposite of melt 
  dadtrans[is.na(dadtrans)] <- 0 #new line changin NAs to 0 
  dadtrans$p.value <- NA
  
  for(i in 1:nrow(dadtrans)){
    
    
    x1 <- dadtrans$A[i]
    x2 <- dadtrans$B[i]
    x3 <- (x1 + x2)/2
    x4 <- x3
    
    x <- matrix(c(x1, x2, x3, x4), nrow = 2)  
    suppressWarnings(fit1 <- fisher.test(x))
    
    dadtrans$p.value[i] <-  fit1$p.value
    
  }
  
  
  dadtrans <- arrange(dadtrans, p.value)
  dadtrans
  
  
  
  ####~~~OFFSPRING TD~~~~###
  
  names(offspringtrans)[names(offspringtrans) == "value"] <- "allele"
  offspringtrans <- dcast(offspringtrans, SNP + SEX ~ allele) #dcast = oppposite of melt 
  offspringtrans[is.na(offspringtrans)] <- 0 
  offspringtrans$p.value <- NA
  
  for(i in 1:nrow(offspringtrans)){
    
    
    x1 <- offspringtrans$A[i]
    x2 <- offspringtrans$B[i]
    x3 <- (x1 + x2)/2
    x4 <- x3
    
    x <- matrix(c(x1, x2, x3, x4), nrow = 2)  
    suppressWarnings(fit1 <- fisher.test(x))
    
    offspringtrans$p.value[i] <-  fit1$p.value
    
  }
  
  
  offspringtrans <- arrange(offspringtrans, p.value)
  offspringtrans
  
  save(mumtrans, dadtrans, offspringtrans, file = paste0("TD.Chromosome.", h, ".RData")) #saving results as r data file
  
}

## ONLY NNEEDS TO RUN ONCE ###
fullmumtrans <- NULL #creating null tables to populate
fulldadtrans <- NULL 
fulloffspringtrans <- NULL
############################################################################

#### creating large dataframe with all results from all chr ####
for(h in 21:27){
  
  load(paste0("TD.Chromosome.", h, ".RData")) #loading results
  
  fullmumtrans <- rbind(fullmumtrans, mumtrans)#rbinding results to one large table 
  
  
}


for(h in 21:27){
  
  load(paste0("TD.Chromosome.", h, ".RData"))
  
  fulldadtrans <- rbind(fulldadtrans, dadtrans)
  
  
}

for(h in 21:27){
  
  load(paste0("TD.Chromosome.", h, ".RData"))
  
  fulloffspringtrans <- rbind(fulloffspringtrans, offspringtrans)
  
  
}



save(fullmumtrans, fulldadtrans, fulloffspringtrans, file = "data.Rdata") #saving all tables as R data



write.table(fullmumtrans,
            file = "Full_Mother_TD_Results.txt",
            row.names = F, quote = F, sep = "\t") #saving tables as txt file 
write.table(fulldadtrans,
            file = "Full_Father_TD_Results.txt",
            row.names = F, quote = F, sep = "\t")
write.table(fulloffspringtrans,
            file = "Full_Offspring_TD_Results.txt",
            row.names = F, quote = F, sep = "\t")






