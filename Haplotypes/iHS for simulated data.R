library(vcfR)
library(rehh)
library(data.table)
library(R.utils)
library(dplyr)
library(plyr)
library(purrr)
library(ggplot2)
library(tidyverse)
library(janitor)
library(stringdist)

setwd("H:/PHD_2ndYR")

scan_output<-list()

for(k in 1:27) {
  if (k==20) next ## for failed simulations on dud machine
  if (k==22) next
  if (k==25) next
  if (k==26) next
  # haplotype file name for each chromosome
  hap_file_input = paste0("Slim/Ash_server_output/neutral/SLiM_raw_output/neutral_sim",k)
  # create internal representation
  hh <- data2haplohh(hap_file = hap_file_input,
                     polarize_vcf = FALSE,
                     vcf_reader = "data.table", 
                     remove_multiple_markers = TRUE
                     )
  # perform scan on a single chromosome (calculate iHH values)
  scan_output[[k]] <- scan_hh(hh)
  scan_output[[k]]$sim<-k
  
}  

#save(scan_output, file="iHS_simulations_scan_output.Rdata")
#load("Haplotype_diversity/iHS_simulations_scan_output.Rdata")

names(scan_output)<-paste0("Sim_num",1:27) #naming the list elements based on the simulation number 

scan_output["Sim_num20"]<-NULL
scan_output["Sim_num22"]<-NULL #remove the simulations that didnt work on dud machine
scan_output["Sim_num25"]<-NULL
scan_output["Sim_num26"]<-NULL 



iHS_list<-map(scan_output, ihh2ihs)# getting iHS values 


## read in locations of hotspots for each simulation 

ROH_hotspots_list<-list()

for(l in 1:27) {
  if (l==20) next ## for failed simulations on dud machine
  if (l==22) next
  if (l==25) next
  if (l==26) next

  ROH_hotspots_list[[l]]<-read.table(paste0("Slim/Ash_server_output/neutral/hom summary files/ROH_out",l,".hom.summary"),header = TRUE, stringsAsFactors=FALSE)%>%
    select(BP,UNAFF)%>%mutate(pall=UNAFF/100)%>%setNames(.,c("POSITION","UNAFF","pall"))
  
  num[l]<-as.numeric(quantile(ROH_hotspots_list[[l]]$pall, c(.99)))# getting hotspot threshold for each sim
  ROH_hotspots_list[[l]]$hotspot<-NA #now accessing the list and adding column to populate
  ROH_hotspots_list[[l]]$hotspot<-ifelse(ROH_hotspots_list[[l]]$pall>num[l], "yes", "no")#filling column based on each number 
} 
  
names(ROH_hotspots_list)<-paste0("Sim_num",1:27) #naming the list elements based on the simulation number 

ROH_hotspots_list["Sim_num20"]<-NULL
ROH_hotspots_list["Sim_num22"]<-NULL #remove the simulations that didnt work on dud machine
ROH_hotspots_list["Sim_num25"]<-NULL
ROH_hotspots_list["Sim_num26"]<-NULL #can maybe gett rid of this as repeated later'


# now to join the two lists together to compare iHS and hotspots 

hotspots_vs_iHS<-list()

for(s in 1:27) {
  if (s==20) next ## for failed simulations on dud machine
  if (s==22) next
  if (s==25) next
  if (s==26) next

  hotspots_vs_iHS[[s]]<-join(iHS_list[[paste0("Sim_num",s)]]$ihs, ROH_hotspots_list[[paste0("Sim_num",s)]])


}

names(hotspots_vs_iHS)<-paste0("Sim_num",1:27)
hotspots_vs_iHS["Sim_num20"]<-NULL
hotspots_vs_iHS["Sim_num22"]<-NULL #remove the simulations that didnt work on dud machine
hotspots_vs_iHS["Sim_num25"]<-NULL
hotspots_vs_iHS["Sim_num26"]<-NULL 



unlisted_hotspots_iHS<-bind_rows(hotspots_vs_iHS, .id = "Sim_num") # unlist for plotting
#.id is the name of the column which fills in whatever the lists are called 


ggplot(unlisted_hotspots_iHS, aes(hotspot, IHS, fill=hotspot)) + geom_boxplot()+
  theme_classic()+
  labs(title = "Haplotype diveristy using iHS")


ggplot(unlisted_hotspots_iHS, aes(hotspot, LOGPVALUE, fill=hotspot)) + geom_boxplot()+
  theme_classic()+
  labs(title = "iHS p-values in ROH hotspots Vs not in neutral simulations (23 iterations")


#########################################################################################################
################## now for the sims with selection and recombination ####################################
#########################################################################################################


scan_output<-list()

for(k in 1:28) {
  if (k==23) next ## for failed simulations on dud machine
  if (k==26) next
 
  # haplotype file name for each chromosome
  hap_file_input = paste0("Slim/Ash_server_output/Sel_d0.1/SLiM_raw output/sel_d0.1_sim",k)
  # create internal representation
  hh <- data2haplohh(hap_file = hap_file_input,
                     polarize_vcf = FALSE,
                     vcf_reader = "data.table", 
                     remove_multiple_markers = TRUE
  )
  # perform scan on a single chromosome (calculate iHH values)
  scan_output[[k]] <- scan_hh(hh)
  scan_output[[k]]$sim<-k
  
}  

#save(scan_output, file="iHS_simulations_selection_scan_output.Rdata")
load("iHS_simulations_selection_scan_output.Rdata")

names(scan_output)<-paste0("Sim_num",1:28) #naming the list elements based on the simulation number 

scan_output["Sim_num1"]<-NULL
scan_output["Sim_num5"]<-NULL
scan_output["Sim_num23"]<-NULL
scan_output["Sim_num26"]<-NULL #remove the simulations that didnt work on dud machine



iHS_list<-map(scan_output, ihh2ihs)# getting iHS values 


## read in locations of hotspots for each simulation 

ROH_hotspots_list<-list()

test<-ROH_hotspots_list[[2]]

for(l in 1:28) {
  if (l==1) next
  if (l==5) next
  if (l==23) next ## for failed simulations on dud machine
  if (l==26) next
  
  ROH_hotspots_list[[l]]<-read.table(paste0("Slim/Ash_server_output/Sel_d0.1/hom_summary_files/sel_d0.1ROH_out",l,".hom.summary"),header = TRUE, stringsAsFactors=FALSE)%>%
    select(BP,UNAFF)%>%mutate(pall=UNAFF/100)%>%setNames(.,c("POSITION","UNAFF","pall"))
  
  num<-as.numeric(quantile(ROH_hotspots_list[[l]]$pall, c(.99),))# getting hotspot threshold for each sim
  ROH_hotspots_list[[l]]$hotspot<-NA #now accessing the list and adding column to populate
  ROH_hotspots_list[[l]]$hotspot<-ifelse(ROH_hotspots_list[[l]]$pall>num, "yes", "no")#filling column based on each number 
} 

names(ROH_hotspots_list)<-paste0("Sim_num",1:28) #naming the list elements based on the simulation number 

ROH_hotspots_list["Sim_num1"]<-NULL
ROH_hotspots_list["Sim_num5"]<-NULL #remove the simulations that didnt work on dud machine
ROH_hotspots_list["Sim_num23"]<-NULL
ROH_hotspots_list["Sim_num26"]<-NULL #can maybe gett rid of this as repeated later'


# now to join the two lists together to compare iHS and hotspots 

hotspots_vs_iHS<-list()

for(s in 1:28) {
  if (s==1) next ## for failed simulations on dud machine
  if (s==5) next
  if (s==23) next
  if (s==26) next
  
  hotspots_vs_iHS[[s]]<-join(iHS_list[[paste0("Sim_num",s)]]$ihs, ROH_hotspots_list[[paste0("Sim_num",s)]])
  
  
}

names(hotspots_vs_iHS)<-paste0("Sim_num",1:28)
hotspots_vs_iHS["Sim_num1"]<-NULL
hotspots_vs_iHS["Sim_num5"]<-NULL #remove the simulations that didnt work on dud machine
hotspots_vs_iHS["Sim_num23"]<-NULL
hotspots_vs_iHS["Sim_num26"]<-NULL 



unlisted_hotspots_iHS<-bind_rows(hotspots_vs_iHS, .id = "Sim_num") # unlist for plotting
#.id is the name of the column which fills in whatever the lists are called 


ggplot(unlisted_hotspots_iHS, aes(hotspot, IHS, fill=hotspot)) + geom_boxplot()+
  theme_classic()+
  labs(title = "iHS in simulations with selection and recombination")


ggplot(unlisted_hotspots_iHS, aes(hotspot, LOGPVALUE, fill=hotspot)) + geom_boxplot()+
  theme_classic()+
  labs(title = "iHS p-values in ROH hotspots Vs not hotspots \nin simulations with selection and recombination (23 iterations)")


