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



iHS_list<-map(scan_output, ihh2ihs)


## read in locations of hotspots for each simulation 

ROH_hotspots_list<-list()

for(l in 1:27) {
  if (l==20) next ## for failed simulations on dud machine
  if (l==22) next
  if (l==25) next
  if (l==26) next

  ROH_hotspots_list[[l]]<-read.table(paste0("Slim/Ash_server_output/neutral/hom summary files/ROH_out",l,".hom.summary"),header = TRUE, stringsAsFactors=FALSE)%>%
    select(BP,UNAFF)%>%mutate(pall=UNAFF/100)%>%setNames(.,c("POSITION","UNAFF","pall"))
  num[l]<-as.numeric(quantile(ROH_hotspots_list[[l]]$pall, c(.99)))
  ROH_hotspots_list[[l]]$hotspot<-NA
  ROH_hotspots_list[[l]]$hotspot<-ifelse(ROH_hotspots_list[[l]]$pall>num[l], "yes", "no")
} 
  
names(ROH_hotspots_list)<-paste0("Sim_num",1:27) #naming the list elements based on the simulation number 

ROH_hotspots_list["Sim_num20"]<-NULL
ROH_hotspots_list["Sim_num22"]<-NULL #remove the simulations that didnt work on dud machine
ROH_hotspots_list["Sim_num25"]<-NULL
ROH_hotspots_list["Sim_num26"]<-NULL 



show_hotspots<-function(sim){
  num<-as.numeric(quantile(sim$pall, c(.99)))
  sim$hotspot<-NA
  sim$hotspot<-ifelse(sim$pall>num, "yes", "no")
  
}


ROH_hotspots_list_y<-map(ROH_hotspots_list, show_hotspots)

testing<-as.data.frame(ROH_hotspots_list["Sim_num1"])

for (line in 1:nrow(testing["Sim_num1"]))
num<-as.numeric(quantile(testing[["Sim_num1"]]$pall, c(.99)))
testing$hotspot<-NA
testing$hotspot<-ifelse(sim$pall>num, "yes", "no")
