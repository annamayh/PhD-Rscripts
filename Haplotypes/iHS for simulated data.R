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

names(scan_output)<-paste0("Sim_num",1:21) #naming the list elements based on the simulation number 

scan_output["Sim_num20"]<-NULL
scan_output["Sim_num22"]<-NULL #remove the simulations that didnt work on dud machine
scan_output["Sim_num25"]<-NULL
scan_output["Sim_num26"]<-NULL ## this sometimes doesnt work??



iHS_list<-map(scan_output, ihh2ihs)


