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

output_haplohh<-list()

for(k in 1:27) {
  if (k==20) next ## for failed simulations on dud machine
  if (k==22) next
  if (k==25) next
  if (k==26) next
  # haplotype file name for each chromosome
  hap_file_input = paste0("Slim/Ash_server_output/neutral/SLiM_raw_output/neutral_sim",k)
  # create internal representation
  output_haplohh[[k]] <- data2haplohh(hap_file = hap_file_input,
                     polarize_vcf = FALSE,
                     vcf_reader = "data.table", 
                     remove_multiple_markers = TRUE
                     )
  # perform scan on a single chromosome (calculate iHH values)
  scan <- scan_hh(hh)
  scan$sim<-k
  
  ihs<-ihh2ihs(scan)
  
  if (k == 1) {
    ihs_all_sims <- ihs
  } else {
    ihs_all_sims <- rbind(ihs_all_sims, ihs)
  }
}




iHS_all_sims<-as.data.frame(ihs_all_sims)

ihs_SNPs<-tibble::rownames_to_column(iHS_genomewide, "SNP")

