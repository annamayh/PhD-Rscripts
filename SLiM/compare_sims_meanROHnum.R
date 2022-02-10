library("purrr")
library(dplyr)
library(plyr)
library(ggplot2)
library(forcats)
library(wesanderson)
library(patchwork)
library(bobfunctions2)

setwd("H:/")




######################################################################################################################################
######## load all simulation files in #####

## neutral model
neutral_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/Rum/Neutral_re_an_noMAF.2/*/*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(neutral_list)){
  neutral_list[[i]]$simulation_num<-i #adding sim number for grouping later
} 

#selection and recombination
sel_reomb_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/Rum/NewH_sel00_r/sel_r_*/ROH_out*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(sel_reomb_list)){
  sel_reomb_list[[i]]$simulation_num<-i #adding sim number for grouping later
} 

## stronger selection coefficient 
sel_s0.05_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/Rum/NewH_strongsel00_r/*/ROH_out*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(sel_s0.05_list)){
  sel_s0.05_list[[i]]$simulation_num<-i #adding sim number for grouping later
}

## stronger selection coefficient 
recomb <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/Rum/Recomb_noMAF/*/ROH_out*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(sel_s0.05_list)){
  sel_s0.05_list[[i]]$simulation_num<-i #adding sim number for grouping later
}


#############################
## servere bottleneck neutral 
##############################
bottleneck <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/btl/Neutral_re_an_noMAF.2/*/*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck)){
  bottleneck[[i]]$simulation_num<-i #adding sim number for grouping later
}

##bottlneck s+r
bottleneck_s_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/btl/NewH_sel00_r/sel_r_*/ROH_out*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_s_r)){
  bottleneck_s_r[[i]]$simulation_num<-i #adding sim number for grouping later
}

## bottleneck s = 0.05
bottleneck_s0.05_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/btl/NewH_strongsel00_r/*/ROH_out*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_s0.05_r)){
  bottleneck_s0.05_r[[i]]$simulation_num<-i #adding sim number for grouping later
}

bottleneck_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/btl/recomb/*/ROH_out*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_s0.05_r)){
  bottleneck_s0.05_r[[i]]$simulation_num<-i #adding sim number for grouping later
}

##########################
## no bottleneck neutral##
##########################
Nobottleneck <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/NObtl/Neutral_noMAF/*/ROH_out*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(Nobottleneck)){
  Nobottleneck[[i]]$simulation_num<-i #adding sim number for grouping later
}

##NObottlneck s+r
NObottleneck_s_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/NObtl/NewH_sel00_r*/*/ROH_out*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(NObottleneck_s_r)){
  NObottleneck_s_r[[i]]$simulation_num<-i #adding sim number for grouping later
}

## NObottleneck s = 0.05

NObottleneck_s0.05_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/NObtl/NewH_strongsel00_r/*/ROH_out*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(NObottleneck_s0.05_r)){
  NObottleneck_s0.05_r[[i]]$simulation_num<-i #adding sim number for grouping later
}

NObottleneck_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Dec_21_current_out/NObtl/Recomb_noMAF/*/ROH_out*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(NObottleneck_s0.05_r)){
  NObottleneck_s0.05_r[[i]]$simulation_num<-i #adding sim number for grouping later
}
###########################################################################################

## overall mean

unlist_f<-function(listed_sims){
  bind_rows(listed_sims, .id = "simulation_num")%>%
    select(simulation_num, NSEG,KB)

}

unRum_neutral_ROH<-neutral_list%>%unlist_f()
unRum_recomb<-recomb%>%unlist_f()
unRum_sel<-sel_reomb_list%>%unlist_f()
unRum_strongsel<-sel_s0.05_list%>%unlist_f()

unbtl_neutral<-bottleneck%>%unlist_f()
unbtl_recomb<-bottleneck_r%>%unlist_f()
unbtl_sel<-bottleneck_s_r%>%unlist_f()
unbtl_strongsel<-bottleneck_s0.05_r%>%unlist_f()

unNObtl_neutral<-Nobottleneck%>%unlist_f()
unNObtl_recomb<-NObottleneck_r%>%unlist_f()
unNObtl_sel<-NObottleneck_s_r%>%unlist_f()
unNObtl_strongsel<-NObottleneck_s0.05_r%>%unlist_f()


mean(unRum_neutral_ROH$KB)
mean(unRum_recomb$KB)
mean(unRum_sel$KB)
mean(unRum_strongsel$KB)

mean(unbtl_neutral$KB)
mean(unbtl_recomb$KB)
mean(unbtl_sel$KB)
mean(unbtl_strongsel$KB)
  
mean(unNObtl_neutral$KB)
mean(unNObtl_recomb$KB)
mean(unNObtl_sel$KB)
mean(unNObtl_strongsel$KB)
