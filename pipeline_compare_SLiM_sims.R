##########################################
### using plink to search for ROH in R ###
##########################################

library("purrr")
library(dplyr)
library(ggplot2)

setwd("H:/") # make sure directory is the same location as the plink application #
## Plink_ROH function is a function to use plink to search for ROH with different input files(vcf/bfile)and output file names 
## file path can be changed directly from the plink input file path or the input file name later

plink_ROH<-function(Input_file, Output_file){
  
  system(paste0("plink --vcf PHD_2ndYR/Slim/Output/",Input_file," ",
                "--autosome --autosome-num 33 --out PHD_2ndYR/Slim/Output/ROH_output/",Output_file," ",
                "--homozyg --homozyg-window-snp 35 --homozyg-snp 40 --homozyg-kb 2500 ",
                "--homozyg-density 70 --homozyg-window-missing 4 ",
                "--homozyg-het 0 ",
                "--maf 0.01 --freq --missing"))
  }


####### NEUTRAL MODEL ##################################################################################

## loop function over multiple files ##

for (i in 1:20){
  
  vcf<-paste0("Neutral_model_outputs/neutral_m2_sim",i) ## all files in model output are called test_neutral_loop1 etc 
  output<-paste0("neutral_m2_sim",i,"_ROH_OUT") ## want output files to be numbered too 
  
  plink_ROH(vcf,output) #run function for all files labeled 1:10
  
  
}
##############################################################################################################
                  #### once this has run once dont need to repeat it (files will be saved) ####
###########################################################################################################

## sys.glob read all files in the form "testing_output_loop*.hom.summary" and puts them into a list 
#super handy 

unaff_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Output/ROH_output/neutral_m2_sim*_ROH_OUT.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)

# function to work out pall (prop of ids w/ ROH at SNP)

for (i in 1:20){
  unaff_list[[i]]$simulation_num<-i #adding sim number for grouping later
  unaff_list[[i]]$pall<-NA
  for(k in 1:nrow(unaff_list[[i]])){
    unaff_list[[i]]$pall[k] <- unaff_list[[i]]$UNAFF[k]/100#pall per SNP
  }
  
}
  


unlisted<-bind_rows(unaff_list, .id = "simulation_num")#unlisting


pall_per_sim<-unlisted%>%select(simulation_num, pall)%>%
  group_by(simulation_num)%>% #grouping by simulation iteration
  dplyr::summarise(Mean_pall=mean(pall))#getting mean pall of simulation

mean(pall_per_sim$Mean_pall)

## can try to get more summary stats e.g. top 1% etc.


top_1_func<-function (x){quantile(x$pall, c(.99))}

top_1_df<-list()
top_1_df<-map(unaff_list,top_1_func)%>%bind_rows(.id = "Value")
names(top_1_df)<-(c("sim","Top_1"))

mean(top_1_df$Top_1)

#### Then can load more simulations e.g. neutral vs selection vs recomb to compare ####


####################################################################################################
################### SELECTION #####################################################################
##################################################################################################

for (i in 1:20){
  
  vcf<-paste0("Selection_model_outputs/Selection_m2_sim",i) ## all files in model output are called test_neutral_loop1 etc 
  output<-paste0("Selection_m2_sim",i,"_ROH_OUT") ## want output files to be numbered too 
  
  plink_ROH(vcf,output) #run function for all files labeled 1:10
  
  
}

unaff_list_sel <- lapply(Sys.glob("PHD_2ndYR/Slim/Output/ROH_output/Selection_m2_sim*_ROH_OUT.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)

# function to work out pall (prop of ids w/ ROH at SNP)

for (i in 1:20){
  unaff_list_sel[[i]]$simulation_num<-i #adding sim number for grouping later
  unaff_list_sel[[i]]$pall<-NA
  for(k in 1:nrow(unaff_list_sel[[i]])){
    unaff_list_sel[[i]]$pall[k] <- unaff_list_sel[[i]]$UNAFF[k]/100#pall per SNP
  }
  
}

unlisted_sel<-bind_rows(unaff_list_sel, .id = "simulation_num")#unlisting


pall_per_sim_sel<-unlisted_sel%>%select(simulation_num, pall)%>%
  group_by(simulation_num)%>% #grouping by simulation iteration
  dplyr::summarise(Mean_pall=mean(pall))#getting mean pall of simulation

mean(pall_per_sim_sel$Mean_pall)



top_1_func<-function (x){quantile(x$pall, c(.99))}

top_1_df_sel<-list()
top_1_df_sel<-map(unaff_list_sel,top_1_func_sel)%>%bind_rows(.id = "Value")
names(top_1_df_sel)<-(c("sim","Top_1"))

mean(top_1_df_sel$Top_1)



################# PLOTS ####################################################################################

ggplot(pall_per_sim, aes(x="Mean_pall", y=Mean_pall))+geom_boxplot()+geom_dotplot(binaxis='y', stackdir='center',
                                                                                  position=position_dodge(1))


ggplot(top_1, aes(x="Top_1", y=Top_1))+geom_boxplot()#+geom_dotplot()

