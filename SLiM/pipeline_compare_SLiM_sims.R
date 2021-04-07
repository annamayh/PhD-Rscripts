##########################################
### using plink to search for ROH in R ###
##########################################

library("purrr")
library(dplyr)
library(ggplot2)

setwd("H:/") # make sure directory is the same location as the plink application #

###FUNCTIONS####

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


## function to add pall values to a list containining all iterations 

add_pall<-function(iterations_list){
  
  for (i in 1:20){
    iterations_list[[i]]$simulation_num<-i #adding sim number for grouping later
    iterations_list[[i]]$pall<-NA
    for(k in 1:nrow(iterations_list[[i]])){
      iterations_list[[i]]$pall[k] <- iterations_list[[i]]$UNAFF[k]/100}#pall per SNP
  }
    
  
}


## function to unlist a list containing all iterations 
unlist<-function(iteration_list){bind_rows(iteration_list, .id = "simulation_num")}

## function to get the mean pall per simulation for the unlisted dataframe
get_mean_pall_per_sim<-function(unlisted_df){
  
  unlisted_df%>%select(simulation_num, pall)%>%
    group_by(simulation_num)%>% #grouping by simulation iteration
    dplyr::summarise(Mean_pall=mean(pall))#getting mean pall of simulation
  
  
}

#func to get top 1% in a datafram
top_1_func<-function (x){quantile(x$pall, c(.99))}

#function using function above to apply to all dataframes in a list of iterations
top_1_make_df<-function(iteration_list){
  
  map(iteration_list,top_1_func)%>%bind_rows(.id = "Value")%>%setNames(.,c("sim","Top_1"))
  
}

####### NEUTRAL MODEL ##################################################################################

## loop function over multiple files ##

for (i in 1:15){
  
  vcf<-paste0("Selection_model_outputs/selection_ben_only_m1_sim",i) ## all files in model output are called test_neutral_loop1 etc 
  output<-paste0("selection_ben_only_m1_sim",i,"_ROH_OUT") ## want output files to be numbered too 
  
  plink_ROH(vcf,output) #run function for all files labeled 1:10
  
  
}
##############################################################################################################
                  #### once this has run once dont need to repeat it (files will be saved) ####
###########################################################################################################

## sys.glob read all files in the form "testing_output_loop*.hom.summary" and puts them into a list 
#super handy 

unaff_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Output/ROH_output/selection_ben_only_m1_sim*_ROH_OUT.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)


for (i in 1:15){
  unaff_list[[i]]$simulation_num<-i #adding sim number for grouping later
  unaff_list[[i]]$pall<-NA
  for(k in 1:nrow(unaff_list[[i]])){
    unaff_list[[i]]$pall[k] <- unaff_list[[i]]$UNAFF[k]/100#pall per SNP
  }
  
}

### use functions above ###

mean_pall_df_neutral<-unaff_list_neutral%>%unlist()%>%get_mean_pall_per_sim()
mean(mean_pall_df_ben$Mean_pall)

###

top_1_df_neu<-top_1_make_df(unaff_list_neutral)
mean(top_1_df_neu$Top_1)








#########################################################################################################
#########################################################################################################

mean_pall_df_ben$model<-"ben_selection"
mean_pall_df_del$model<-"del_selection"
mean_pall_df_neutral$model<-"neutral"

plot_table<-rbind(mean_pall_df_ben,mean_pall_df_del,mean_pall_df_neutral)

###
top_1_df_neu$model<-"Neutral"
top_1_df_del$model<-"del_selection"
top_1_df_ben$model<-"ben_selection"

plot_table_1per<-rbind(top_1_df_neu,top_1_df_del,top_1_df_ben)



################# PLOTS ####################################################################################

library(viridis)

ggplot(plot_table, aes(x=model, y=Mean_pall, fill=model))+
  geom_boxplot()+
  geom_jitter(color="black", size=1.5, alpha=0.9) 


ggplot(plot_table_1per, aes(x=model, y=Top_1, fill=model))+
  geom_boxplot()+
  geom_jitter(color="black", size=1.5, alpha=0.9) 


