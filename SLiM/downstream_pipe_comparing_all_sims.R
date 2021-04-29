
library("purrr")
library(dplyr)
library(plyr)
library(ggplot2)

setwd("H:/")

## function to add pall values to a list containining all iterations 
add_pall<-function(iteration_list){
  iteration_list%>%mutate(pall=UNAFF/100)
}

add_pall_2000<-function(iteration_list){
  iteration_list%>%mutate(pall=UNAFF/2000)
}


## 
add_simulation_num<-function(list){
  for (i in 1:length(list)){
    list[[i]]$simulation_num<-i} 
  
}

## function to unlist a list containing all iterations 
unlist<-function(iteration_list){bind_rows(iteration_list, .id = "simulation_num")}

## function to get the mean pall per simulation for the unlisted dataframe
unlist_then_get_mean_pall_per_sim<-function(listed_sims){
  bind_rows(listed_sims, .id = "simulation_num")%>%
  select(simulation_num, pall)%>%
    group_by(simulation_num)%>% #grouping by simulation iteration
    dplyr::summarise(Mean_pall=mean(pall))#getting mean pall of simulation
 }



######################################################################################################################################
######## load all simulation files in #####

## neutral model
neutral_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Output/ROH_output/neutral_m6_sim*_ROH_OUT.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(neutral_list)){
  neutral_list[[i]]$simulation_num<-i #adding sim number for grouping later
} 

#recombination model
recomb_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Sel_recomb/sel_recomb_*/sel_recomb_ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(recomb_list)){
  recomb_list[[i]]$simulation_num<-i #adding sim number for grouping later
} 

### selection only
sel_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Sel_d0.1/sel_d0.1_*/sel_d0.1ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(sel_list)){
  sel_list[[i]]$simulation_num<-i #adding sim number for grouping later
} 

#selection and recombination
sel_reomb_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Sel_recomb_d0.1/sel_recomb_d0.1_*/sel_recomb_d0.1ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(sel_reomb_list)){
  sel_reomb_list[[i]]$simulation_num<-i #adding sim number for grouping later
} 

## stronger selection coefficient 
sel_s0.05_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Sel_s0.05d0.1/sel_s0.05d0.1_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(sel_s0.05_list)){
  sel_s0.05_list[[i]]$simulation_num<-i #adding sim number for grouping later
}

## servere bottleneck 
bottleneck <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Neutral_bottleneck/Neutral_btl_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck)){
  bottleneck[[i]]$simulation_num<-i #adding sim number for grouping later
}

## no bottleneck
Nobottleneck <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Neutral_NObottleneck/Neutral_nobtl_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(Nobottleneck)){
  Nobottleneck[[i]]$simulation_num<-i #adding sim number for grouping later
}
###########################################################################################

pall_neutral<-neutral_list%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_recomb<-recomb_list%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_sel<-sel_list%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_sel_recomb<-sel_reomb_list%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_sel_s0.05<-sel_s0.05_list%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_bottleneck<-bottleneck%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_nobtl<-Nobottleneck%>%map(add_pall_2000)%>%unlist_then_get_mean_pall_per_sim


pall_neutral$Model<-"Neutral"
pall_recomb$Model<-"Varied Recomb"
pall_sel$Model<-"Selection"
pall_sel_recomb$Model<-"Varied recomb + Selection"
pall_sel_s0.05$Model<-"Higher selection coefficient"
pall_nobtl$Model<-"Neutral No bottle"
pall_bottleneck$Model<-"Neutral servere bottle"


plot_pall_mean<-rbind(pall_neutral,pall_recomb,pall_sel,pall_sel_recomb,pall_sel_s0.05,pall_nobtl,pall_bottleneck)

ggplot(plot_pall_mean, aes(x=Model, y=Mean_pall, fill=Model))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_brewer(palette = "Pastel2")+
  geom_jitter(color="black", size=1.5, alpha=0.9) 



###################################################################################################################################


#func to get top 1% in a datafram
top_1_func<-function (x){quantile(x$pall, c(.99))}

#function using function above to apply to all dataframes in a list of iterations
top_1_make_df<-function(iteration_list){
  
  map(iteration_list,top_1_func)%>%bind_rows(.id = "Value")%>%setNames(.,c("sim","Top_1"))
  
}

per_sim_neutral$model<-"Neutral"
top_1_neutral<-list()
top_1_neutral<-map(neutral_list,top_1_func)%>%bind_rows(.id = "Value")
names(top_1_neutral)<-(c("simulation_num","Top_1"))
top_1_neutral$model<-"Neutral"
diff_neutral<-join(top_1_neutral, per_sim_neutral)%>%mutate(diff=Top_1-Mean_pall)






####################################

####



plot_pall_mean<-rbind(per_sim_sel_s0.05,per_sim_neutral_n200,per_sim_neutral_n50,per_sim_sel_recomb_d0.1,per_sim_sel_d0.1,per_sim_recomb,per_sim_neutral)

ggplot(plot_pall_mean, aes(x=model, y=Mean_pall, fill=model))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_brewer(palette = "Pastel2")+
  geom_jitter(color="black", size=1.5, alpha=0.9) 


plot_top1_mean<-rbind(n200_list,n50_list,top_1_sel_r_s0.05_list,top_1_sel_s0.05_list,top_1_sel_recomb_d0.1,top_1_sel_d0.1,top_1_recomb,top_1_neutral)

ggplot(plot_top1_mean, aes(x=model, y=Top_1, fill=model))+
  geom_boxplot()+
  theme_bw()+
  geom_jitter(color="black", size=1.5, alpha=0.9)+
  scale_fill_brewer(palette = "Pastel2")
  ggtitle("Top 1% p")

  
  
  
plot_diff<-rbind(diff_neutral,diff_recomb, diff_sel_d0.1, diff_sel_d0.5, diff_recomb_d0.1, diff_recomb_d0.5)

ggplot(plot_diff, aes(x=model, y=diff, fill=model))+
  geom_boxplot()+
  theme_bw()+
  geom_jitter(color="black", size=1.5, alpha=0.9)+
  scale_fill_brewer(palette = "Pastel2")


