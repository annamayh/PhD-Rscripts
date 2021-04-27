
library("purrr")
library(dplyr)
library(plyr)
library(ggplot2)

setwd("H:/")


## function to add pall values to a list containining all iterations 

add_pall<-function(iteration_list){
  iteration_list%>%mutate(pall=UNAFF/100)
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




########################################################################################
#### loding diff models #####

#neutral model (done before ash server)
neutral_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Output/ROH_output/neutral_m6_sim*_ROH_OUT.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)

neutral_list<-map(neutral_list, add_pall)

for (i in 1:15){
  neutral_list[[i]]$simulation_num<-i #adding sim number for grouping later
}

unlisted_neutral<-bind_rows(neutral_list, .id = "simulation_num")#unlisting
mean(unlisted_neutral$pall)

per_sim_neutral<-unlisted_neutral%>%select(simulation_num, pall)%>%
  group_by(simulation_num)%>% #grouping by simulation iteration
  dplyr::summarise(Mean_pall=mean(pall))

per_sim_neutral$model<-"Neutral"

top_1_neutral<-list()
top_1_neutral<-map(neutral_list,top_1_func)%>%bind_rows(.id = "Value")
names(top_1_neutral)<-(c("simulation_num","Top_1"))
mean(top_1_neutral$Top_1)
top_1_neutral$model<-"Neutral"

diff_neutral<-join(top_1_neutral, per_sim_neutral)%>%mutate(diff=Top_1-Mean_pall)




#recombination only
recomb_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Sel_recomb/sel_recomb_*/sel_recomb_ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
recomb_list<-map(recomb_list, add_pall)
for (i in 1:14){
  recomb_list[[i]]$simulation_num<-i #adding sim number for grouping later
 }

unlisted_recomb<-bind_rows(recomb_list, .id = "simulation_num")#unlisting
mean(unlisted_recomb$pall)

per_sim_recomb<-unlisted_recomb%>%select(simulation_num, pall)%>%
  group_by(simulation_num)%>% #grouping by simulation iteration
  dplyr::summarise(Mean_pall=mean(pall))

per_sim_recomb$model<-"Recomb"

top_1_recomb<-list()
top_1_recomb<-map(recomb_list,top_1_func)%>%bind_rows(.id = "Value")
names(top_1_recomb)<-(c("simulation_num","Top_1"))
mean(top_1_recomb$Top_1)
top_1_recomb$model<-"Recomb"

diff_recomb<-join(top_1_recomb, per_sim_recomb)%>%mutate(diff=Top_1-Mean_pall)


#both selections d=0.1
sel_d0.1_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Sel_d0.1/sel_d0.1_*/sel_d0.1ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
sel_d0.1_list<-map(sel_d0.1_list, add_pall)

for (i in 1:13){
  sel_d0.1_list[[i]]$simulation_num<-i #adding sim number for grouping later
 }

unlisted_sel_d0.1<-bind_rows(sel_d0.1_list, .id = "simulation_num")#unlisting
mean(unlisted_sel_d0.1$pall)

per_sim_sel_d0.1<-unlisted_sel_d0.1%>%select(simulation_num, pall)%>%
  group_by(simulation_num)%>% #grouping by simulation iteration
  dplyr::summarise(Mean_pall=mean(pall))

per_sim_sel_d0.1$model<-"Sel_d0.1"

top_1_sel_d0.1<-list()
top_1_sel_d0.1<-map(sel_d0.1_list,top_1_func)%>%bind_rows(.id = "Value")
names(top_1_sel_d0.1)<-(c("simulation_num","Top_1"))
mean(top_1_sel_d0.1$Top_1)
top_1_sel_d0.1$model<-"Sel_d0.1"

diff_sel_d0.1<-join(top_1_sel_d0.1, per_sim_sel_d0.1)%>%mutate(diff=Top_1-Mean_pall)




#selection plus recobination, d=0.1
sel_reomb_d0.1_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Sel_recomb_d0.1/sel_recomb_d0.1_*/sel_recomb_d0.1ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)

sel_reomb_d0.1_list<-map(sel_reomb_d0.1_list, add_pall)

for (i in 1:13){
  sel_reomb_d0.1_list[[i]]$simulation_num<-i #adding sim number for grouping later
}

unlisted_sel_recomb_d0.1<-bind_rows(sel_reomb_d0.1_list, .id = "simulation_num")#unlisting
mean(unlisted_sel_recomb_d0.1$pall)

per_sim_sel_recomb_d0.1<-unlisted_sel_recomb_d0.1%>%select(simulation_num, pall)%>%
  group_by(simulation_num)%>% #grouping by simulation iteration
  dplyr::summarise(Mean_pall=mean(pall))

per_sim_sel_recomb_d0.1$model<-"Sel_recomb_d0.1"

top_1_sel_recomb_d0.1<-list()
top_1_sel_recomb_d0.1<-map(sel_reomb_d0.1_list,top_1_func)%>%bind_rows(.id = "Value")
names(top_1_sel_recomb_d0.1)<-(c("simulation_num","Top_1"))
mean(top_1_sel_recomb_d0.1$Top_1)
top_1_sel_recomb_d0.1$model<-"Sel_recomb_d0.1"

diff_recomb_d0.1<-join(top_1_sel_recomb_d0.1, per_sim_sel_recomb_d0.1)%>%mutate(diff=Top_1-Mean_pall)


#selection plus recombination d=0.5
sel_s0.05_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Sel_s0.05d0.1/sel_s0.05d0.1_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
sel_s0.05_list<-map(sel_s0.05_list, add_pall)
for (i in 1:23){
  sel_s0.05_list[[i]]$simulation_num<-i #adding sim number for grouping later
  }

unlisted_sel_s0.05_list<-bind_rows(sel_s0.05_list, .id = "simulation_num")#unlisting
mean(unlisted_sel_s0.05_list$pall)

per_sim_sel_s0.05<-unlisted_sel_s0.05_list%>%select(simulation_num, pall)%>%
  group_by(simulation_num)%>% #grouping by simulation iteration
  dplyr::summarise(Mean_pall=mean(pall))

per_sim_sel_s0.05$model<-"Sel_s0.05"

top_1_sel_s0.05_list<-list()
top_1_sel_s0.05_list<-map(sel_s0.05_list,top_1_func)%>%bind_rows(.id = "Value")
names(top_1_sel_s0.05_list)<-(c("simulation_num","Top_1"))
mean(top_1_sel_s0.05_list$Top_1)
top_1_sel_s0.05_list$model<-"Sel_s0.05"

#####
sel_r_s0.05_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Sel_r_s0.05d0.1/sel_r_s0.05d0.1_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
sel_r_s0.05_list<-map(sel_r_s0.05_list, add_pall)
for (i in 1:20){
  sel_r_s0.05_list[[i]]$simulation_num<-i #adding sim number for grouping later
}

unlisted_sel_r_s0.05_list<-bind_rows(sel_r_s0.05_list, .id = "simulation_num")#unlisting
mean(unlisted_sel_r_s0.05_list$pall)

per_sim_sel_r_s0.05<-unlisted_sel_r_s0.05_list%>%select(simulation_num, pall)%>%
  group_by(simulation_num)%>% #grouping by simulation iteration
  dplyr::summarise(Mean_pall=mean(pall))

per_sim_sel_s0.05$model<-"Sel_r_s0.05"

top_1_sel_r_s0.05_list<-list()
top_1_sel_r_s0.05_list<-map(sel_r_s0.05_list,top_1_func)%>%bind_rows(.id = "Value")
names(top_1_sel_r_s0.05_list)<-(c("simulation_num","Top_1"))
mean(top_1_sel_r_s0.05_list$Top_1)
top_1_sel_r_s0.05_list$model<-"Sel_r_s0.05"
####


neutral_n200_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Neutral_Ne200/neutral_Ne200_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
neutral_n200_list<-map(neutral_n200_list, add_pall)
for (i in 1:23){
  neutral_n200_list[[i]]$simulation_num<-i #adding sim number for grouping later
}

unlisted_neutral_n200<-bind_rows(neutral_n200_list, .id = "simulation_num")#unlisting
mean(unlisted_neutral_n200$pall)

per_sim_neutral_n200<-unlisted_neutral_n200%>%select(simulation_num, pall)%>%
  group_by(simulation_num)%>% #grouping by simulation iteration
  dplyr::summarise(Mean_pall=mean(pall))

per_sim_neutral_n200$model<-"N200"

n200_list<-list()
n200_list<-map(neutral_n200_list,top_1_func)%>%bind_rows(.id = "Value")
names(n200_list)<-(c("simulation_num","Top_1"))
mean(n200_list$Top_1)
n200_list$model<-"n200"

#####


neutral_n50_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Neutral_Ne50/neutral_Ne50_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
neutral_n50_list<-map(neutral_n50_list, add_pall)
for (i in 1:25){
  neutral_n50_list[[i]]$simulation_num<-i #adding sim number for grouping later
}

unlisted_neutral_n50<-bind_rows(neutral_n50_list, .id = "simulation_num")#unlisting
mean(unlisted_neutral_n50$pall)

per_sim_neutral_n50<-unlisted_neutral_n50%>%select(simulation_num, pall)%>%
  group_by(simulation_num)%>% #grouping by simulation iteration
  dplyr::summarise(Mean_pall=mean(pall))

per_sim_neutral_n50$model<-"n50"

n50_list<-list()
n50_list<-map(neutral_n50_list,top_1_func)%>%bind_rows(.id = "Value")
names(n50_list)<-(c("simulation_num","Top_1"))
mean(n50_list$Top_1)
n50_list$model<-"n50"



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


