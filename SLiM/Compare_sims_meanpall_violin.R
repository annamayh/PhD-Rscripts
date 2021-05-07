
library("purrr")
library(dplyr)
library(plyr)
library(ggplot2)
library(forcats)
library(wesanderson)

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


#func to get top 1% in a datafram
top_1_func<-function (x){quantile(x$pall, c(.99))}

#function using function above to apply to all dataframes in a list of iterations
get_top1_pall<-function(iteration_list){
  
  map(iteration_list,top_1_func)%>%bind_rows(.id = "Value")%>%
    setNames(.,c("simulation_num","Top_1"))
    
  
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
sel_s0.05_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Sel_r_s0.05d0.1/sel_r_s0.05d0.1_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(sel_s0.05_list)){
  sel_s0.05_list[[i]]$simulation_num<-i #adding sim number for grouping later
}

## servere bottleneck neutral 
bottleneck <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Neutral_bottleneck/Neutral_btl_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck)){
  bottleneck[[i]]$simulation_num<-i #adding sim number for grouping later
}

## bottleneck selection
bottleneck_sel <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Btl_sel/Btl_sel_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_sel)){
  bottleneck_sel[[i]]$simulation_num<-i #adding sim number for grouping later
}

## bottleneck recomb
bottleneck_recomb <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Btl_recomb/Btl_recomb_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_recomb)){
  bottleneck_recomb[[i]]$simulation_num<-i #adding sim number for grouping later
}

##bottlneck s+r
bottleneck_s_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Btl_sel_recomb/Btl_sel_r_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_s_r)){
  bottleneck_s_r[[i]]$simulation_num<-i #adding sim number for grouping later
}

## bottleneck s = 0.05
bottleneck_s0.05_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Btl_sel0.05_recomb/Btl_sel0.05_r_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_s0.05_r)){
  bottleneck_s0.05_r[[i]]$simulation_num<-i #adding sim number for grouping later
}

## no bottleneck neutral 
Nobottleneck <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Neutral_NObottleneck/Neutral_nobtl_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(Nobottleneck)){
  Nobottleneck[[i]]$simulation_num<-i #adding sim number for grouping later
}

 
## NO bottleneck selection
NObottleneck_sel <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/NOBtl_sel/NOBtl_sel_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(NObottleneck_sel)){
  NObottleneck_sel[[i]]$simulation_num<-i #adding sim number for grouping later
}

## NO bottleneck recomb
NObottleneck_recomb <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/NOBtl_recomb/NOBtl_recomb_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(NObottleneck_recomb)){
  NObottleneck_recomb[[i]]$simulation_num<-i #adding sim number for grouping later
}

##NObottlneck s+r
NObottleneck_s_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/NOBtl_sel_r/NOBtl_sel_r_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(NObottleneck_s_r)){
  NObottleneck_s_r[[i]]$simulation_num<-i #adding sim number for grouping later
}

## NObottleneck s = 0.05
NObottleneck_s0.05_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/NOBtl_sel0.05_r/NOBtl_sel0.05_r_*/ROH_out_sim*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(NObottleneck_s0.05_r)){
  NObottleneck_s0.05_r[[i]]$simulation_num<-i #adding sim number for grouping later
}
###########################################################################################

## rum pop history simulations
pall_neutral<-neutral_list%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_recomb<-recomb_list%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_sel<-sel_list%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_sel_recomb<-sel_reomb_list%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_sel_s0.05<-sel_s0.05_list%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim

top1_neutral<-neutral_list%>%map(add_pall)%>%get_top1_pall()
top1_recomb<-recomb_list%>%map(add_pall)%>%get_top1_pall()
top1_sel<-sel_list%>%map(add_pall)%>%get_top1_pall()
top1_sel_recomb<-sel_reomb_list%>%map(add_pall)%>%get_top1_pall()
top1_sel_s0.05<-sel_s0.05_list%>%map(add_pall)%>%get_top1_pall()

pall_neutral<-join(pall_neutral,top1_neutral)%>%mutate(diff=Top_1-Mean_pall)
pall_recomb<-join(pall_recomb,top1_recomb)%>%mutate(diff=Top_1-Mean_pall)
pall_sel<-join(pall_sel,top1_sel)%>%mutate(diff=Top_1-Mean_pall)
pall_sel_recomb<-join(pall_sel_recomb,top1_sel_recomb)%>%mutate(diff=Top_1-Mean_pall)
pall_sel_s0.05<-join(pall_sel_s0.05,top1_sel_s0.05)%>%mutate(diff=Top_1-Mean_pall)

pall_neutral$Model<-"Neutral"
pall_recomb$Model<-"Varied\nRecomb"
pall_sel$Model<-"Selection"
pall_sel_recomb$Model<-"Varied\nrecomb\n+ selection"
pall_sel_s0.05$Model<-"Higher\nselection\ncoefficients"

pall_neutral$Pop_History<-"Rum"
pall_recomb$Pop_History<-"Rum"
pall_sel$Pop_History<-"Rum"
pall_sel_recomb$Pop_History<-"Rum"
pall_sel_s0.05$Pop_History<-"Rum"

### bottleneck simulations
pall_bottleneck<-bottleneck%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_bottleneck_sel<-bottleneck_sel%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_bottleneck_recomb<-bottleneck_recomb%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_bottleneck_s_r<-bottleneck_s_r%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim
pall_bottleneck_s0.05_r<-bottleneck_s0.05_r%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim

Btl_top1_neutral<-bottleneck%>%map(add_pall)%>%get_top1_pall()
Btl_top1_sel<-bottleneck_sel%>%map(add_pall)%>%get_top1_pall()
Btl_top1_recomb<-bottleneck_recomb%>%map(add_pall)%>%get_top1_pall()
Btl_top1_sel_recomb<-bottleneck_s_r%>%map(add_pall)%>%get_top1_pall()
Btl_top1_sel_s0.05<-bottleneck_s0.05_r%>%map(add_pall)%>%get_top1_pall()

pall_bottleneck<-join(pall_bottleneck,Btl_top1_neutral)%>%mutate(diff=Top_1-Mean_pall)
pall_bottleneck_sel<-join(pall_bottleneck_sel,Btl_top1_sel)%>%mutate(diff=Top_1-Mean_pall)
pall_bottleneck_recomb<-join(pall_bottleneck_recomb,Btl_top1_recomb)%>%mutate(diff=Top_1-Mean_pall)
pall_bottleneck_s_r<-join(pall_bottleneck_s_r,Btl_top1_sel_recomb)%>%mutate(diff=Top_1-Mean_pall)
pall_bottleneck_s0.05_r<-join(pall_bottleneck_s0.05_r,Btl_top1_sel_s0.05)%>%mutate(diff=Top_1-Mean_pall)


pall_bottleneck$Model<-"Neutral"
pall_bottleneck_sel$Model<-"Selection"
pall_bottleneck_recomb$Model<-"Varied\nRecomb"
pall_bottleneck_s_r$Model<-"Varied\nrecomb\n+ selection"
pall_bottleneck_s0.05_r$Model<-"Higher\nselection\ncoefficients"

pall_bottleneck$Pop_History<-"Severe bottleneck"
pall_bottleneck_sel$Pop_History<-"Severe bottleneck"
pall_bottleneck_recomb$Pop_History<-"Severe bottleneck"
pall_bottleneck_s_r$Pop_History<-"Severe bottleneck"
pall_bottleneck_s0.05_r$Pop_History<-"Severe bottleneck"

##### no bottleneck simulations
pall_nobtl<-Nobottleneck%>%map(add_pall_2000)%>%unlist_then_get_mean_pall_per_sim
pall_nobtl_sel<-NObottleneck_sel%>%map(add_pall_2000)%>%unlist_then_get_mean_pall_per_sim
pall_nobtl_recomb<-NObottleneck_recomb%>%map(add_pall_2000)%>%unlist_then_get_mean_pall_per_sim
pall_nobtl_s_r<-NObottleneck_s_r%>%map(add_pall_2000)%>%unlist_then_get_mean_pall_per_sim
pall_nobtl_s0.05_r<-NObottleneck_s0.05_r%>%map(add_pall_2000)%>%unlist_then_get_mean_pall_per_sim

NOBtl_top1_neutral<-Nobottleneck%>%map(add_pall)%>%get_top1_pall()
NOBtl_top1_sel<-NObottleneck_sel%>%map(add_pall)%>%get_top1_pall()
NOBtl_top1_recomb<-NObottleneck_recomb%>%map(add_pall)%>%get_top1_pall()
NOBtl_top1_sel_recomb<-NObottleneck_s_r%>%map(add_pall)%>%get_top1_pall()
NOBtl_top1_sel_s0.05<-NObottleneck_s0.05_r%>%map(add_pall)%>%get_top1_pall()

pall_nobtl<-join(pall_nobtl,NOBtl_top1_neutral)%>%mutate(diff=Top_1-Mean_pall)
pall_nobtl_sel<-join(pall_nobtl_sel,NOBtl_top1_sel)%>%mutate(diff=Top_1-Mean_pall)
pall_nobtl_recomb<-join(pall_nobtl_recomb,NOBtl_top1_recomb)%>%mutate(diff=Top_1-Mean_pall)
pall_nobtl_s_r<-join(pall_nobtl_s_r,NOBtl_top1_sel_recomb)%>%mutate(diff=Top_1-Mean_pall)
pall_nobtl_s0.05_r<-join(pall_nobtl_s0.05_r,NOBtl_top1_sel_s0.05)%>%mutate(diff=Top_1-Mean_pall)

pall_nobtl$Model<-"Neutral"
pall_nobtl_sel$Model<-"Selection"
pall_nobtl_recomb$Model<-"Varied\nRecomb"
pall_nobtl_s_r$Model<-"Varied\nrecomb\n+ selection"
pall_nobtl_s0.05_r$Model<-"Higher\nselection\ncoefficients"

pall_nobtl$Pop_History<-"No bottleneck"
pall_nobtl_sel$Pop_History<-"No bottleneck"
pall_nobtl_recomb$Pop_History<-"No bottleneck"
pall_nobtl_s_r$Pop_History<-"No bottleneck"
pall_nobtl_s0.05_r$Pop_History<-"No bottleneck"





plot_pall<-rbind(pall_neutral,pall_recomb,pall_sel,pall_sel_recomb,pall_sel_s0.05,
                      pall_bottleneck,pall_bottleneck_sel,pall_bottleneck_recomb,pall_bottleneck_s_r,pall_bottleneck_s0.05_r,
                      pall_nobtl,pall_nobtl_sel,pall_nobtl_recomb,pall_nobtl_s_r,pall_nobtl_s0.05_r)


plot_pall%>%
  mutate(Model = fct_relevel(Model, "Neutral","Selection","Varied \nRecomb","Varied \nrecomb \n+ selection","Higher \nselection \ncoefficients"))%>%
  mutate(Pop_History = fct_relevel(Pop_History,"Rum","Severe bottleneck","No bottleneck"))%>%
  ggplot(aes(x=Model, y=Mean_pall, colour=Pop_History))+
  geom_boxplot()+
  theme_bw()+
  geom_jitter(size=1.5, alpha=0.9)+
  scale_fill_brewer(palette = "Pastel2")+
  xlab("Simulation Model")+
  ylab("Mean proportion of individuals with ROH at a SNP per simulation")+
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold")
    
  )+
  facet_wrap(~Pop_History)
  


plot_pall_mean%>%
  mutate(Model = fct_relevel(Model, "Neutral","Selection","Varied\nRecomb","Varied\nrecomb\n+ selection","Higher\nselection\ncoefficients"))%>%
  mutate(Pop_History = fct_relevel(Pop_History,"Rum","Severe bottleneck","No bottleneck"))%>%
  ggplot(aes(x=Model, y=Mean_pall, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  xlab("Simulation Model")+
  ylab("Mean proportion of individuals with ROH at a SNP per simulation")+
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold")
    
  )+
  facet_wrap(~Pop_History)+
  theme(legend.position = "none")

###################################################################################################################################


## plotting top 1%
plot_pall%>%
  mutate(Model = fct_relevel(Model, "Neutral","Selection","Varied\nRecomb","Varied\nrecomb\n+ selection","Higher\nselection\ncoefficients"))%>%
  mutate(Pop_History = fct_relevel(Pop_History,"Rum","Severe bottleneck","No bottleneck"))%>%
  ggplot(aes(x=Model, y=Top_1, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  xlab("Simulation Model")+
  ylab("Top 1% proportion of individuals with ROH \nat a SNP (hotspot threshold) per simulation")+
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold")
    
  )+
  facet_wrap(~Pop_History)+
  theme(legend.position = "none")


#### plotting difference between top 1% and mean
plot_pall%>%
  mutate(Model = fct_relevel(Model, "Neutral","Selection","Varied\nRecomb","Varied\nrecomb\n+ selection","Higher\nselection\ncoefficients"))%>%
  mutate(Pop_History = fct_relevel(Pop_History,"Rum","Severe bottleneck","No bottleneck"))%>%
  ggplot(aes(x=Model, y=diff, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  xlab("Simulation Model")+
  ylab("Difference between Top 1% and mean proportion of individuals \nwith ROH at a SNP per simulation")+
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold")
    
  )+
  facet_wrap(~Pop_History)+
  theme(legend.position = "none")
