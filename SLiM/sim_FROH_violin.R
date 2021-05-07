library("purrr")
library(dplyr)
library(plyr)
library(ggplot2)
library(forcats)
library(wesanderson)

setwd("H:/")

get_FROH<-function(iteration_list){
  iteration_list%>%mutate(FROH=KB/100000)
}

mean_FROH_persim<-function(list_sims){bind_rows(list_sims, .id = "simulation_num")%>%
  select(simulation_num, FROH)%>%
  group_by(simulation_num)%>% #grouping by simulation iteration
  dplyr::summarise(Mean_FROH=mean(FROH))}

## load files

neutral_list_FROH <- lapply(Sys.glob("PHD_2ndYR/Slim/Output/ROH_output/neutral_m6_sim*_ROH_OUT.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(neutral_list_FROH)){
  neutral_list_FROH[[i]]$simulation_num<-i #adding sim number for grouping later
} 

recomb_list_FROH <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Sel_recomb/sel_recomb_*/sel_recomb_ROH_out*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(recomb_list_FROH)){
  recomb_list_FROH[[i]]$simulation_num<-i #adding sim number for grouping later
} 

sel_list_FROH <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Sel_d0.1/sel_d0.1_*/sel_d0.1ROH_out*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(sel_list_FROH)){
  sel_list_FROH[[i]]$simulation_num<-i #adding sim number for grouping later
} 

sel_reomb_list_FROH <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Sel_recomb_d0.1/sel_recomb_d0.1_*/sel_recomb_d0.1ROH_out*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(sel_reomb_list_FROH)){
  sel_reomb_list_FROH[[i]]$simulation_num<-i #adding sim number for grouping later
} 

## stronger selection coefficient 
sel_s0.05_list_FROH <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Sel_r_s0.05d0.1/sel_r_s0.05d0.1_*/ROH_out_sim*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(sel_s0.05_list_FROH)){
  sel_s0.05_list_FROH[[i]]$simulation_num<-i #adding sim number for grouping later
}

## servere bottleneck neutral 
bottleneck_FROH <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Neutral_bottleneck/Neutral_btl_*/ROH_out_sim*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_FROH)){
  bottleneck_FROH[[i]]$simulation_num<-i #adding sim number for grouping later
}

## bottleneck selection
bottleneck_sel_FROH <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Btl_sel/Btl_sel_*/ROH_out_sim*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_sel_FROH)){
  bottleneck_sel_FROH[[i]]$simulation_num<-i #adding sim number for grouping later
}

## bottleneck recomb
bottleneck_recomb_FROH <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Btl_recomb/Btl_recomb_*/ROH_out_sim*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_recomb_FROH)){
  bottleneck_recomb_FROH[[i]]$simulation_num<-i #adding sim number for grouping later
}

##bottlneck s+r
bottleneck_s_r_FROH <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Btl_sel_recomb/Btl_sel_r_*/ROH_out_sim*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_s_r_FROH)){
  bottleneck_s_r_FROH[[i]]$simulation_num<-i #adding sim number for grouping later
}

## bottleneck s = 0.05
bottleneck_s0.05_r_FROH <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Btl_sel0.05_recomb/Btl_sel0.05_r_*/ROH_out_sim*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_s0.05_r_FROH)){
  bottleneck_s0.05_r_FROH[[i]]$simulation_num<-i #adding sim number for grouping later
}

## no bottleneck neutral 
Nobottleneck_FROH <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Neutral_NObottleneck/Neutral_nobtl_*/ROH_out_sim*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(Nobottleneck_FROH)){
  Nobottleneck_FROH[[i]]$simulation_num<-i #adding sim number for grouping later
}


## NO bottleneck selection
NObottleneck_sel_FROH <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/NOBtl_sel/NOBtl_sel_*/ROH_out_sim*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(NObottleneck_sel_FROH)){
  NObottleneck_sel_FROH[[i]]$simulation_num<-i #adding sim number for grouping later
}

## NO bottleneck recomb
NObottleneck_recomb_FROH <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/NOBtl_recomb/NOBtl_recomb_*/ROH_out_sim*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(NObottleneck_recomb_FROH)){
  NObottleneck_recomb_FROH[[i]]$simulation_num<-i #adding sim number for grouping later
}

##NObottlneck s+r
NObottleneck_s_r_FROH <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/NOBtl_sel_r/NOBtl_sel_r_*/ROH_out_sim*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(NObottleneck_s_r_FROH)){
  NObottleneck_s_r_FROH[[i]]$simulation_num<-i #adding sim number for grouping later
}

## NObottleneck s = 0.05
NObottleneck_s0.05_r_FROH <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/NOBtl_sel0.05_r/NOBtl_sel0.05_r_*/ROH_out_sim*.hom.indiv"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(NObottleneck_s0.05_r_FROH)){
  NObottleneck_s0.05_r_FROH[[i]]$simulation_num<-i #adding sim number for grouping later
}




FROH_neutral<-neutral_list_FROH%>%map(get_FROH)%>%mean_FROH_persim()
FROH_recomb<-recomb_list_FROH%>%map(get_FROH)%>%mean_FROH_persim()
FROH_sel<-sel_list_FROH%>%map(get_FROH)%>%mean_FROH_persim()
FROH_s_recomb<-sel_reomb_list_FROH%>%map(get_FROH)%>%mean_FROH_persim()
FROH_s0.05<-sel_s0.05_list_FROH%>%map(get_FROH)%>%mean_FROH_persim()

FROH_neutral$Model<-"Neutral"
FROH_recomb$Model<-"Varied\nRecomb"
FROH_sel$Model<-"Selection"
FROH_s_recomb$Model<-"Varied\nrecomb\n+ selection"
FROH_s0.05$Model<-"Higher\nselection\ncoefficients"

FROH_neutral$Pop_History<-"Rum"
FROH_recomb$Pop_History<-"Rum"
FROH_sel$Pop_History<-"Rum"
FROH_s_recomb$Pop_History<-"Rum"
FROH_s0.05$Pop_History<-"Rum"

FROH_Btlneutral<-bottleneck_FROH%>%map(get_FROH)%>%mean_FROH_persim()
FROH_Btlsel<-bottleneck_sel_FROH%>%map(get_FROH)%>%mean_FROH_persim()
FROH_Btlrecomb<-bottleneck_recomb_FROH%>%map(get_FROH)%>%mean_FROH_persim()
FROH_Btls_recomb<-bottleneck_s_r_FROH%>%map(get_FROH)%>%mean_FROH_persim()
FROH_Btls0.05<-bottleneck_s0.05_r_FROH%>%map(get_FROH)%>%mean_FROH_persim()

FROH_Btlneutral$Model<-"Neutral"
FROH_Btlsel$Model<-"Selection"
FROH_Btlrecomb$Model<-"Varied\nRecomb"
FROH_Btls_recomb$Model<-"Varied\nrecomb\n+ selection"
FROH_Btls0.05$Model<-"Higher\nselection\ncoefficients"

FROH_Btlneutral$Pop_History<-"Severe bottleneck"
FROH_Btlsel$Pop_History<-"Severe bottleneck"
FROH_Btlrecomb$Pop_History<-"Severe bottleneck"
FROH_Btls_recomb$Pop_History<-"Severe bottleneck"
FROH_Btls0.05$Pop_History<-"Severe bottleneck"

FROH_NObtlneutral<-Nobottleneck_FROH%>%map(get_FROH)%>%mean_FROH_persim()
FROH_NObtlsel<-NObottleneck_sel_FROH%>%map(get_FROH)%>%mean_FROH_persim()
FROH_NObtlrecomb<-NObottleneck_recomb_FROH%>%map(get_FROH)%>%mean_FROH_persim()
FROH_NObtls_recomb<-NObottleneck_s_r_FROH%>%map(get_FROH)%>%mean_FROH_persim()
FROH_NObtls0.05<-NObottleneck_s0.05_r_FROH%>%map(get_FROH)%>%mean_FROH_persim()

FROH_NObtlneutral$Model<-"Neutral"
FROH_NObtlsel$Model<-"Selection"
FROH_NObtlrecomb$Model<-"Varied\nRecomb"
FROH_NObtls_recomb$Model<-"Varied\nrecomb\n+ selection"
FROH_NObtls0.05$Model<-"Higher\nselection\ncoefficients"

FROH_NObtlneutral$Pop_History<-"No bottleneck"
FROH_NObtlsel$Pop_History<-"No bottleneck"
FROH_NObtlrecomb$Pop_History<-"No bottleneck"
FROH_NObtls_recomb$Pop_History<-"No bottleneck"
FROH_NObtls0.05$Pop_History<-"No bottleneck"


####

plot_FROH<-rbind(FROH_neutral,FROH_recomb,FROH_sel,FROH_s_recomb,FROH_s0.05,
                 FROH_Btlneutral,FROH_Btlsel,FROH_Btlrecomb,FROH_Btls_recomb,FROH_Btls0.05,
                 FROH_NObtlneutral,FROH_NObtlsel,FROH_NObtlrecomb,FROH_NObtls_recomb,FROH_NObtls0.05)



plot_FROH%>%
  mutate(Model = fct_relevel(Model, "Neutral","Selection","Varied\nRecomb","Varied\nrecomb\n+ selection","Higher\nselection\ncoefficients"))%>%
  mutate(Pop_History = fct_relevel(Pop_History,"Rum","Severe bottleneck","No bottleneck"))%>%
  ggplot(aes(x=Model, y=Mean_FROH, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  xlab("Simulation Model")+
  ylab("Mean FROH per simulation")+
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold")
    
  )+
  facet_wrap(~Pop_History)+
  theme(legend.position = "none")

