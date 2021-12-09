
library("purrr")
library(dplyr)
library(plyr)
library(ggplot2)
library(forcats)
library(wesanderson)
library(patchwork)

setwd("H:/")

######################################################################################################################################
######## load all simulation files in #####
## neutral model
neutral_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Ne7500/Rum/Neutral/neutral*/*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(neutral_list)){
  neutral_list[[i]]$simulation_num<-i
  neutral_list[[i]]$Model<-"Neutral"
  neutral_list[[i]]$PopHistory<-"Rum"#adding sim number for grouping later
} 
#selection and recombination
sel_reomb_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Ne7500/Rum/NewH_sel_r/sel_r_*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(sel_reomb_list)){
  sel_reomb_list[[i]]$simulation_num<-i
  sel_reomb_list[[i]]$Model<-"Varied\nrecombination\n+ selection"
  sel_reomb_list[[i]]$PopHistory<-"Rum"
} 
## stronger selection coefficient 
sel_s0.05_list <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Ne7500/Rum/NewH_strongsel_r_noMAF/ROH_noMAF*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(sel_s0.05_list)){
  sel_s0.05_list[[i]]$simulation_num<-i
  sel_s0.05_list[[i]]$Model<-"Varied\nrecombination\n+stronger\nselection"
  sel_s0.05_list[[i]]$PopHistory<-"Rum"#adding sim number for grouping later
}
#############################
## servere bottleneck neutral 
##############################
bottleneck <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Ne7500/btl/Neutral/neutral*/*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck)){
  bottleneck[[i]]$simulation_num<-i #adding sim number for grouping later
  bottleneck[[i]]$Model<-"Neutral"
  bottleneck[[i]]$PopHistory<-"Severe bottleneck"
}
##bottlneck s+r
bottleneck_s_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Ne7500/btl/NewH_sel_r/sel_r_*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_s_r)){
  bottleneck_s_r[[i]]$simulation_num<-i #adding sim number for grouping later
  bottleneck_s_r[[i]]$Model<-"Varied\nrecombination\n+ selection"
  bottleneck_s_r[[i]]$PopHistory<-"Severe bottleneck"
}
## bottleneck s = 0.05
bottleneck_s0.05_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Ne7500/btl/NewH_strongsel_r_noMAF/ROH_noMAF*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(bottleneck_s0.05_r)){
  bottleneck_s0.05_r[[i]]$simulation_num<-i #adding sim number for grouping later
  bottleneck_s0.05_r[[i]]$Model<-"Varied\nrecombination\n+stronger\nselection"
  bottleneck_s0.05_r[[i]]$PopHistory<-"Severe bottleneck"
}
##########################
## no bottleneck neutral##
##########################
Nobottleneck <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Ne7500/NObtl/NObtl/Neutral/neutral_*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(Nobottleneck)){
  Nobottleneck[[i]]$simulation_num<-i #adding sim number for grouping later
  Nobottleneck[[i]]$Model<-"Neutral"
  Nobottleneck[[i]]$PopHistory<-"No bottleneck"
}
##NObottlneck s+r
setwd("C:/Users/s1881212/Ash_server/Output/ROH_output")
NObottleneck_s_r <- lapply(Sys.glob("NewH_sel_r/sel_r*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(NObottleneck_s_r)){
  NObottleneck_s_r[[i]]$simulation_num<-i #adding sim number for grouping later
  NObottleneck_s_r[[i]]$Model<-"Varied\nrecombination\n+ selection"
  NObottleneck_s_r[[i]]$PopHistory<-"No bottleneck"
}
## NObottleneck s = 0.05
setwd("H:/")
NObottleneck_s0.05_r <- lapply(Sys.glob("PHD_2ndYR/Slim/Ash_server_output/Ne7500/NObtl/NObtl/NewH_strongsel_r_noMAF/ROH_noMAF*/ROH_out*.hom.summary"), read.table, header = TRUE, stringsAsFactors=FALSE)
for (i in 1:length(NObottleneck_s0.05_r)){
  NObottleneck_s0.05_r[[i]]$simulation_num<-i #adding sim number for grouping later
  NObottleneck_s0.05_r[[i]]$Model<-"Varied\nrecombination\n+stronger\nselection"
  NObottleneck_s0.05_r[[i]]$PopHistory<-"No bottleneck"
}

  
  
################################################################################################################
###### FUNCTIONS ######################
#######################################

add_pall<-function(iteration_list){
  iteration_list%>%mutate(pall=UNAFF/100)
}

add_pall_no_btl<-function(iteration_list){
  iteration_list%>%mutate(pall=UNAFF/7500)
}

## function to get the mean pall per simulation for the unlisted dataframe
unlist_then_get_mean_pall_per_sim<-function(listed_sims){
  bind_rows(listed_sims, .id = "simulation_num")%>%
  select(simulation_num, pall, Model, PopHistory)%>%
    group_by(simulation_num)%>% #grouping by simulation iteration
    dplyr::summarise(Model=first(Model),PopHistory = first(PopHistory),Mean_pall=mean(pall))#getting mean pall of simulation
 }


max_prop_ROH<-function(listed_sims){
  bind_rows(listed_sims, .id = "simulation_num")%>%
    select(simulation_num, pall, Model, PopHistory)%>%
    group_by(simulation_num)%>% #grouping by simulation iteration
    dplyr::summarise(Model=first(Model),PopHistory = first(PopHistory),Max_pall=max(pall))
 
}

#func to get top 1% in a datafram
top_1_func<-function (x){quantile(x$pall, c(.01,.99))}

#function using function above to apply to all dataframes in a list of iterations
get_top1_pall<-function(iteration_list){
  
  map(iteration_list,top_1_func)%>%bind_rows(.id = "Value")%>%
    setNames(.,c("simulation_num","Bottom_1","Top_1"))
 }


###########################################################################################

#save.image(file="PHD_2ndYR/Slim/List_SLiM_output_new_sims.RData") #saved on 2/12/21

#################################################
########### LOAD IN DATA ########################
#load("PHD_2ndYR/Slim/List_all_SLiM_output.RData")
#################################################

###   SET UP FOR PLOT 1 ##########
### getting mean pall for plot ######

models_minus<-list(neutral_list,sel_reomb_list,sel_s0.05_list,bottleneck,bottleneck_s_r,bottleneck_s0.05_r)
names(models_minus)<-c("Rum_neutral","Rum_sel_reomb","Rum_s0.05_recomb","btl_neutral","btl_sel_recomb","btl_sel0.05_recomb")

#need to do pop size =100 and =7500 separately 
get_propROH<-function(list){
  list%>%map(add_pall)%>%unlist_then_get_mean_pall_per_sim}

prop_wROH_minus<-models_minus%>%map(get_propROH)

### now models with 7500 Ne 
nobtl_models<-list(Nobottleneck, NObottleneck_s_r, NObottleneck_s0.05_r)
names(nobtl_models)<-c("Nobtl_neutral", "Nobtl_sel_recomb", "Nobtl_sel0.05_recomb")

get_propROH_btl<-function(list){
  list%>%map(add_pall_no_btl)%>%unlist_then_get_mean_pall_per_sim}

prop_wROH_btl<-nobtl_models%>%map(get_propROH_btl)

## merging 2 lists together 
prop_wROH_all<-c(prop_wROH_minus, prop_wROH_btl)




### SET UP FOR PLOT 2 ##########
##working out maximum pall value#####

get_max_propROH_minus<-function(list){
  list%>%map(add_pall)%>%max_prop_ROH()
}
max_minus<-models_minus%>%map(get_max_propROH_minus)

get_max_propROH_no_btl<-function(list){
  list%>%map(add_pall_no_btl)%>%max_prop_ROH()
  }

max_nobtl<-nobtl_models%>%map(get_max_propROH_no_btl)

max_all<-c(max_nobtl, max_minus)




######################################################################################


## PLOT 1 ####
mean_pall_plot<-do.call(rbind.data.frame, prop_wROH_all)

##inset plot 
no_bottle<-rbind(pall_nobtl,pall_nobtl_s_r,pall_nobtl_s0.05_r)

## This function allows us to specify which facet to annotate
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}

mean_ann_text <- data.frame(Model = "Varied\nrecombination\n+ selection",Mean_pall = 0.0755,lab = "Rum actual value",
                       PopHistory = factor("No bottleneck",levels = c("Rum","Severe bottleneck","No bottleneck")))

mean_ann_text_kint <- data.frame(Model = "Varied\nrecombination\n+ selection",Mean_pall = 0.0458,lab = "Kintyre actual value",
                            PopHistory = factor("No bottleneck",levels = c("Rum","Severe bottleneck","No bottleneck")))

Mean_plot<-mean_pall_plot%>%
  mutate(Model = fct_relevel(Model, "Neutral","Varied\nrecombination\n+ selection","Varied\nrecombination\n+stronger\nselection"))%>%
  mutate(PopHistory = fct_relevel(PopHistory,"Severe bottleneck","Rum","No bottleneck"))%>%
  ggplot(aes(x=Model, y=Mean_pall, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  xlab("Simulation Model")+
  ylab("Mean ROH density per simulation")+
  labs(tag = "A")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    axis.text.y = element_text(face="bold"),
    axis.text.x = element_text(face="bold")
    
  )+
  facet_wrap(~PopHistory)+
  theme(legend.position = "none")+
  stat_summary(fun=mean, geom="point",size=2)+
  geom_hline(yintercept = 0.064, linetype=4, color = "black", size = 1)+
  geom_hline(yintercept = 0.034, linetype=5, color = "black", size = 1)+
  geom_text(data = mean_ann_text,label = "Rum actual value")+
  geom_text(data = mean_ann_text_kint,label = "Kintyre actual value")
  annotation_custom2(grob=ggplotGrob(inset_mean), 
                     data = data.frame(Model="No bottleneck", Mean_pall=0.2),
                     ymin = 0.2, ymax=0.4, xmin="Varied\nRecomb", xmax="Varied\nrecomb\n+stronger\nselection")


## inset plot 
no_bottle%>%
  mutate(Model = fct_relevel(Model, "Neutral","Varied\nrecombination\n+ selection","Varied\nrecombination\n+stronger\nselection"))%>%
  ggplot(aes(x=Model, y=Mean_pall, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_classic()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text = element_text(size=20)
    
  )+
  theme(legend.position = "none")+
  stat_summary(fun=mean, geom="point",size=2)
  

###########################################################################################################
## PLOT 2 ######
### no fucking clue 

max_pall_plot<-do.call(rbind.data.frame, max_all)

## This function allows us to specify which facet to annotate
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}

mean_ann_text <- data.frame(Model = "Varied\nrecombination\n+ selection",Mean_pall = 0.0755,lab = "Rum actual value",
                            PopHistory = factor("No bottleneck",levels = c("Rum","Severe bottleneck","No bottleneck")))

mean_ann_text_kint <- data.frame(Model = "Varied\nrecombination\n+ selection",Mean_pall = 0.0458,lab = "Kintyre actual value",
                                 PopHistory = factor("No bottleneck",levels = c("Rum","Severe bottleneck","No bottleneck")))
max_pall_plot$Max_pall

Max_plot<-max_pall_plot%>%
  mutate(Model = fct_relevel(Model, "Neutral","Varied\nrecombination\n+ selection","Varied\nrecombination\n+stronger\nselection"))%>%
  mutate(PopHistory = fct_relevel(PopHistory,"Severe bottleneck","Rum","No bottleneck"))%>%
  ggplot(aes(x=Model,y=Max_pall, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  xlab("Simulation Model")+
  ylab("Mean ROH density per simulation")+
  labs(tag = "A")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    axis.text.y = element_text(face="bold"),
    axis.text.x = element_text(face="bold")
    
  )+
  facet_wrap(~PopHistory)+
  theme(legend.position = "none")+
  stat_summary(fun=mean, geom="point",size=2)+
  geom_hline(yintercept = 0.064, linetype=4, color = "black", size = 1)+
  geom_hline(yintercept = 0.034, linetype=5, color = "black", size = 1)+
  geom_text(data = mean_ann_text,label = "Rum actual value")+
  geom_text(data = mean_ann_text_kint,label = "Kintyre actual value")
annotation_custom2(grob=ggplotGrob(inset_mean), 
                   data = data.frame(Model="No bottleneck", Mean_pall=0.2),
                   ymin = 0.2, ymax=0.4, xmin="Varied\nRecomb", xmax="Varied\nrecomb\n+stronger\nselection")


## inset plot 
no_bottle%>%
  mutate(Model = fct_relevel(Model, "Neutral","Varied\nrecombination\n+ selection","Varied\nrecombination\n+stronger\nselection"))%>%
  ggplot(aes(x=Model, y=Mean_pall, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_classic()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text = element_text(size=20)
    
  )+
  theme(legend.position = "none")+
  stat_summary(fun=mean, geom="point",size=2)



###################################################################################################################################

top_ann_text <- data.frame(Model = "Varied\nrecombination\n+ selection",Top_1 = 0.17,lab = "Rum actual value",
                            Pop_History = factor("No bottleneck",levels = c("Rum","Severe bottleneck","No bottleneck")))

top_ann_text_kint <- data.frame(Model = "Varied\nrecombination\n+ selection",Top_1 = 0.105,lab = "Kintyre actual value",
                                 Pop_History = factor("No bottleneck",levels = c("Rum","Severe bottleneck","No bottleneck")))

## plotting top 1%
Top_1<-plot_pall%>%
  mutate(Model = fct_relevel(Model, "Neutral","Varied\nrecombination\n+ selection","Varied\nrecombination\n+stronger\nselection"))%>%
  mutate(Pop_History = fct_relevel(Pop_History,"Severe bottleneck","Rum","No bottleneck"))%>%
  ggplot(aes(x=Model, y=Top_1, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  xlab("Simulation Model")+
  ylab("ROH hotspot threshold")+
  labs(tag = "B")+
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold")
    
  )+
  facet_wrap(~Pop_History)+
  theme(legend.position = "none")+
  stat_summary(fun=mean, geom="point",size=2)+ # adds the mean dot 
  geom_hline(yintercept = 0.149, linetype=4, color = "black", size = 1)+
  geom_hline(yintercept = 0.083, linetype=5, color = "black", size = 1)+
  geom_text(data = top_ann_text,label = "Rum actual value")+
  geom_text(data = top_ann_text_kint,label = "Kintyre actual value")

no_bottle%>%
  mutate(Model = fct_relevel(Model, "Neutral","Varied\nrecombination\n+ selection","Varied\nrecombination\n+stronger\nselection"))%>%
  ggplot(aes(x=Model, y=Top_1, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_classic()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text = element_text(size=20)
    
  )+
  theme(legend.position = "none")+
  stat_summary(fun=mean, geom="point",size=2)


Mean_plot/Top_1

####################################################################################################

#####################################################################################################

#####################################################################################################
### plotting bottom 1%

plot_pall%>%
  mutate(Model = fct_relevel(Model, "Neutral","Selection","Varied\nRecomb","Varied\nrecomb\n+ selection","Higher\nselection\ncoefficients"))%>%
  mutate(Pop_History = fct_relevel(Pop_History,"Rum","Severe bottleneck","No bottleneck"))%>%
  ggplot(aes(x=Model, y=Bottom_1, fill=Model))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=5))+
  xlab("Simulation Model")+
  ylab("Bottom 1 % with ROH at a SNP per simulation")+
  theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(face="bold")
    
  )+
  facet_wrap(~Pop_History)+
  theme(legend.position = "none")+
  stat_summary(fun.y=mean, geom="point",size=2)


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
  theme(legend.position = "none")+
  stat_summary(fun=mean, geom="point",size=2)

