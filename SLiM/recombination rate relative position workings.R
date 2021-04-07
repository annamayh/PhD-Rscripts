
library("purrr")
library(dplyr)
library(ggplot2)


setwd("H:/")


sup_table<-read.table("PHD_2ndYR/Slim/RR_TableS4.txt", header=T, stringsAsFactors = F)

##putting into list per chr
sup_table_list<-list()
for (i in 1:33){
sup_table_list[[i]]<- sup_table%>%filter(CEL.LG==i)  
 }

# working out relative positions
for (j in 1:33){
sup_table_list[[j]]$relative_pos<-NA
for(k in 1:nrow(sup_table_list[[j]])){
  sup_table_list[[j]]$relative_pos[k] <- sup_table_list[[j]]$Window[k]/nrow(sup_table_list[[j]])}#pall per SNP
}

#rounding relative positions to nearest 0.1 
for (m in 1:33){
  sup_table_list[[m]]$relative_pos_rounded<-NA
  for(s in 1:nrow(sup_table_list[[m]])){
    sup_table_list[[m]]$relative_pos_rounded[s] <- plyr::round_any(sup_table_list[[m]]$relative_pos[s], 0.1)
    
    
  }}



relative_pos_table<-bind_rows(sup_table_list)%>%#unlisting
  subset(CEL.LG!=5)%>% #removing LG5
  select(relative_pos_rounded, cM)%>%na.omit()%>%#selecting specfic columns and removing NA values
  group_by(relative_pos_rounded)%>% #grouping by the relative pos to get mean
  dplyr::summarise(mean_block_cM=mean(cM))

#checking in plot
ggplot(relative_pos_table, aes(x=relative_pos_rounded, y=mean_block_cM))+
  geom_line() 

