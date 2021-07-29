library(plyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)

setwd("H:/PHD_2ndYR/Recombination ROH/Data_files")
CMpall_main <- read.table("KINTYRE_Pall_CM_filt_moreSNPs.txt", header = TRUE)#reading in cM filtered SNPs kintyre data with ROH density values attched (pall)
CMpall<-read.table("Pall_CM_filt_moreSNPs.txt", header = TRUE)#filtered SNPs for Rum pop using cM map, with ROH denisty values

#getting means and hotspot/coldspot thresholds for all
se <- function(x) sqrt(var(x)/length(x))

#stats for Rum pop using genetic map
mean(CMpall$pall) 
se(CMpall$pall)
quantile(CMpall$pall, c(.01, .99)) ## 
length(which(CMpall$pall > 0.1486795)) ##
length(which(CMpall$pall <= 0)) ## 
hotspots_CM_Rum<-subset(CMpall, pall > 0.1486795) %>%rename(pall_Rum=pall)%>%select(pall_Rum, SNP, CHR)
unique(hotspots_CM_Rum$CHR)
coldspots_CM_Rum<-subset(CMpall, pall <=0) 
unique(coldspots_CM_Rum$CHR)

#stats for Kint pop using genetic map
mean(CMpall_main$pall) 
quantile(CMpall_main$pall, c(.01, .99)) ## 
length(which(CMpall_main$pall > 0.08280255)) ##
length(which(CMpall_main$pall <= 0)) ## 
hotspots_CM_kint<-subset(CMpall_main, pall > 0.08280255)%>%rename(pall_kint=pall)%>%select(pall_kint, SNP, CHR)
unique(hotspots_CM_kint$CHR)
coldspots_CM_kint<-subset(CMpall_main, pall <= 0) 
unique(coldspots_CM_kint$CHR)


## shared hotspots between two pops using cM map

CM_hot_shared_pops<-join(hotspots_CM_kint,hotspots_CM_Rum)#%>%na.omit()
#### figures to compare maps using cM ####

main_sub <- subset(CMpall_main, CHR == "18")
main_sub$order<-1:nrow(main_sub)
mainland_CM<-ggplot(main_sub, aes(order, pall)) +
  geom_point() +
  theme_classic()+
  labs(x = "SNP order", y = "P") +
  geom_hline(yintercept = 0.03199692,  color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0.08280255, linetype="dashed",colour = "red")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

rum_sub <- subset(CMpall, CHR == "18")
rum_sub$order<-1:nrow(rum_sub)
rum_CM <- ggplot(rum_sub, aes(order, pall)) +
  geom_point() +
  theme_classic()+
  labs(x = "SNP order", y = "P") +
  geom_hline(yintercept = 0.1486795, linetype="dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0.06085681, colour = "red")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



##### correlation between cm pall for Rum and kintyre ###
CM_rum<-CMpall%>%select(CHR,SNP,pall)%>%dplyr::rename(pall_Rum = pall)
CM_kint<-CMpall_main%>%select(CHR,SNP,pall)%>%dplyr::rename(pall_kint = pall)
CM_corr<-join(CM_rum,CM_kint)
cor.test(CM_corr$pall_Rum, CM_corr$pall_kint, method="pearson")

#plot(CM_corr$pall_Rum, CM_corr$pall_kint)
#abline(lm(CM_corr$pall_kint ~ CM_corr$pall_Rum, col = "blue"))


#############################################
##### BP maps###########################
BP_palls<-read.table("Rum_Pall_BP_filt.txt", header = TRUE)
BP_pall_main <- read.table("Kintyre_Pall_BP_filt.txt", header = TRUE)

# Rum stats#
mean(BP_palls$pall)
se(BP_palls$pall)
quantile(BP_palls$pall, c(.01, .99)) ## 
length(which(BP_palls$pall > 0.146667754)) ##
length(which(BP_palls$pall < 0.000652103)) ## 
hotspots_BP_Rum<-subset(BP_palls, pall > 0.146667754)%>%rename(pall_Rum_BP=pall)%>%select(pall_Rum_BP, SNP, CHR)
unique(hotspots_BP_Rum$CHR)
coldspots_BP_Rum<-subset(BP_palls, pall <= 0) 
unique(coldspots_BP_Rum$CHR)

#Kintyre stats#
mean(BP_pall_main$pall)
se(BP_pall_main$pall)
quantile(BP_pall_main$pall, c(.01, .99)) ## 
length(which(BP_pall_main$pall > 0.0955414)) ##
length(which(BP_pall_main$pall <= 0)) ## 
hotspots_BP_kint<-subset(BP_pall_main, pall > 0.0955414)%>%rename(pall_kint_BP=pall)%>%select(pall_kint_BP, SNP, CHR)
unique(hotspots_BP_kint$CHR)
coldspots_BP_kint<-subset(BP_pall_main, pall <= 0) 
unique(coldspots_BP_kint$CHR)



Rum_hotspots_shared<- join(hotspots_CM_Rum,hotspots_BP_Rum)%>%na.omit()

BP_hotspots_shared<-join(hotspots_BP_Rum,hotspots_BP_kint)%>%na.omit()

###################################################################################################################


main_sub_BP <- subset(BP_pall_main, CHR == "18")
mainland_BP<-ggplot(main_sub_BP, aes(BP, pall)) +
  geom_point() +
  theme_classic()+
  labs(x = "SNP order", y = "P") +
  geom_hline(yintercept = 0.02929397,  color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0.0955414, linetype="dashed",colour = "red") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

rum_sub_BP <- subset(BP_palls, CHR == "18")
rum_BP <- ggplot(rum_sub_BP, aes(BP, pall)) +
  geom_point() +
  theme_classic()+
  labs(x = "SNP order", y = "P") +
  geom_hline(yintercept = 0.146667754, linetype="dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0.05770727, colour = "red")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  


BP_rum<-BP_palls%>%select(CHR,SNP,pall)%>%dplyr::rename(pall_Rum = pall)
BP_kint<-BP_pall_main%>%select(CHR,SNP,pall)%>%dplyr::rename(pall_kint = pall)

BP_corr<-join(BP_rum,BP_kint)

cor.test(BP_corr$pall_Rum, BP_corr$pall_kint, method="pearson")

plot(BP_corr$pall_Rum, BP_corr$pall_kint)
abline(lm(BP_corr$pall_kint ~ BP_corr$pall_Rum, col = "blue"))

#############################
#### using patchwork ####


row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Rum population", angle = 90) + theme_void() 
row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Kintyre population", angle = 90) + theme_void() 
col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Using cM map") + theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Using BP map") + theme_void() 
title <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Distrbution of ROH across LG 18 using different methods and populations",size=5) + theme_void() 


layoutplot <- "
#iiiii#
#cccddd
aeeefff
aeeefff
bggghhh
bggghhh
"## layout of plot in rectangle using above labels specified ... # puts space

plotlist <- list(i=title,a = row1, b = row2, c = col1, d = col2, e= rum_CM, f=rum_BP, g=mainland_CM, h=mainland_BP)

wrap_plots(plotlist, guides = 'collect', design = layoutplot, title(main="Q"))




