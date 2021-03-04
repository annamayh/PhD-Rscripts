library(plyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)

setwd("H:/PHD_2ndYR/Recombination ROH")
CMpall_main <- read.table("KINTYRE_Pall_CM_filt_moreSNPs.txt", header = TRUE)
CMpall<-read.table("Pall_CM_filt_moreSNPs.txt", header = TRUE)

mean(CMpall$pall) 
quantile(CMpall$pall, c(.01, .99)) ## 
length(which(CMpall_main$pall > 0.08280255)) ##
length(which(CMpall_main$pall <= 0)) ## 


#### CM maps ####

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
  geom_hline(yintercept = 0.1471829, linetype="dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0.06085681, colour = "red")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


#############################################
##### BP maps###########################

BP_pall_main <- read.table("Kintyre_Pall_BP_filt.txt", header = TRUE)
BP_palls<-read.table("Rum_Pall_BP_filt.txt", header = TRUE)


mean(BP_palls$pall) 
quantile(BP_palls$pall, c(.01, .99)) ## 


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
  geom_hline(yintercept = 0.1436224, linetype="dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0.05770727, colour = "red")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  


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




