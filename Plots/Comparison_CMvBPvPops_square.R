library(plyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)

setwd("M:/PHD_2ndYR/Recombination_ROH/Data_files")
CMpall_main <- read.table("KINTYRE_Pall_CM_filt_moreSNPs.txt", header = TRUE)
CMpall<-read.table("Pall_CM_filt_moreSNPs.txt", header = TRUE)

mean(CMpall_main$pall) 
quantile(CMpall_main$pall, c(.01, .99)) ## 
length(which(CMpall_main$pall > 0.08280255)) ##
hotspots_kint_gen<-subset(CMpall_main, pall>0.08280255)%>%mutate(pall_cm_kint=pall)%>%select(CHR,pall_cm_kint,SNP)


mean(CMpall$pall) 
quantile(CMpall$pall, c(.01, .99)) ## 
length(which(CMpall$pall > 0.1486795 )) ##
hotspots_rum_gen<-subset(CMpall, pall>0.1486795)%>%mutate(pall_cm_rum=pall)%>%select(CHR,pall_cm_rum,SNP)



#### CM maps ####

main_sub <- subset(CMpall_main, CHR == "18")
main_sub$order<-1:nrow(main_sub)
mainland_CM<-ggplot(main_sub, aes(order, pall)) +
  geom_point() +
  theme_classic()+
  labs(x = "SNP order", y = "ROH density") +
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
  labs(x = "SNP order", y = "ROH denisty") +
  geom_hline(yintercept = 0.1471829, linetype="dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0.06085681, colour = "red")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


#############################################
##### BP maps###########################


BP_pall_main <- read.table("Kintyre_BP_UpdatedMbSANGER_filtSNPS.txt", header = TRUE)
BP_palls<-read.table("Rum_Pall_BP_filt_UPDATED_SANGER_POS.txt", header = TRUE)

mean(BP_pall_main$pall) 
quantile(BP_pall_main$pall, c(.01, .99)) ## 
length(which(BP_pall_main$pall > 0.08917197  )) ##
hotspots_kint<-subset(BP_pall_main, pall>0.08917197)%>%mutate(pall_bp_kint=pall)%>%select(CHR,pall_bp_kint,SNP)

mean(BP_palls$pall) 
quantile(BP_palls$pall, c(.01, .99)) ## 
length(which(BP_palls$pall > 0.159197913  )) ##
length(which(BP_palls$pall <= 0)) ##  
hotspots_rum<-subset(BP_palls, pall>0.159197913)%>%mutate(pall_bp_rum=pall)%>%select(CHR,pall_bp_rum,SNP)
unique(hotspots_rum[c("CHR")])

main_sub_BP <- subset(BP_pall_main, CHR == "18")
mainland_BP<-ggplot(main_sub_BP, aes(BP, pall)) +
  geom_point() +
  theme_classic()+
  labs(x = "SNP order", y = "ROH denisty") +
  geom_hline(yintercept = 0.0386007, linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0.08917197, linetype="dashed",colour = "red") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

rum_sub_BP <- subset(BP_palls, CHR == "18")
rum_BP <- ggplot(rum_sub_BP, aes(BP, pall)) +
  geom_point() +
  theme_classic()+
  labs(x = "SNP order", y = "ROH denisty") +
  geom_hline(yintercept = 0.159197913, linetype="dashed", color = "red") +
  geom_hline(yintercept = 0.07320247, colour = "red")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  #annotate("rect", xmin = 3, xmax = 4.2, ymin = 12, ymax = 21,
   #        alpha = .2)
  


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



rum<-join(hotspots_rum, hotspots_rum_gen)%>%na.omit()
kint<-join(hotspots_kint,hotspots_kint_gen)%>%na.omit()
cm<-join(hotspots_kint_gen,hotspots_rum_gen)%>%na.omit()
bp<-join(hotspots_kint,hotspots_rum)%>%na.omit()
