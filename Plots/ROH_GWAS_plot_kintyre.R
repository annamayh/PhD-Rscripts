## GWAS like plot of ROH density for Kintyre pop 

library(plyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(forcats)

setwd("H:/PHD_2ndYR")
CM_kint<-read.table("Recombination ROH/Data_files/KINTYRE_Pall_CM_filt_moreSNPs.txt", header = TRUE,stringsAsFactors = FALSE)
CM_kint$CHR<-as.factor(CM_kint$CHR)
CM_kint$CHR<-ordered(CM_kint$CHR, levels = c("5", "18", "20", "9","11","12","19","15",
                                           "30","21","23","1","14","33","25","13","17","29","28", "4",
                                           "27","22", "24","8","3","31","6","7","2","16","32","10","26"))


CM_kint<-CM_kint%>%arrange(CHR)
CM_kint$Order <- 1:nrow(CM_kint)


axis.set <- CM_kint %>% 
  group_by(CHR) %>% 
  summarize(center = (max(Order) + min(Order)) / 2)

CM_kint$CHR<-as.numeric(CM_kint$CHR)

CM_plot_kint<-ggplot(CM_kint, aes(Order, pall, col = as.factor(CHR%% 2))) +
  scale_colour_manual(values = c("coral", "cornflowerblue")) +
  geom_point() +
  theme_classic()+
  labs(tag = "Kintyre", y="ROH density", title="Using genetic (cM) map")+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 8, vjust = 0.5),
        axis.text.y = element_text(size = 12, vjust = 0.5)) +
  geom_hline(yintercept = 0.083, linetype="dashed", color = "red") +
  geom_hline(yintercept = 0.034, colour = "red")+
  scale_x_continuous(limits = c(0,25200), expand = c(0, 0), label = axis.set$CHR, breaks = axis.set$center)



###### BP genomwide plot using same script ################

BP_kint<-read.table("Recombination ROH/Data_files/Kintyre_Pall_BP_filt.txt", header = TRUE,stringsAsFactors = FALSE)
BP_kint$CHR<-as.factor(BP_kint$CHR)
BP_kint$CHR<-ordered(BP_kint$CHR, levels = c("5", "18", "20", "9","11","12","19","15",
                                           "30","21","23","1","14","33","25","13","17","29","28", "4",
                                           "27","22", "24","8","3","31","6","7","2","16","32","10","26"))


BP_kint<-BP_kint%>%arrange(CHR)
BP_kint$Order <- 1:nrow(BP_kint)


axis.set_BP <- BP_kint %>% 
  group_by(CHR) %>% 
  summarize(center = (max(Order) + min(Order)) / 2)

BP_kint$CHR<-as.numeric(BP_kint$CHR)

BP_plot_kint<-ggplot(BP_kint, aes(Order, pall, col = as.factor(CHR%% 2))) +
  scale_colour_manual(values = c("coral", "cornflowerblue")) +
  geom_point() +
  theme_classic()+
  labs(x="Linkage group (Size ordered)", y="ROH density", title = "Using physical (bp) map")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8, vjust = 0.5),
        axis.text.y = element_text(size = 12, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0.0955, linetype="dashed", color = "red") +
  geom_hline(yintercept = 0.032, colour = "red")+
  scale_x_continuous(limits = c(0,25700), expand = c(0, 0), label = axis.set_BP$CHR, breaks = axis.set_BP$center)



CM_plot_kint + BP_plot_kint + plot_layout(ncol=1)


# save plot as W=980 H=440

#idk why its not working 
ggsave("Manuscript_draft/images/GWAS_ROH_Kintyre.png", plot = GWAS_plot, 
       scale = 1, 
       width = 970,
       height = 440, 
       units = "px"
       )






