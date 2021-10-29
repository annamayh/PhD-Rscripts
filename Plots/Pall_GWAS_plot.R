library(plyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(forcats)

setwd("H:/PHD_2ndYR")

CMpall<-read.table("Recombination ROH/Data_files/Pall_CM_filt_moreSNPs.txt", header = TRUE,stringsAsFactors = FALSE)
CMpall$CHR<-as.factor(CMpall$CHR)
CMpall$CHR<-ordered(CMpall$CHR, levels = c("5", "18", "20", "9","11","12","19","15",
                                          "30","21","23","1","14","33","25","13","17","29","28", "4",
                                          "27","22", "24","8","3","31","6","7","2","16","32","10","26"))


CMpall<-CMpall%>%arrange(CHR)
CMpall$Order <- 1:nrow(CMpall)


axis.set <- CMpall %>% 
  group_by(CHR) %>% 
  summarize(center = (max(Order) + min(Order)) / 2)

CMpall$CHR<-as.numeric(CMpall$CHR)

CM_g<-ggplot(CMpall, aes(Order, pall, col = as.factor(CHR%% 2))) +
  scale_colour_manual(values = c("coral", "cornflowerblue")) +
  geom_point() +
  theme_classic()+
  labs(tag = "Rum", x="Linkage group (Size ordered)",y="ROH density", title="Using genetic (cM) map")+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 8, vjust = 0.5),
        axis.text.y = element_text(size = 12, vjust = 0.5)) +
  geom_hline(yintercept = 0.1486795, linetype="dashed", color = "red") +
  geom_hline(yintercept = 0.06409285, colour = "red")+
  scale_x_continuous(limits = c(0,25798), expand = c(0, 0), label = axis.set$CHR, breaks = axis.set$center)



###### BP genomwide plot using same script ################

BPpall<-read.table("Recombination ROH/Data_files/Rum_Pall_BP_filt.txt", header = TRUE,stringsAsFactors = FALSE)
BPpall$CHR<-as.factor(BPpall$CHR)
BPpall$CHR<-ordered(BPpall$CHR, levels = c("5", "18", "20", "9","11","12","19","15",
                                           "30","21","23","1","14","33","25","13","17","29","28", "4",
                                           "27","22", "24","8","3","31","6","7","2","16","32","10","26"))


BPpall<-BPpall%>%arrange(CHR)
BPpall$Order <- 1:nrow(BPpall)


axis.set_BP <- BPpall %>% 
  group_by(CHR) %>% 
  summarize(center = (max(Order) + min(Order)) / 2)

BPpall$CHR<-as.numeric(BPpall$CHR)

BP_g<-ggplot(BPpall, aes(Order, pall, col = as.factor(CHR%% 2))) +
  scale_colour_manual(values = c("coral", "cornflowerblue")) +
  geom_point() +
  theme_classic()+
  labs(x="Linkage group (Size ordered)", y="ROH density", title = "Using physical (BP) map")+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 8, vjust = 0.5),
        axis.text.y = element_text(size = 12, vjust = 0.5)) +
  geom_hline(yintercept = 0.146667754, linetype="dashed", color = "red") +
  geom_hline(yintercept = 0.06237447, colour = "red")+
  scale_x_continuous(limits = c(0,26318), expand = c(0, 0), label = axis.set_BP$CHR, breaks = axis.set_BP$center)

quantile(CMpall$pall, c(.01, .99)) ## 
mean(CMpall$pall)
mean(BPpall$pall)


CM_g + BP_g + plot_layout(ncol=1)



