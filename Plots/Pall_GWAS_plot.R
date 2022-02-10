library(plyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(forcats)

setwd("M:/")

CMpall<-read.table("PHD_2ndYR/Recombination_ROH/Data_files/Pall_CM_filt_moreSNPs.txt", header = TRUE,stringsAsFactors = FALSE)
CMpall$CHR<-as.factor(CMpall$CHR)
CMpall$CHR<-ordered(CMpall$CHR, levels = c("5", "18", "20", "9","11","12","19","15",
                                          "30","21","23","1","14","33","25","13","17","29","28", "4",
                                          "27","22", "24","8","3","31","6","7","2","16","32","10","26"))


CMpall<-CMpall%>%arrange(CHR)%>%mutate(perc=pall*100)
CMpall$Order <- 1:nrow(CMpall)


axis.set <- CMpall %>% 
  group_by(CHR) %>% 
  summarize(center = (max(Order) + min(Order)) / 2)

CMpall$CHR<-as.numeric(CMpall$CHR)

CM_g<-ggplot(CMpall, aes(Order, perc, col = as.factor(CHR%% 2))) +
  scale_colour_manual(values = c("coral", "cornflowerblue")) +
  geom_point() +
  theme_classic()+
  labs(x="Chromosome (Size ordered)",y="ROH density \n(% of ids with a ROH)", title="Using genetic map positions")+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 8, vjust = 0.5),
        axis.text.y = element_text(size = 12, vjust = 0.5)) +
  geom_hline(yintercept = 14.9, linetype="dashed", color = "red") +
  geom_hline(yintercept = 6.5, colour = "red")+
  scale_x_continuous(limits = c(0,25798), expand = c(0, 0), label = axis.set$CHR, breaks = axis.set$center)+
  scale_y_continuous(limits = c(0,20))



###### BP genomwide plot using same script ################

BPpall<-read.table("PhD_3rdYR/Data_files/Rum_Pall_BP_filt_UPDATED_SANGER_POS.txt", header = TRUE,stringsAsFactors = FALSE)
BPpall$CHR<-as.factor(BPpall$CHR)
BPpall$CHR<-ordered(BPpall$CHR, levels = c("5", "18", "20", "9","11","12","19","15",
                                           "30","21","23","1","14","33","25","13","17","29","28", "4",
                                           "27","22", "24","8","3","31","6","7","2","16","32","10","26"))


BPpall<-BPpall%>%arrange(CHR)%>%mutate(perc=pall*100)
BPpall$Order <- 1:nrow(BPpall)


axis.set_BP <- BPpall %>% 
  group_by(CHR) %>% 
  summarize(center = (max(Order) + min(Order)) / 2)

BPpall$CHR<-as.numeric(BPpall$CHR)

BP_g<-ggplot(BPpall, aes(Order, perc, col = as.factor(CHR%% 2))) +
  scale_colour_manual(values = c("coral", "cornflowerblue")) +
  geom_point() +
  theme_classic()+
  labs(tag = "A", y="ROH density \n(%of ids with a ROH)", title = "Using physical map positions")+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 8, vjust = 0.5),
        axis.text.y = element_text(size = 12, vjust = 0.5)) +
  geom_hline(yintercept = 15.9, linetype="dashed", color = "red") +
  geom_hline(yintercept = 6.23, colour = "red")+
  scale_x_continuous(limits = c(0,28875), expand = c(0, 0), label = axis.set_BP$CHR, breaks = axis.set_BP$center)+
  scale_y_continuous(limits = c(0,20))

quantile(BPpall$pall, c(.01, .99)) ## 
quantile(CMpall$pall, c(.01, .99)) ## 
mean(CMpall$pall)
mean(BPpall$pall)


BP_g +CM_g + plot_layout(ncol=1)



