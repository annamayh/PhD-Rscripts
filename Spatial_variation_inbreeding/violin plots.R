
library(tidyverse)
library(patchwork)


setwd("H:/")

surv_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/survival_loc.txt", sep = ",", header = TRUE)%>%
  select(-BirthWt)%>%na.omit()

table(surv_loc_df$Reg)
# IM  LA  MG  NG  SG  SI 
# 382 362 303 505 260 428 


calf_surving=surv_loc_df%>%
  filter(juvenile_survival==1)


mum=surv_loc_df%>%
  distinct(Reg,MumCode, .keep_all = TRUE)%>%
  ggplot(aes(x=MumFROH, y=Reg)) +
  geom_violin()+
  geom_point(position = "jitter", alpha=0.2) +
  theme_classic()+
  xlim(0, 0.35)
  

calf=surv_loc_df%>%
  ggplot(aes(x=FROH, y=Reg, fill=Reg)) +
  geom_violin(alpha=0.8)+
  geom_point(position = "jitter", alpha=0.2) +
  theme_classic()+
  scale_fill_manual(values = c("#a20025","#00aba9","#60a917","chocolate1", "#647687","#f0a30a" ))+
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size=12),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size=10))+
  xlim(0, 0.35)+
  geom_vline(xintercept=0.1, linetype="dashed", alpha=0.5)+
  geom_vline(xintercept=0.2, linetype="dashed",alpha=0.5)+
  geom_vline(xintercept=0.3, linetype="dashed",alpha=0.5)+
  labs(x=expression(F["ROH"]), y="Spatial region", title = "All calves")
  

calf_surv=surv_loc_df%>%
  filter(juvenile_survival==1)%>%
  ggplot(aes(x=FROH, y=Reg, fill=Reg)) +
  geom_violin(alpha=0.8)+
  geom_point(position = "jitter", alpha=0.2) +
  theme_classic()+
  scale_fill_manual(values = c("#a20025","#00aba9","#60a917","chocolate1", "#647687","#f0a30a" ))+
  theme(axis.title.y=element_blank(),
         axis.text.y = element_text(size=12),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size=10))+
  xlim(0, 0.35)+
  geom_vline(xintercept=0.1, linetype="dashed",alpha=0.5)+
  geom_vline(xintercept=0.2, linetype="dashed",alpha=0.5)+
  geom_vline(xintercept=0.3, linetype="dashed",alpha=0.5)+
  labs(x=expression(F["ROH"]), title="Calves surviving juvenile stage")


calf_vio=calf+calf_surv
calf_vio

ggsave(calf_vio,
       file = "PhD_4th_yr/Spatial_var_inbreeding/juvenile_surv_violin.png",
       width = 10,
       height = 7)



#################################################################################################
### winter survival ############################################################################

winter_surv_loc_df=read.table("PhD_4th_yr/Spatial_var_inbreeding/winter_survival_loc.txt", sep = ",", header = TRUE)%>%
 na.omit()


calf_surving_winter=winter_surv_loc_df%>%
  filter(winter_survival==1)

calf_wint=calf_surving_winter%>%
  ggplot(aes(x=FROH, y=Reg, fill=Reg)) +
  geom_violin(alpha=0.8)+
  geom_point(position = "jitter", alpha=0.2) +
  theme_classic()+
  scale_fill_manual(values = c("#a20025","#00aba9","#60a917","chocolate1", "#647687","#f0a30a" ))+
  theme(legend.position = "none",
        axis.title.y = element_blank(), axis.text.y = element_text(size=12),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size=10))+
  xlim(0, 0.35)+
  geom_vline(xintercept=0.1, linetype="dashed", alpha=0.5)+
  geom_vline(xintercept=0.2, linetype="dashed",alpha=0.5)+
  geom_vline(xintercept=0.3, linetype="dashed",alpha=0.5)+
  labs(x=expression(F["ROH"]), y="Spatial region", title = "Calves surviving winter")




calf_vi=calf+calf_wint+calf_surv
calf_vi

ggsave(calf_vi,
       file = "PhD_4th_yr/Spatial_var_inbreeding/juv_wint_surv_violin.png",
       width = 12,
       height = 7)

