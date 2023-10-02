library(tidyverse)
library(ggstance)
library(RODBC)
library(ggcorrplot)
library(patchwork)

setwd("H:/")


FROH_id<-read.table("PhD_4th_yr/2023_ROH_search/2021_sleuthed_052023.hom.indiv", header=T, stringsAsFactors = F)%>%
dplyr::rename(Code=IID)%>%mutate(FROH=KB/2591865)%>%filter(nchar(Code)==5)

## pull out most and lest inbred indivs
tail_inbred=FROH_id%>%dplyr::arrange(FROH)%>%head(20)
tail_inbred_names=unique(tail_inbred$Code)

head_inbred=FROH_id%>%dplyr::arrange(desc(FROH))%>%head(20)
head_inbred_names=unique(head_inbred$Code)


top_bottom_nam=c(head_inbred_names,tail_inbred_names)

#basic stats 

sum(FROH_id$FROH>0)

# > 3195/3198
# [1] 0.9990619
# 


sum(FROH_id$FROH>=0.0625)

# > 1478/3198
# [1] 0.4621639


sum(FROH_id$FROH>=0.125)
# > 88/3198
# [1] 0.0275172


mean(FROH_id$FROH)
sd(FROH_id$FROH)





##histogram of all ids

ggplot(FROH_id, aes(FROH))+
  #geom_histogram(aes(y=..density..),bins = 60,colour="black", fill="white")+
  geom_density(alpha=0.6,fill="#FF6666")+
  theme_bw()+
  geom_vline(aes(xintercept=mean(FROH)))+
  geom_rug()
  # geom_density_ridges(jittered_points = TRUE, 
  #                     position = position_points_jitter(height = 0),
  #                     point_shape = '|', point_size = 3, 
  #                     point_alpha = 1, alpha = 0.7)
  # 

ggplot(FROH_id, aes(x=NSEG, y=KB))+
  geom_point()



### seperate out ROH into lengths 

FROH_full<-read.table("PhD_4th_yr/2023_ROH_search/2021_sleuthed_052023.hom", header=T, stringsAsFactors = F)



FROH_full$sizeKB <- NA
FROH_full$sizeKB <- ifelse(FROH_full$KB<=5000, "Short (2.5-5Mbp)", 
                     ifelse(FROH_full$KB>5000 & FROH_full$KB<=16000, "Medium (5-16Mbp)", 
                            ifelse(FROH_full$KB>16000, "Long (>16mbp)","NA"
                                  )))


table(FROH_full$sizeKB)

## stacked chromosome plot, number ROH per chr 

num_grouped_chr=FROH_full%>%
  dplyr::select(CHR,sizeKB)%>%
  group_by(CHR,sizeKB)%>%
  tally()

num_grouped_chr$CHR=as.factor(num_grouped_chr$CHR)
num_grouped_chr$CHR<-ordered(num_grouped_chr$CHR, levels = c("5", "18", "20", "9","11","12","19","15",
                                        "30","21","23","1","14","33","25","13","17","29","28", "4",
                                       "27","22", "24","8","3","31","6","7","2","16","32","10","26"))

num_grouped_chr_plot=ggplot(num_grouped_chr, aes(x=CHR, y=n, fill=sizeKB))+
  geom_bar(position = "stack", stat = "identity") +
  theme_classic()+
  scale_fill_manual(values = c("darksalmon","mediumaquamarine", "deepskyblue4"))+
  labs(x="Chromosome (Size ordered)", y="Number of ROH", fill="ROH size category")+
  theme(legend.position=c(0.9,0.9))
num_grouped_chr_plot

ggsave(num_grouped_chr_plot, 
       file="PhD_4th_yr/Chapter_1/Figs/ROH_chr_size_ordered_new_cols.png", 
       width = 10, 
       height = 8)





##############################################################################


# stacked ID plot 
num_grouped_id=FROH_full%>%
  dplyr::select(IID,sizeKB)%>%
  group_by(IID,sizeKB)%>%
  tally()

sample=sample(unique(num_grouped_id$IID), 150)

sub_num_grouped_id=num_grouped_id%>%filter(IID%in%sample)


ROH_size_ids_plot=ggplot(sub_num_grouped_id, aes(x=IID, y=n, fill=sizeKB))+
  geom_bar(position = "stack", stat = "identity") +
  theme_classic()+
  scale_fill_manual(values = c("darksalmon","mediumaquamarine", "deepskyblue4"))+
  labs(x="Individual", y="Number of ROH", fill="ROH size category")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        )

# 
# 
# ggsave(ROH_size_ids_plot, 
#        file="PhD_4th_yr/Chapter_1/Figs/ROH_ids_size_new_cols.png", 
#        width = 12, 
#        height = 6)
# 
# 

FROH_split=FROH_full%>%
  dplyr::select(IID, KB, sizeKB)%>%
  group_by(IID, sizeKB)%>%
  summarise(total_KB=sum(KB))%>%
  ungroup()%>%
  complete(IID,sizeKB, fill = list(total_KB=0))%>%
  dplyr::rename(Code=IID)%>%
  mutate(FROH_size=total_KB/2591865)

FROH_split$sizeKB=as.factor(FROH_split$sizeKB)

#boxplot of categories
sizes_boxplot=ggplot(FROH_split, aes(sizeKB, FROH_size, fill=sizeKB))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = c("darksalmon","mediumaquamarine", "deepskyblue4"))+
  coord_flip()+
  #geom_jitter(shape=16, position=position_jitter(0.2))+
  theme(legend.position="none",
        axis.title.y = element_blank()) +
  labs(y="FROH per size category")

sizes_boxplot

# ggsave(sizes_boxplot, 
#        file="PhD_4th_yr/Chapter_1/Figs/ROH_boxplot_new_cols.png", 
#        width = 7, 
#        height = 6)
# 


### FROH evolution ####

db<-"C:\\Users\\s1881212\\Documents\\Deer_database_2022/RedDeer2.05.1.accdb" #open connection
con<-odbcConnectAccess2007(db)
year<-sqlFetch(con, "tbllife")%>%dplyr::select(Code, BirthYear)

odbcClose(con)



FROH_evol=FROH_id%>%dplyr::select(Code, FROH)%>%
  inner_join(year)%>%
  mutate(year_cont=BirthYear-min(BirthYear)) # counting the years from first year of study to fit year as continuous variable 


table_birthyr=as.data.frame(table(FROH_evol$BirthYear))%>%mutate(Var1=as.numeric(as.character(Var1)))%>%
  filter(Var1>"1980")

min(table_birthyr$Freq)
max(table_birthyr$Freq)



FROH_evol_means=FROH_evol%>%dplyr::select(BirthYear, FROH)%>%
  filter(BirthYear>"1980")%>%
  group_by(BirthYear)%>%
  summarise(mean_FROH=mean(FROH))

FROH_evol$Code=as.factor(FROH_evol$Code)
FROH_evol$BirthYear=as.factor(FROH_evol$BirthYear)


ggplot(FROH_evol_means, aes(x=BirthYear, y=mean_FROH))+
  geom_point()+
  geom_smooth(method = lm, se=TRUE)+
  theme_bw()


## FROH split evolution



FROH_split_evol_means2=FROH_split%>%  
  inner_join(year)%>%
  dplyr::select(-Code, -total_KB)%>%
  filter(BirthYear>"1980")%>%
  group_by(BirthYear, sizeKB)%>%
  summarise(mean_FROH_size=mean(FROH_size))

FROH_evol_size=ggplot(FROH_split_evol_means2, aes(x=BirthYear, y=mean_FROH_size, colour=sizeKB))+
  facet_wrap(~sizeKB, scales = "free_y")+
  geom_point()+
  geom_smooth(method = lm, se=TRUE)+
  theme_bw()+
  scale_colour_manual(values = c("darksalmon","mediumaquamarine", "deepskyblue4"))+
  labs(y = expression(paste("Mean F"[ROH])))
 

# ggsave(FROH_evol_size, 
#        file="PhD_4th_yr/Chapter_1/Figs/FROH_evol_sizes.png", 
#        width = 12, 
#        height = 4)
# 
# 










## correlation of FROH with other inbreeding coefficients 

Fgrm_full<-read.table("PhD_4th_yr/test_Fgrm_for_loeske/Fgrm_full_subset.ibc", header=T, stringsAsFactors = F)%>%
  rename(Code=IID)


Fped_full=read.csv("PhD_4th_yr/Chapter_1/Pedigree_SEQ.csv", header=T, stringsAsFactors = F)%>%
  dplyr::select(ID, Inbreeding)%>%
  rename(Code=ID, Fped=Inbreeding)
  
  

## split FROH per size to check correlations etc.

FROH_only=FROH_id%>%dplyr::select(Code, FROH)
Fgrm_only=Fgrm_full%>%dplyr::select(Code, Fhat3)

ibc_corr=FROH_split%>%dplyr::select(-total_KB)%>%
  pivot_wider(names_from = sizeKB, values_from = FROH_size)%>%
  left_join(FROH_only)%>%
  left_join(Fgrm_only)%>%
  left_join(Fped_full)%>%
  dplyr::select(-Code)%>%
  na.omit()%>%
  rename("FROH Short"="Short (2.5-5Mbp)", "FROH Medium"="Medium (5-16Mbp)", "FROH Long" = "Long (>16mbp)", Fgrm=Fhat3)


ibc_corr_reorder=ibc_corr[,c(4,3,2,1,5,6)]


corr <- round(cor(ibc_corr_reorder), 1)
head(corr[, 1:6])

corr_matrix=ggcorrplot(corr, lab=TRUE)


ggsave(corr_matrix, 
       file="PhD_4th_yr/Chapter_1/Figs/corr_matrix.png", 
       width = 5, 
       height = 5)




vz_fgrm=ggplot(ibc_corr, aes(x=FROH, y=Fgrm))+
  geom_point(color="black",
             alpha=0.3,
             size=2,
             stroke = 1)+
  geom_smooth(method = lm, se=TRUE)+
  theme_bw()+
  labs(x="FROH", y="Fgrm")+
  ylim(-0.1, 0.4)+
  geom_rug(col="steelblue",alpha=0.1, size=1.5)



vs_fped=ggplot(ibc_corr, aes(x=FROH, y=Fped))+
  geom_point(  color="black",
               alpha=0.3,
               size=2,
               stroke = 1)+
  geom_smooth(method = lm, se=TRUE)+
  theme_bw()+
  labs(x="FROH", y="Fped")+
 ylim(-0.1, 0.4)+
  geom_rug(col="steelblue",alpha=0.1, size=1.5)



corr_plots=vz_fgrm+vs_fped


ggsave(corr_plots, 
       file="PhD_4th_yr/Chapter_1/Figs/corr_plots.png", 
       width = 10, 
       height = 5)


