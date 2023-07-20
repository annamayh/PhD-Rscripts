
#### summary of chr FROHs for data chapter 3 #####

library(tidyverse)
library(ggcorrplot)
library(ggridges)



setwd("H:/")

FROH_full<-read.table("PhD_4th_yr/2023_ROH_search/2021_sleuthed_052023.hom", header=T, stringsAsFactors = F)%>%
  dplyr::rename(Code=IID)%>%
  filter(nchar(Code)==5) #removing IDs with non-sensical ID codes


KB_per_chr_per_ID<-FROH_full%>%
  group_by(Code, CHR)%>%
  summarise(KB_chr=sum(KB))%>% ##for every id and every chr summing the KB 
  ungroup %>% 
  complete(Code, CHR, fill = list(KB_chr = 0)) #completes for all chromosomes 


deermap <- read.csv("PhD_3rdYR/Data_files/Genome_assembly_mCerEla1.1.csv", header = T, stringsAsFactors = F)%>%
  filter(!CHR %in% c("All","All_auto","X","unplaced"))%>%## reading in chromosome sizes from genome assembly 
  mutate(CHR=as.numeric(CHR))
  
froh_per_chr<-full_join(KB_per_chr_per_ID,deermap, by = join_by(CHR))%>% 
  mutate(chr_froh=KB_chr/length_Kb)%>%
  select(Code, CHR, chr_froh)%>% 
  reshape2::dcast(Code~CHR) %>% 
  mutate(FROHsum = rowSums(.[2:34]))%>%
  mutate(FROH_sum_div=FROHsum/33)


colnames(froh_per_chr) <- c("Code", paste0("FROH_chr", 1:33),"FROHsum","FROH_sum_div")



# write.table(froh_per_chr,
#             file = "PhD_4th_yr/Inbreeding_depression_models/FROH_per_Chr.txt",
#             row.names = F, quote = F, sep = ",",na = "NA") #saving tables as txt file 
# 
# 
# 



corr=froh_per_chr%>%
  select(-Code,-FROHsum,-FROH_sum_div)%>%cor()

cor(froh_per_chr$FROH_chr8, froh_per_chr$FROH_chr22)

max(corr[corr<1])
min(corr)
# > max(corr[corr<1])
# [1] 0.1354674
# > min(corr)
# [1] -0.01329034

corr_matrix=ggcorrplot(corr)+
  #scale_fill_gradient2(low = "white", high = "red", breaks=c(0, 1), limit=c(-0.1, 0.15))
  scale_fill_gradientn(colours = c("blue", "white", "seagreen"),limit=c(-0.12, 0.14))
  

corr_matrix

ggsave(corr_matrix,
       file = "PhD_4th_yr/Inbreeding_depression_chapter/plots/corr_matrix_chrFROH.png",
       width = 10,
       height = 7, 
       bg = "white")




covm=froh_per_chr%>%
  select(-Code,-FROHsum,-FROH_sum_div)%>%cov()

max(covm)
min(covm)




  
froh_per_chr_2<-full_join(KB_per_chr_per_ID,deermap, by = join_by(CHR))%>% 
    mutate(chr_froh=KB_chr/length_Kb)%>%
    select(Code, CHR, chr_froh)

mean(froh_per_chr_2$chr_froh)
max(froh_per_chr_2$chr_froh) ##0.94 on chr 3 CUC84 and ZAP13


boxplot_chrfroh=ggplot(froh_per_chr_2, aes(x = chr_froh, y = as.factor(CHR), fill = as.factor(CHR))) +
  geom_boxplot(alpha=0.5) +
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 18))+
  labs(x=expression(F["ROHchr"]), y="Chromosome")

boxplot_chrfroh

ggsave(boxplot_chrfroh,
       file = "PhD_4th_yr/Inbreeding_depression_chapter/plots/boxplot_chrFROH.png",
       width = 6,
       height = 9)



