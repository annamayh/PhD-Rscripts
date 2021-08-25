library(vcfR)
library(rehh)
library(data.table)
library(R.utils)
library(dplyr)
library(plyr)
library(purrr)
library(ggplot2)
library(tidyverse)
library(janitor)
library(stringdist)

### working out iHS for emperical dataset 

setwd("H:/PHD_2ndYR")

## getting correct map input file ####
marker_all<-NULL

for (k in 1:33){
  markers_oneLG<-read.table(paste0("Haplotype_diversity/AlphaPeel_files/all_chr/Deer_AlphaPeel_marker_file_chr", k, ".txt"), stringsAsFactors = FALSE, header=TRUE)%>%
    dplyr::rename(SNP.Name=SNP.name)
  marker_all <- rbind(marker_all, markers_oneLG)
}

filtered_maps<-join(map_file,marker_all)%>%select(-marker.order,-cMPosition.SexAveraged)%>%select(SNP.Name,CEL.LG,Est.Mb.pos,A1,A2)%>%na.omit()

#write.table(filtered_maps,
 #           file = "Deer_data/map_positions_SNPalphapeel.txt",
  #          row.names = F, quote = F, sep = "\t")


## now getting correct format for haplotype files ####

rum_roh<-read.table("Haplotype_diversity/old_files/Froh_2019.txt", header = TRUE) %>%
  select(IID) ## IDs used in ROH analysis 
rum_roh$in_dataset<-"yes" #adding vector to filter by later
phased_lab<-plyr::join(phased, rum_roh)%>%filter(in_dataset == "yes")%>% 
  select(-IID, -in_dataset) #filtering out IDs not used in ROH analysis 

for (l in 1:33){
phased_LG<-read.table(paste0("Haplotype_diversity/AlphaPeel_files/all_chr/AlphaPeel_out_chr", l, ".called_phase.0.95"), stringsAsFactors = FALSE)%>%
  dplyr::rename(IID=V1)

phased_lab<-plyr::join(phased_LG, rum_roh)%>%filter(in_dataset == "yes")%>% 
  select(-in_dataset) #filtering out IDs not used in ROH analysis 

phased_lab[phased_lab==9]<-"."

#write.table(phased_lab,
 #           file = paste0("Haplotype_diversity/Phased_formatted_for_Rehh/Phased_formatted_LG",l,".txt"),
  #          row.names = T, col.names= T, quote = F, sep = "\t")
}







#########################################################################################################################
###################### CALCULATING iHS #################################################################################
########################################################################################################################



### another approach 
output_haplohh<-list()

for (i in 1:33){
  
  hap_file_input = paste0("Haplotype_diversity/Phased_formatted_for_Rehh/Phased_formatted_LG", i,".txt")
  # create internal representation
  output_haplohh[[i]] <- data2haplohh(hap_file = hap_file_input,
                                      map_file = paste0("Deer_data/map_positions_SNPalphapeel.txt"),
                                      chr.name = i,
                                      remove_multiple_markers = TRUE,
                                      allele_coding = "none", 
                                      min_perc_geno.mrk = 90) # creating list of all LG containing output from haplohh function
  
}

scan_output_wg<-map(output_haplohh, scan_hh)%>%unlist()
wgscan.ihs <- ihh2ihs(scan_output_wg)  


       
##########################################################################################################################       
## demo code taken from rehh ############################################################################################
#########################################################################################################################
for(i in 1:33) {
  # haplotype file name for each chromosome
  hap_file_input = paste0("Haplotype_diversity/Phased_formatted_for_Rehh/Phased_formatted_LG", i,".txt")
  # create internal representation
  hh <- data2haplohh(hap_file = hap_file_input,
                     map_file = paste0("Deer_data/map_positions_SNPalphapeel.txt"),
                     chr.name = i,
                     remove_multiple_markers = TRUE,
                     allele_coding = "none", 
                     min_perc_geno.mrk = 90)
  # perform scan on a single chromosome (calculate iHH values)
  scan <- scan_hh(hh)
  # concatenate chromosome-wise data frames to
  # a data frame for the whole genome
  # (more efficient ways certainly exist...)
  if (i == 1) {
    wgscan <- scan
  } else {
    wgscan <- rbind(wgscan, scan)
  }
}
# calculate genome-wide iHS values


#save(hh, scan, wgscan, file="iHS_Rum.Rdata")
# saved output of loop so can just load it in 

load("Haplotype_diversity/iHS_Rum.Rdata")

wgscan.ihs <- ihh2ihs(wgscan)       
       
       
iHS_genomewide<-as.data.frame(wgscan.ihs$ihs)


ihs_SNPs<-tibble::rownames_to_column(iHS_genomewide, "SNP")


## read in ROH values 

ROH<-read.table("Recombination ROH/CM_ROH_output/ROH_output_trueCM.hom.summary", header=T)%>%mutate(pall=UNAFF/3046)


ROH_vs_iHS<-join(ihs_SNPs,ROH)%>%na.omit()

Top_1<-as.numeric(quantile(ROH_vs_iHS$pall, c(.99)))
ROH_vs_iHS$hotspot<-NA
ROH_vs_iHS$hotspot<-ifelse(ROH_vs_iHS$pall>Top_1, "yes", "no")


#######################

iHS_hotspots<-ROH_vs_iHS%>%filter(hotspot=="yes")
  
  
(quantile(ROH_vs_iHS$IHS, c(.01,.99)))

ggplot(ROH_vs_iHS, aes(hotspot, IHS, fill=hotspot)) + geom_boxplot()+
  theme_classic()+
  labs(title = "Haplotype diveristy using iHS")


ggplot(ROH_vs_iHS, aes(hotspot, LOGPVALUE, fill=hotspot)) + geom_boxplot()+
  theme_classic()+
  labs(title = "Haplotype diveristy using iHS")


ROH_vs_iHS$Order<-1:nrow(ROH_vs_iHS)

ROH_vs_iHS$CHR<-as.factor(ROH_vs_iHS$CHR)
ggplot(ROH_vs_iHS, aes(Order, IHS, col = (CHR%%2))) +
  scale_colour_manual(values = c("coral", "cornflowerblue")) +
  geom_point() +
  theme_classic()


axis.set <- ROH_vs_iHS %>% 
  group_by(CHR) %>% 
  summarize(center = (max(Order) + min(Order)) / 2)


ROH_vs_iHS$CHR<-as.numeric(ROH_vs_iHS$CHR)

ggplot(ROH_vs_iHS, aes(Order, LOGPVALUE, col = as.factor(CHR%%2))) +
  scale_colour_manual(values = c("black", "gray")) +
  geom_point() +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(0,36000), expand = c(0, 0), label = axis.set$CHR, breaks = axis.set$center)+
  geom_point(data=iHS_hotspots,
             aes(x=Order,y=LOGPVALUE), color="red",size=3)

manhattanplot(wgscan.ihs,
              main = "iHS")

manhattanplot(wgscan.ihs,
              pval = TRUE,
              threshold = 4,
              main = "p-value of iHS (CGU cattle breed)")
