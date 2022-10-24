### For all 33 chr finding haplotype diveristy for 10 SNP windows ###


library(dplyr)
library(tidyverse)
library(janitor)
library(stringdist)

setwd("H:/PHD_2ndYR/Haplotype_diversity")


HRiD=list()


for (k in 1:33){
  
  phased<-read.table(paste0("AlphaPeel_files/all_chr/AlphaPeel_out_chr", k, ".called_phase.0.95"),stringsAsFactors = FALSE)
  #Lucys phased data for chr 15, each row is a SNP, 1317 rows 
  markers<-read.table(paste0("AlphaPeel_files/all_chr/Deer_AlphaPeel_marker_file_chr", k, ".txt"), stringsAsFactors = FALSE)%>%
    select(V2) #reading in marker names used OR V3 = BP 
  names(phased) <- c("IID", markers$V2[-1])
  
  ##### counting number of haplotypes in defined region #####
  ###########################################################
  
  
  Hap_div_chr <- list()
  
  for (i in seq(from = 1, to = ncol(phased), by = 10)) {
    until <- ifelse((i+19) > ncol(phased),  ncol(phased), i+19) #40 SNP windows 
    if(until==i){break} #if window is only 1 snp long 
    window_temp <- phased[, c(i:(until))]
    if (ncol(window_temp)<20){break}
    window_unite <- window_temp %>% 
      unite("string_haplo", na.rm = TRUE, remove = FALSE) %>% 
      filter_all(all_vars(!grepl("9", .)))%>% #removing the non-called 9s
      group_by(string_haplo)%>% #grouping by unique haplotype
      tally() #counting number of unique haplotypes
    window_unite$ns<-NA
    for (j in 1:nrow(window_unite)){
      window_unite$ns[j]<-window_unite$n[j]*((window_unite$n[j])-1)
      
    }
    
    top<-sum(window_unite$ns)
    bottom<-sum(window_unite$n)*(sum(window_unite$n)-1)
    D<-1-(top/bottom)
    Hap_div_chr <- c(Hap_div_chr, D)
  }
  
  hap_div_df<-as.data.frame(unlist(Hap_div_chr))
  hap_div_df$window_number<-(1:nrow(hap_div_df))
  names(hap_div_df)[1]<-"haplotype_div"
  hap_div_df$CHR<-k

  
  hap_div_df$HRiD_focal=NA
  for (col in 2:nrow(hap_div_df)){

    D_minus_1=as.numeric(hap_div_df[col-1,1])
    D_plus_1=as.numeric(hap_div_df[col+1,1])
    hap_div_df$HRiD_focal[col]=(D_minus_1+D_plus_1)/(2*as.numeric(hap_div_df[col,1]))
          }
  
  
  HRiD[[k]]=hap_div_df
  
  print(paste0("Finsihed Chromosome ",k))
  
}


full_HRiD_df<-do.call(rbind.data.frame, HRiD)


mean_HRiD=as.numeric(mean(full_HRiD_df$HRiD_focal,na.rm = TRUE))
sd_HRiD <- sd(full_HRiD_df$HRiD_focal,na.rm = TRUE)


full_HRiD_df$Z_value=NA
full_HRiD_df$P_value=NA
full_HRiD_df$LogPValue=NA

for (col in 2:nrow(full_HRiD_df)){
  
  full_HRiD_df$Z_value[col] <- (full_HRiD_df$HRiD_focal[col] - mean_HRiD)/sd_HRiD
  full_HRiD_df$P_value[col] <- pnorm(-(full_HRiD_df$Z_value[col]))           # One-sided test 
  full_HRiD_df$LogPValue[col] <- -log10(full_HRiD_df$P_value[col])
    }



setwd("H:/")

#save(full_HRiD_df,HRiD, file = "PhD_4th_yr/Heredity_ROH_density_manuscript/HRiD/hap_div_drop_scan.RData")

#load("PhD_4th_yr/Heredity_ROH_density_manuscript/HRiD/hap_div_drop_scan.RData")



top=as.numeric(quantile(full_HRiD_df$HRiD_focal, c(0.998),na.rm = TRUE))
top #1.055362
z_val=(top - mean_HRiD)/sd_HRiD
p_va=pnorm(-(z_val))           # One-sided test 
log_p=-log10(p_va) #4.571651


full_HRiD_df$Order <- 1:nrow(full_HRiD_df)

full_HRiD_df$CHR<-as.numeric(full_HRiD_df$CHR)

axis.set <- full_HRiD_df %>% 
  group_by(CHR) %>% 
  summarize(center = (max(Order) + min(Order)) / 2)


ggplot(full_HRiD_df, aes(Order, LogPValue, col = as.factor(CHR%% 2), group=CHR)) +
  scale_colour_manual(values = c("gray36", "gray77")) +
  geom_line() +
  #geom_point()+
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 8, vjust = 0.5),
        axis.text.y = element_text(size = 12, vjust = 0.5))+
  scale_x_continuous(expand = c(0, 0), label = axis.set$CHR, breaks = axis.set$center)+
  labs(x="Chromosome",y="-Log(p-value)",title="Haplotype diversity drop")+
  geom_hline(yintercept = log_p, linetype="dashed", color = "red") 
  
  



chr_sub=full_HRiD_df%>%filter(CHR==15)

ggplot(chr_sub, aes(window_number, haplotype_div)) +
  geom_line()+
  theme_bw()


ggplot(full_HRiD_df, aes(window_number, haplotype_div)) +
  geom_line()+
  theme_bw()+
  facet_wrap(~CHR, scales = "free_x")#+
  geom_rect(data = subset(full_HRiD_df,CHR==15), 
            alpha = 0.6, xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf) 
