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

setwd("H:/")


for (k in 1:29){
  if (k==3)next
  if (k==5) next
  if (k==9) next
  if (k==13) next ## for failed simulations on dud machine
  if (k==19) next
  if (k==24) next
  
  hh <- data2haplohh(hap_file = paste0("PHD_2ndYR/Slim/Ash_server_output/LG15/LG15_",k,"/slim_out_sim_",k),
                   polarize_vcf = FALSE,
                   vcf_reader = "data.table", 
                   remove_multiple_markers = TRUE)


  phased_df<-as.data.frame(hh@haplo)

  Haplotype_diversity <- list()

  for (i in seq(from = 1, to = ncol(phased_df), by = 100)) {
    until <- ifelse((i+99) > ncol(phased_df),  ncol(phased_df), i+99) #10 SNP windows 
    if(until==i){break} #if window is only 1 snp long 
    window_temp <- phased_df[, c(i:(until))]
    if (ncol(window_temp)<100){break}
   window_unite <- window_temp %>% 
      unite("string_haplo", na.rm = TRUE, remove = FALSE) %>% 
      filter_all(all_vars(!grepl("9", .)))%>% #removing the non-called 9s
      group_by(string_haplo)%>% #grouping by unique haplotype
     tally() #counting number of unique haplotypes
   window_unite$ns<-NA
   for (j in 1:nrow(window_unite)){
      window_unite$ns[j]<-window_unite$n[j]*((window_unite$n[j])-1)
    
   }
  
   tops<-sum(window_unite$ns)
   bottom<-sum(window_unite$n)*(sum(window_unite$n)-1)
    D<-1-(tops/bottom)
   Haplotype_diversity <- c(Haplotype_diversity, D)
  }

  window_haplotype_div<-as.data.frame(unlist(Haplotype_diversity))
  window_haplotype_div$window_number<-(1:nrow(window_haplotype_div))
  names(window_haplotype_div)[1]<-"haplotype_div"
  window_haplotype_div$SIM<-k
  save(window_haplotype_div, Haplotype_diversity, phased_df, file = paste0("PHD_2ndYR/Slim/Haplotypes/window_haplotype_div/LG15_NOSEL_Simulated_Haplotype_div_sim_100snpwind", k, ".RData"))


}


Genomewide_hap_diversity <- NULL

for(z in 1:28){
  if (z==3) next
  if (z==5) next
  if (z==9) next
  if (z==13) next ## for failed simulations on dud machine
  if (z==19) next
  if (z==24) next
  
  load(paste0("PHD_2ndYR/Slim/Haplotypes/window_haplotype_div/LG15_NOSEL_Simulated_Haplotype_div_sim_100snpwind", z, ".RData"))
  Genomewide_hap_diversity <- rbind(Genomewide_hap_diversity, window_haplotype_div)
 
}



write.table(Genomewide_hap_diversity,
            file = "PHD_2ndYR/Slim/Haplotypes/window_haplotype_div/simulatedLG15_NOSEL_hapdiv_100snpwin.txt",
            row.names = F, quote = F, sep = "\t")



################
### now load in location of hotspots 

top_1_func<-function (x){quantile(x$pall, c(.99))}

hotspots_full<-NULL

### for local saved sims
#for (l in 1:15){

#neutral_1<-read.table(paste0("PHD_2ndYR/Slim/Output/ROH_output/neutral_m6_sim",l,"_ROH_OUT.hom.summary"),header = TRUE, stringsAsFactors=FALSE)%>%
 # select(BP,UNAFF)%>%mutate(pall=UNAFF/100)%>%setNames(.,c("POSITION","UNAFF","pall"))
#num<-as.numeric(top_1_func(neutral_1))
#neutral_1$hotspot<-NA
#neutral_1$hotspot<-ifelse(neutral_1$pall>num, "yes", "no")
#hotspots<-neutral_1%>%mutate(SNP_num = row_number())%>%
#  mutate(window_num=SNP_num/10)%>%mutate(window_number=plyr::round_any(window_num,1,ceiling))
#hotspots$SIM<-l
#hotspots_full <- rbind(hotspots_full, hotspots)

#}


## for server sims 


for (l in 1:28){
  if (l==3) next
  if (l==5) next
  if (l==9) next
  if (l==13) next ## for failed simulations on dud machine
  if (l==19) next
  if (l==24) next
  
  neutral_1<-read.table(paste0("PHD_2ndYR/Slim/Ash_server_output/LG15/LG15_",l,"/ROH_out",l,".hom.summary"),header = TRUE, stringsAsFactors=FALSE)%>%
    select(BP,UNAFF)%>%mutate(pall=UNAFF/100)%>%setNames(.,c("POSITION","UNAFF","pall"))
  
  num<-as.numeric(top_1_func(neutral_1))
  neutral_1$hotspot<-NA
  neutral_1$hotspot<-ifelse(neutral_1$pall>num, "yes", "no")
  
  hotspots<-neutral_1%>%mutate(SNP_num = row_number())%>%
    mutate(window_num=SNP_num/100)%>%mutate(window_number=plyr::round_any(window_num,1,ceiling))
  
  hotspots$SIM<-l
  hotspots_full <- rbind(hotspots_full, hotspots)
  
}

hapdiv_all_sims<-plyr::join(hotspots_full,Genomewide_hap_diversity)



######################################################################################################################
ggplot(hapdiv_all_sims, aes(hotspot, haplotype_div, fill=hotspot)) + geom_boxplot()+
  theme_classic()+
  labs(title = "Haplotype diveristy over 100 snp windows for neutral simulation of LG15")
  #geom_point(aes(x = hotspot, y = haplotype_div, fill=factor(hotspot)), position = "jitter", alpha = 0.1)
