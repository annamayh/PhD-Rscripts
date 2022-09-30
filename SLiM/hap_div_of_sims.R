#library(tidyverse)
library(tidyverse)
library(janitor)

setwd("H:/")



sim=read.table("PhD_4th_yr/Heredity_ROH_density_manuscript/sims_less_SNPs/Rum_neutral_SLiM_output_sim_1")%>%
  select(-V1,-V3,-V4, -V5, -V6, -V7, -V8, -V9)#%>%#removing other random columns you get from vcf files


sims_snps=sim%>%mutate(V2 = ifelse(duplicated(V2), paste0(V2,"_2"),V2))%>% ## making snp names unique for duplicated names
  t()%>%as_tibble()  ##trnaspose so SNPs are columns and IDs are rows

haps=sims_snps%>%row_to_names(row_number = 1)
phased_haps=haps%>%separate_rows(1:nrow(haps), sep="|")

phased_sims<-phased_haps[!(phased_haps$`18`=="" | phased_haps$`18`=="|"),]

phased_lab=phased_sims
## hap div function taken from own github and modified slightly
Haplotype_diversity <- list()

for (i in seq(from = 1, to = ncol(phased_lab), by = 100)) {
  until <- ifelse((i+99) > ncol(phased_lab),  ncol(phased_lab), i+99) #10 SNP windows 
  if(until==i){break} #if window is only 1 snp long 
  window_temp <- phased_lab[, c(i:(until))]
  if (ncol(window_temp)<100){break}
  window_unite <- window_temp %>% 
    unite("string_haplo", na.rm = TRUE, remove = FALSE) %>% 
    group_by(string_haplo)%>% #grouping by unique haplotype
    tally() #counting number of unique haplotypes
  window_unite$ns<-NA
  for (j in 1:nrow(window_unite)){
    window_unite$ns[j]<-window_unite$n[j]*((window_unite$n[j])-1)
    
  }
  
  tops<-sum(window_unite$ns)
  bottom<-sum(window_unite$n)*(sum(window_unite$n)-1)
  D<-1-(tops/bottom)
  start_snp=colnames(phased_lab[i])## adding in start snp tag
  end_snp=colnames(phased_lab[until]) #end snp tag
  focal=matrix(c(D,start_snp,end_snp ), nrow=1)#combining into matrix to rbind into columns
  Haplotype_diversity <- as.data.frame(rbind(Haplotype_diversity,focal))
}

Haplotype_diversity$window_number<-(1:nrow(Haplotype_diversity)) #maybe dont need the window number but we'll see

names(Haplotype_diversity)=c("hap_div","start_snp","end_snp","window_num")#sorting out column names 






ROH=read.table("PhD_4th_yr/Heredity_ROH_density_manuscript/sims_less_SNPs/ROH_scan_sim_1_full_no_mins.hom.summary", header=T, stringsAsFactors = F)

hot_threshold=as.numeric(quantile(ROH$UNAFF, .99))
##set uo if else for yes no over threshold

hotspot_snp_bp=ROH%>%filter(UNAFF>=hot_threshold)%>%select(BP) ##filtering for only snps at hotspots




## now need to 
