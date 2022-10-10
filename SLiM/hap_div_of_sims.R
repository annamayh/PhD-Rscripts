#library(tidyverse)
library(tidyverse)
library(janitor)

setwd("H:/")

#doesnt work for sim num 22

all_sims_hap_div=list()



sim_num=1

for (sim_num in 1:25){
  
      # if (sim_num==7)next #'for failed sims on dud machine'
      # if (sim_num==15)next
      # if (sim_num==22)next
      # if (sim_num==24)next
      # 
        
    
    sim=read.table(paste0("PhD_4th_yr/Heredity_ROH_density_manuscript/hap_div_sims/Rum_neutral_output/SLiM_output_sim_",sim_num))%>%
      select(-V1,-V3,-V4, -V5, -V6, -V7, -V8, -V9)#%>%#removing other random columns you get from vcf files
    
    rand=sample(1:nrow(sim), 1634)## 1634 = average snp density for searched ROH in real dataset
    
    sims_sub=sim[(rand),] #subsetting the number of snps to give the average num of snps in real data
    
    sims_snps=sims_sub%>%#mutate(V2 = ifelse(duplicated(V2), paste0(V2,"_2"),V2))%>% ## making snp names unique for duplicated names
      t()%>%as_tibble()  ##trnaspose so SNPs are columns and IDs are rows
    
    haps=sims_snps%>%row_to_names(row_number = 1)
    phased_haps=haps%>%separate_rows(1:nrow(haps), sep="|") #seperating into phased haps
    phased_sims<-phased_haps[!(phased_haps[,1]=="" | phased_haps[,1]=="|"),] #removing the rows with the spaces cos didnt work properly 
    phased_lab=phased_sims
    
    ## hap div function taken from own github and modified slightly
    Haplotype_diversity <- list()
    
    for (i in seq(from = 1, to = ncol(phased_lab), by = 10)) {
          until <- ifelse((i+9) > ncol(phased_lab),  ncol(phased_lab), i+9) #10 SNP windows 
          if(until==i){break} #if window is only 1 snp long 
          window_temp <- phased_lab[, c(i:(until))]
          if (ncol(window_temp)<10){break}
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
      all_snps=colnames(phased_lab[i:until])## all snps in the window
      focal=matrix(c(D, all_snps), nrow=1)#combining into matrix to rbind into columns
      Haplotype_diversity <-(rbind(Haplotype_diversity,focal))
    }
    
    #######
    ## read in ROH to fund hotspots 
    ROH=read.table(paste0("PhD_4th_yr/Heredity_ROH_density_manuscript/hap_div_sims/Rum_neutral_output/ROH_scan_sim",sim_num,".hom.summary"), header=T, stringsAsFactors = F)
    hot_threshold=as.numeric(quantile(ROH$UNAFF, .99))
    ##set uo if else for yes no over threshold
    hotspot_snp_bp=ROH%>%filter(UNAFF>=hot_threshold)%>%pull(BP)%>%as.character()##filtering for only snps at hotspots
    ## need to search each row for the wheher it includes a hotspot snp or not
    hotspots_to_search=paste(hotspot_snp_bp, collapse = "|")
    
    view=as.data.frame(hotspot_snp_bp)
    
    Hap_div_df=as.data.frame(Haplotype_diversity) #duplicating to use in loop later
    Hap_div_df$hotspot=NA
    
    
    for (j in 1:nrow(Hap_div_df)){
    
        row=Hap_div_df[j,] ## pull out one rofrom haplotype div table
        match=grepl(hotspots_to_search, row)#searching whether any hotspot snps match the first row of the table
        Hap_div_df$hotspot[[j]]=ifelse(length(match[match==TRUE]), "yes","no") # adding in column whih states whether window has a hotspot snp in it 
    }
    
    
    ##making it pretty df
    pretty_hap_div=Hap_div_df%>%select(V1, hotspot)
    pretty_hap_div$window_number<-(1:nrow(pretty_hap_div))
    names(pretty_hap_div)[1]<-"D"
    pretty_hap_div$D=as.numeric(pretty_hap_div$D)
    
    all_sims_hap_div[[sim_num]]=pretty_hap_div

}



############################################################################################################################
########## plot ############################################################################################################
########################################################################################


all_sims_hap_div_2=all_sims_hap_div[-c(7,15,22)]

all_sims_unlisted=do.call(rbind.data.frame, all_sims_hap_div_2)


ggplot(all_sims_unlisted, aes(hotspot, D, fill=hotspot)) + geom_boxplot()+
  theme_classic()+
  scale_x_discrete(labels = c("Remainder of genome","ROH Hotspots" ))+
 # theme(legend.position = "none", 
  #      axis.title.x=element_blank())+
  ylab("Window haplotype diversity")+
  scale_fill_brewer(palette="Dark2")+
  geom_point(aes(x = hotspot, y = D, fill=factor(hotspot)), position = "jitter", alpha = 0.05)


library(wesanderson)


ggplot(all_sims_unlisted,aes(x=hotspot, y=D, fill=hotspot))+
  geom_violin(position="dodge", alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values = wes_palette("Darjeeling1", n=2))+
  stat_summary(fun=mean, geom="point",size=2)


