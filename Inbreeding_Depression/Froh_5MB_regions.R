library(tidyverse)
library(tibble)


setwd("H:/")

FROH<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom", header=T, stringsAsFactors = F)%>%
 dplyr::rename(Code=IID)%>%mutate(Start_Mb=POS1/1000000)%>%
  mutate(End_Mb=POS2/1000000)%>%select(Code,CHR,Start_Mb,End_Mb)%>%
  mutate(length_Mb=End_Mb-Start_Mb)

names=as.data.frame(unique(FROH$Code))%>%rename("Code"="unique(FROH$Code)")#for joining later


all_chrs_froh_split=list()


for (chr in 1:33){

    focal_chr=FROH%>%filter(CHR==chr)
    froh_split_full=list()
    kb_in_roh=list()
    focal_chr_end=max(focal_chr$End_Mb)
    focal_chr_last_window=plyr::round_any(focal_chr_end, 10, f=floor)
    
    for (row in 1:nrow(focal_chr)){
    
      kb_in_roh=list()
      
      for (i in seq(from = 1, to = 140, by = 10)){  #longest chr is 140Mb 
      window_start_diff=focal_chr$Start_Mb[[row]]-i
      ROH_length=focal_chr$length_Mb[[row]]
      window_end_diff=focal_chr$End_Mb[[row]]-i
      
      # if start of window is before when roh starts, skip to next window 
      if(window_start_diff>10){next
      #if window encompasses the end of the chr
      }else if(i==focal_chr_last_window){
        kb_in_roh[[i]]=window_end_diff/(focal_chr_end-i)
      
      #if ROH starts within window and ROH only spans one window 
      }else if(window_start_diff<10&window_start_diff>0&window_end_diff<10&ROH_length<10){
        kb_in_roh[[i]]=ROH_length/10
      
      # if ROH starts within window but spans multiple windows 
      } else if(window_start_diff<10&window_start_diff>0&window_end_diff>10){
        kb_in_roh[[i]]=(10-window_start_diff)/10
      
        #If ROH started in window before but ends in the current window 
      }else if(window_start_diff<0&window_start_diff<0&window_end_diff<10&window_end_diff>0) {
        second_window=i-focal_chr$End_Mb[row]
        kb_in_roh[[i]]=(-second_window/10)  
      
      #if ROH starts in window and ends in next/other window
      }else if(window_start_diff<0&window_end_diff>10&ROH_length>10) {
        kb_in_roh[[i]]=1} else{
          break}
        }
    
      names(kb_in_roh) <- seq_along(kb_in_roh)#converts window starts to name of list 
      froh_split=as.data.frame(kb_in_roh%>%
                                 compact()%>% #removes unused list objects 
                                 unlist())%>% # unlisting object
                                rownames_to_column() # makes list names a column (where window starts) 
      
      ID=focal_chr$Code[[row]] ##getiing ID of deer
      froh_split$Code=ID #Making ID a column
    
      colnames(froh_split) <- c("Window_start","Kb_in_roh","Code")
      
      froh_split_full[[row]]=froh_split
      
      }
    
    
    focal_chr_froh_split=do.call(rbind.data.frame, froh_split_full)%>% 
          rbind(1)%>% #this is a stupid way of getting around the fact that some chr dont have any ROH in the first window. 
          rbind(11)%>% #chr 20 does have any ROH for first 2 columns
          group_by(Code,Window_start)%>% ## for times when are 2 ROH in same window w/ space inbtw
          summarise(Kb_in_roh_all = sum(Kb_in_roh), .groups = "keep")%>%
          ungroup()%>%
          complete(Window_start, Code)%>% #completing for windows with no ROH
          cbind(CHR=chr)%>%# add chr name for ease
          mutate(Window_start=as.numeric(Window_start))%>% #making window number numeric so in right order
          arrange(Code,Window_start)%>%
          subset(Code!="1" & Code!="11")
   
    
    focal_chr_froh_split$label=paste("Froh_c",focal_chr_froh_split$CHR,"w",focal_chr_froh_split$Window_start, sep = "_")
    
    
    #Getting into df with code as column and froh split as diff columns 
    M=focal_chr_froh_split%>%select(-CHR, -Window_start)%>%
      pivot_wider(names_from = label, values_from = Kb_in_roh_all)%>%
      dplyr::full_join(names) %>%##joining all names together fo ids with no ROH on chr
      plyr::arrange(Code)
      
    M[is.na(M)] <- 0
    
    all_chrs_froh_split[[chr]]=M
    
    print(paste0("Finished chr ", chr))
      
}    


#sub=all_chrs_froh_split[[3]]
#sub2=all_chrs_froh_split[[4]]


all_chr_split_df=do.call(cbind.data.frame, all_chrs_froh_split)%>% #problem is here!! with the cbinding ... names are out of order or smth
  subset(., select = which(!duplicated(names(.))))%>%  ###removing multiple cols of Code 
  mutate(FROH_split_sum = rowSums(.[2:263]))%>%
  mutate(FROH_split_sum_div = FROH_split_sum/33)
  



save(all_chr_split_df, file="PhD_4th_yr/Inbreeding_depression_models/Froh_split_regions/all_chr_split10Mb_df.RData")




