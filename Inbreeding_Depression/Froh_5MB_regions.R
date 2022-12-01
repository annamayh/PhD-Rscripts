library(tidyverse)
library(tibble)


setwd("H:/")

FROH<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom", header=T, stringsAsFactors = F)%>%
 dplyr::rename(Code=IID)%>%mutate(Start_Mb=POS1/1000000)%>%mutate(End_Mb=POS2/1000000)%>%select(Code,CHR,Start_Mb,End_Mb)%>%
  mutate(length_Mb=End_Mb-Start_Mb)

chr=1


#all_chrs_froh_split=list()


for (chr in 1:33){

    
    
    focal_chr=FROH%>%filter(CHR==chr)
    froh_split_full=list()
    
    
    for (row in 1:nrow(focal_chr)){
    
      kb_in_roh=list()
      
      
      for (i in seq(from = 0, to = 140, by = 10)){  #longest chr is 140Mb 
      
      window_start_diff=focal_chr$Start_Mb[[row]]-i
      ROH_length=focal_chr$length_Mb[[row]]
      window_end_diff=focal_chr$End_Mb[[row]]-i
      
      # if start of window is before when roh starts, skip to next window 
      if(window_start_diff>10){next
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
          group_by(Code,Window_start)%>% ## for times when are 2 ROH in same window w/ space inbtw
          summarise(Kb_in_roh_all = sum(Kb_in_roh), .groups = "keep")%>%
          ungroup()%>%
          complete(Window_start, Code)%>% #completeing for windows with no ROH
          cbind(CHR=chr)%>%# add chr name for ease
          arrange(Code)
   
    focal_chr_froh_split[is.na(focal_chr_froh_split)] <- 0
          
    #all_chrs_froh_split[[chr]]=focal_chr_froh_split
    
    print(paste0("Finished chr ", chr))
      
}    

