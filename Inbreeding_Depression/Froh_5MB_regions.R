library(tidyverse)


setwd("H:/")

FROH<-read.table("PhD_3rdYR/Data_files/ROH_output/ROH_search_UpdatedMb_01_2022.hom", header=T, stringsAsFactors = F)%>%
 dplyr::rename(Code=IID)%>%mutate(Start_Mb=POS1/1000000)%>%mutate(End_Mb=POS2/1000000)%>%select(Code,CHR,Start_Mb,End_Mb)

##need length of Chr from map to know when to stop