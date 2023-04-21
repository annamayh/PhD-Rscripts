library(tidyverse)
library(tibble)





setwd("H:/")

FROH<-read.table("PhD_4th_yr/2023_ROH_search/2021_calves_ROH_UpdatedSortedMb_032023.hom", header=T, stringsAsFactors = F)%>%
  dplyr::rename(Code=IID)%>%
  filter(nchar(Code)==5)%>%
  mutate(Start_Mb=POS1/1000000)%>%
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


sub=all_chrs_froh_split[[3]]
#sub2=all_chrs_froh_split[[4]]


all_chr_split_df=do.call(cbind.data.frame, all_chrs_froh_split)%>% #problem is here!! with the cbinding ... names are out of order or smth
  subset(., select = which(!duplicated(names(.))))%>%  ###removing multiple cols of Code 
  mutate(FROH_split_sum = rowSums(.[2:263]))%>%
  mutate(FROH_split_sum_div = FROH_split_sum/33)



#############################################################################################
###### now add survival data ##################################################
##################################################################################################
library(RODBC)



db<-"C:\\Users\\s1881212\\Documents\\Deer_database_2022/RedDeer2.05.accdb" #open connection
con<-odbcConnectAccess2007(db)

mum_birthyr<-sqlFetch(con, "tbllife")%>%dplyr::select(Code,BirthYear)%>%dplyr::rename(mum_birthyear=BirthYear,MumCode=Code)
mum_stat<-sqlFetch(con, "sys_HindStatusAtConception") 
mum_stat=mum_stat%>%dplyr::select(-CalfBirthYear)%>%dplyr::rename(MumCode=Mum, Code=Calf)
birth_wt<-sqlFetch(con, "sys_BirthWt")%>%dplyr::select(BirthWt,Code)

life<-sqlFetch(con, "tbllife") %>% 
  dplyr::select(Code, BirthDay,BirthMonth, BirthYear,DeathDay,DeathMonth,DeathYear, Sex, MumCode, DeathType)
# drop_na(BirthYear,DeathYear,DeathMonth,BirthMonth)%>%
# join(mum_birthyr)%>%mutate(mum_age=BirthYear-mum_birthyear)%>%
# join(mum_stat)%>%join(birth_wt)

ped<-sqlFetch(con, "sys_Pedigree")%>%dplyr::select(Code, MumCode,Sire)

cen=sqlFetch(con, "sys_LastCensusSighting")

odbcClose(con)

##################################
#### calculating juvenile survival###
######################################


##based on last census sighting 

year=life

census_sep=cen%>%separate(LastSeen,into=c("YearLast","MonthLast","DayLast"),sep="-")%>%
  dplyr::select(Code, YearLast)

census_sep$YearLast=as.numeric(census_sep$YearLast)
year_plus_census=year%>%right_join(census_sep)%>%
  mutate(yrs_since_last_Seen=(2023-YearLast))

year_plus_census_check=year_plus_census%>%filter(is.na(DeathYear))%>%filter(!is.na(BirthYear))%>%
  filter(BirthYear==YearLast)

## for all ids with missing death years fill in estimated death year when havent been seen for at least 3 yrs 
year_plus_census_filt=year_plus_census%>%filter(is.na(DeathYear))%>%
  mutate(DeathYear_est=YearLast+3)%>%
  filter(yrs_since_last_Seen>1)%>% #has to be more than 3 yrs scince its been seen to have an estimated death
  filter(DeathYear_est<2023)%>% 
  dplyr::select(Code,DeathYear_est,YearLast)

year_plus_est_death=year%>%left_join(year_plus_census_filt) ##adding in an estimated death yr and year last seen for ids w/ no death year (not including ids born in the last 3 yrs )


year_plus_est_death$DeathMonth=as.numeric(year_plus_est_death$DeathMonth)
year_plus_est_death$DeathYear=as.numeric(year_plus_est_death$DeathYear)
year_plus_est_death$DeathYear_est=as.numeric(year_plus_est_death$DeathYear_est)
year_plus_est_death$BirthYear=as.numeric(year_plus_est_death$BirthYear)
year_plus_est_death$YearLast=as.numeric(year_plus_est_death$YearLast)

# year_plus_est_death_check=year_plus_est_death%>%filter(is.na(DeathYear))%>%filter(!is.na(BirthYear))%>%
#   filter(BirthYear==YearLast)

year_juve=year_plus_est_death%>%
  mutate(juvenile_survival = case_when(
    ((is.na(DeathYear)) & (BirthYear==YearLast))~"0",
    
    ## when there is no recorded death year use the estimated death year (which is when id has not been seen for 3 years...)
    # not sure any juveniles will actually be recorded as 0 in this case as theyre killed off after 3 year anyway 
    (is.na(DeathYear) & (DeathYear_est<(BirthYear+2)))~"0",
    #if died in same year born then juvenile surv=0
    (DeathYear==BirthYear) ~ "0", 
    #if died the year after it was born juvenle surv=0
    (DeathYear==(BirthYear+1)) ~ "0", 
    #if died before May in the second year of life then did not survival to juvenile (cut off is May)
    (DeathYear==(BirthYear+2))&(DeathMonth<5) ~ "0", 
    #any ids with no recorded death month and a death year 2 yrs after birth is given NA as we dont know if id died before or after 2nd birhday (only 5 ids)
    (DeathYear==(BirthYear+2))&(is.na(DeathMonth)) ~ "NA", 
    
  ) )%>%
  # all ids that didnt get a 0 for juvenile survival gets a 1 for successfully surviving 
  mutate(juvenile_survival =replace_na(juvenile_survival,"1"))



# 
# testing=year_juve%>%mutate(DeathYear_all=coalesce(DeathYear,DeathYear_est))%>%
#   mutate(diff=DeathYear_all-BirthYear)%>%filter(diff>=2)

table(year_juve$DeathType)


juvenile_surv=year_juve%>%filter(DeathType!= "S"| is.na(DeathType))%>% #filter out ids that were shot
  filter(DeathType!= "A" | is.na(DeathType)) %>% #%>% ##filter out ids that were an accidental death
  filter(DeathType!= "D" | is.na(DeathType))%>% ## keeping NA death type
  dplyr::select(Code, BirthYear, Sex, MumCode, juvenile_survival)

table(juvenile_surv$juvenile_survival)

# want to keep ids that were shot but also made it to adulthood 
juvenile_surv_shot=year_juve%>%filter(DeathType== "S"&juvenile_survival=="1")%>%
  select(Code, BirthYear, Sex, MumCode, juvenile_survival)


juvenile_surv=juvenile_surv%>%rbind(juvenile_surv_shot)
table(juvenile_surv$juvenile_survival)





juv_surv_10MB=juvenile_surv%>%inner_join(all_chr_split_df)%>%
  left_join(mum_stat)%>%#left_join(birth_wt)%>%
  left_join(mum_birthyr)%>%mutate(mum_age=BirthYear-mum_birthyear)%>%
  mutate(mum_age_sq=mum_age^2)##2871 ids in df




write.table(juv_surv_10MB,
            file = "PhD_4th_yr/Inbreeding_depression_models/survival/10MBsplit_juvenile_survival_df.txt",
            row.names = F, quote = F, sep = ",",na = "NA") #saving tables as txt file 


