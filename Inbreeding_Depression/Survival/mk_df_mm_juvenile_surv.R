library(tidyverse)
library(RODBC)




db<-"C:\\Users\\s1881212\\Documents\\Deer_database_2022/RedDeer2.05.accdb" #open connection
con<-odbcConnectAccess2007(db)

mum_birthyr<-sqlFetch(con, "tbllife")%>%dplyr::select(Code,BirthYear)%>%dplyr::rename(mum_birthyear=BirthYear,MumCode=Code)
mum_stat<-sqlFetch(con, "sys_HindStatusAtConception") %>% dplyr::select(-CalfBirthYear)%>%dplyr::rename(MumCode=Mum, Code=Calf)
birth_wt<-sqlFetch(con, "sys_BirthWt")%>%dplyr::select(BirthWt,Code)

life<-sqlFetch(con, "tbllife") %>% 
  dplyr::select(Code, BirthDay,BirthMonth, BirthYear,DeathDay,DeathMonth,DeathYear, Sex, MumCode, DeathType)%>%
  filter(DeathType!= c("S","A","D"))#%>%
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
  select(Code, YearLast)

census_sep$YearLast=as.numeric(census_sep$YearLast)
year_plus_census=year%>%full_join(census_sep)%>%
  mutate(yrs_since_last_Seen=(2022-YearLast))


## for all ids with missing death years fill in estimated death year when havent been seen for at least 3 yrs 
year_plus_census_filt=year_plus_census%>%filter(is.na(DeathYear))%>%
  filter(yrs_since_last_Seen>3)%>%
  mutate(DeathYear_est=YearLast+3)%>%
  select(Code,DeathYear_est)

year_plus_est_death=year%>%left_join(year_plus_census_filt)


year_plus_est_death$DeathMonth=as.numeric(year_plus_est_death$DeathMonth)
year_plus_est_death$DeathYear=as.numeric(year_plus_est_death$DeathYear)
year_plus_est_death$DeathYear_est=as.numeric(year_plus_est_death$DeathYear_est)
year_plus_est_death$BirthYear=as.numeric(year_plus_est_death$BirthYear)

year_juve=year_plus_est_death%>%
  mutate(juvenile_survival = case_when(
    ## when there is no recorded death year use the estimated death year (which is when id has not been seen for 3 years...)
    # not sure any juveniles will actually be recorded as 0 in this case as theyre killed off after 3 year anyway 
      (is.na(DeathYear) & (DeathYear_est<(BirthYear+2)))~"0",
      #if died in same year born then juvenile surv=0
      (DeathYear==BirthYear) ~ "0", 
      #if died the year after it was born juvenle surv=0
      (DeathYear==(BirthYear+1)) ~ "0", 
      #if died before May in the second year of life then did not survival to juvenile (cut off is May)
     (DeathYear==(BirthYear+2))&(DeathMonth<5) ~ "0", 
      #any ids with no recorded death month and a death year 2 yrs after birth is given NA as we dont know if i died before or after 2nd birhday (only 5 ids)
     (DeathYear==(BirthYear+2))&(is.na(DeathMonth)) ~ "NA", 
      
      ) )%>%
  # all ids that didnt get a 0 for juvenile survival gets a 1 for successfully surviving 
  mutate(juvenile_survival =replace_na(juvenile_survival,"1"))
  
  

  
  testing=year_juve%>%mutate(DeathYear_all=coalesce(DeathYear,DeathYear_est))%>%
    mutate(diff=DeathYear_all-BirthYear)%>%filter(diff>=2)

  
juvenile_surv=year_juve%>%filter(DeathType!= "S")%>% #filter out ids that were shot
  filter(DeathType!= "A")%>% ##filter out ids that were an accidental death
  filter(DeathType!= "D")%>%
  select(Code, BirthYear, Sex, MumCode, juvenile_survival)



## read in FROH values

setwd("H:/")

FROH<-read.table("PhD_4th_yr/2023_ROH_search/2021_calves_ROH_UpdatedSortedMb_032023.hom", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,CHR,KB) %>% dplyr::rename(Code=IID)%>%filter(nchar(Code)==5)

KB_perLG<-FROH%>%dplyr::group_by(Code, CHR)%>%
  dplyr::summarise(KB_chr=sum(KB))%>%
  ungroup %>% 
  complete(Code, CHR, fill = list(KB_chr = 0)) #completes for all chromosomes 


deermap <- read.csv("PhD_3rdYR/Data_files/Genome_assembly_mCerEla1.1.csv", header = T, stringsAsFactors = F)%>%
  dplyr::filter(!CHR %in% c("All","All_auto","X","unplaced"))

froh_per_chr<-plyr::join(KB_perLG,deermap)%>% mutate(chr_froh=KB_chr/length_Kb)%>%
  dplyr::select(-length_Mb,-length)%>% reshape2::dcast(Code~CHR) %>% 
  mutate(FROHsum = rowSums(.[2:34]))%>%mutate(FROH_sum_div=FROHsum/33)


colnames(froh_per_chr) <- c("Code", paste0("FROH_chr", 1:33),"FROHsum","FROH_sum_div")

FROH_full<-read.table("PhD_4th_yr/2023_ROH_search/2021_calves_ROH_UpdatedSortedMb_032023.hom.indiv", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,KB) %>% dplyr::rename(Code=IID)%>%mutate(FROH=KB/2591865)%>%filter(nchar(Code)==5)

#combine all datasets to make a dataframe 



juvenile_surv_df=juvenile_surv%>%
  left_join(mum_stat)%>%left_join(birth_wt)%>%
  left_join(mum_birthyr)%>%mutate(mum_age=BirthYear-mum_birthyear)%>%
  right_join(froh_per_chr)%>%## only joining when ids have FROH values
  mutate(mum_age_sq=mum_age^2)
##some NAs in mum stat and birth weight




write.table(juvenile_surv_df,
            file = "PhD_4th_yr/Inbreeding_depression_models/survival/juvenile_survival_df.txt",
            row.names = F, quote = F, sep = ",",na = "NA") #saving tables as txt file 


write.table(froh_per_chr,
            file = "PhD_4th_yr/Inbreeding_depression_models/Froh_per_chr.txt",
            row.names = F, quote = F, sep = "\t")



