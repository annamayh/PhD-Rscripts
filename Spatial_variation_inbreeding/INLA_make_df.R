library(RODBC)
library(tidyverse)



db<-"C:\\Users\\s1881212\\Documents\\Deer_database_2022/RedDeer2.05.accdb" #open connection
con<-odbcConnectAccess2007(db)
ped<-sqlFetch(con, "sys_Pedigree")%>%dplyr::select(Code, BirthYear, MumCode)

life<-sqlFetch(con, "tbllife") %>% 
  dplyr::select(Code,Sex)

census<-sqlFetch(con, "tblCensus")%>%dplyr::select(Date,Code,Northing,Easting)

odbcClose(con)

head(census)



census_yr=census%>%separate(Date,into=c("Year","month","Day"),sep="-")


census_yr$month=as.numeric(census_yr$month)
census_yr$half_of_year=NA
census_yr$half_of_year=ifelse(census_yr$month<5, "first", "last")


ave_N_E_yr_split=census_yr%>%
  group_by(Code,Year,half_of_year)%>% #grouping by ID, year of census and whether its the first or last half of the yr
  summarise(N_mean=mean(Northing), E_mean=mean(Easting))%>% #mean per id for first and last half of yr
  ungroup()

ave_N_E_yr_split$Year=as.numeric(ave_N_E_yr_split$Year)

ave_N_E_calf_yr=ave_N_E_yr_split%>% #creating a new variable 'pseudo year' which is the year from May-April of the calfs first year of life
  mutate(pseudo_calf_yr=case_when(half_of_year=="first" ~ Year-1, # when its the first half of the yr is taken as the year previous as this is where the calf will be 
                             half_of_year=="last" ~ Year))%>% # 
  select(-half_of_year,-Year)%>%
  group_by(Code, pseudo_calf_yr)%>% #grouping by ID and 'calf year'
  summarise(N_calf_yr=mean(N_mean, na.rm=TRUE), E_calf_yr=mean(E_mean, na.rm = TRUE)) %>%#taking mean location of calf year
  ungroup()%>%
  na.omit()




#now to join location info to calf and add in FROH 

ave_N_E_calf_yr$pseudo_calf_yr=as.factor(ave_N_E_calf_yr$pseudo_calf_yr)
ped$BirthYear=as.factor(ped$BirthYear)


mum_loc_calf_yr=ave_N_E_calf_yr%>%
  rename(MumCode=Code, BirthYear=pseudo_calf_yr)%>% #so we get the average location of the mum during the birth year of the calf
  inner_join(ped)


## read in FROH values ##

setwd("H:/")

FROH_full<-read.table("PhD_4th_yr/2023_ROH_search/2021_calves_ROH_UpdatedSortedMb_032023.hom.indiv", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,KB) %>% dplyr::rename(Code=IID)%>%mutate(FROH=KB/2591865)%>%
  filter(nchar(Code)==5)%>% #removing IDs with non-sensical ID codes
  select(-KB)

FROH_loc_full_df=mum_loc_calf_yr%>%right_join(FROH_full)%>%
  rename(N=N_calf_yr, E=E_calf_yr)%>%
  na.omit() %>%#removng all ids where mum / mum loctation no known ~250 ids from dataset
  mutate(year_cont=BirthYear-min(BirthYear)) 



write.table(FROH_loc_full_df,
            file = "PhD_4th_yr/Spatial_var_inbreeding/INLA/Location_plus_FROH.txt",
            row.names = F, quote = F, sep = ",",na = "NA") #saving tables as txt file 




