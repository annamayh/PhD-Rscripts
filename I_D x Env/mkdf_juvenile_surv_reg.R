library(RODBC)
library(tidyverse)



db<-"C:\\Users\\s1881212\\Documents\\Deer_database_2022/RedDeer2.05.1.accdb" #open connection
con<-odbcConnectAccess2007(db)

mum_stat<-sqlFetch(con, "sys_HindStatusAtConception") %>%dplyr::select(-CalfBirthYear)%>%dplyr::rename(MumCode=Mum, Code=Calf)
#ped<-sqlFetch(con, "sys_Pedigree")%>%dplyr::select(Code, BirthYear, MumCode)
census<-sqlFetch(con, "tblCensus")%>%dplyr::select(Date,Code,Northing,Easting)

mum_birthyr<-sqlFetch(con, "tbllife")%>%dplyr::select(Code,BirthYear)%>%dplyr::rename(mum_birthyear=BirthYear,MumCode=Code)
birth_wt<-sqlFetch(con, "sys_BirthWt")%>%dplyr::select(BirthWt,Code)

life<-sqlFetch(con, "tbllife") %>% 
  dplyr::select(Code, BirthDay,BirthMonth, BirthYear,DeathDay,DeathMonth,DeathYear, Sex, MumCode, DeathType)

cen=sqlFetch(con, "sys_LastCensusSighting")
odbcClose(con)
#####################
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
life$BirthYear=as.factor(life$BirthYear)


mum_loc_calf_yr=ave_N_E_calf_yr%>%
  rename(MumCode=Code, BirthYear=pseudo_calf_yr)%>% #so we get the average location of the mum during the birth year of the calf
  inner_join(life)%>%
  select(Code, MumCode, BirthYear, N_calf_yr, E_calf_yr)


#################################################################################################
#### survival data

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


#year_plus_est_death$DeathMonth=as.numeric((year_plus_est_death$DeathMonth))
year_plus_est_death$DeathYear=as.numeric(as.character(year_plus_est_death$DeathYear))
year_plus_est_death$DeathYear_est=as.numeric(as.character(year_plus_est_death$DeathYear_est))
year_plus_est_death$BirthYear=as.numeric(as.character(year_plus_est_death$BirthYear))
year_plus_est_death$YearLast=as.numeric(as.character(year_plus_est_death$YearLast))

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


##### ##################
## work out DOB cont ###
########################

# first birth date is 30th/4
# Last is 10th Oct
day_seq=read.csv("PhD_4th_yr/Inbreeding_depression_models/date_seq.csv", header = TRUE)
head(day_seq)
DOBs=year_juve%>%select(Code, BirthDay, BirthMonth)%>%
  left_join(day_seq)%>% ## sequences of day starting at 30th April 
  select(Code, Day_seq)



table(year_juve$DeathType)


juvenile_surv=year_juve%>%filter(DeathType!= "S"| is.na(DeathType))%>% #filter out ids that were shot
  filter(DeathType!= "A" | is.na(DeathType)) %>% #%>% ##filter out ids that were an accidental death
  filter(DeathType!= "D" | is.na(DeathType))%>% ## keeping NA death type
  dplyr::select(Code, BirthYear, Sex, MumCode, juvenile_survival)

table(juvenile_surv$juvenile_survival)

# want to keep ids that were shot but also made it to adulthood 
juvenile_surv_shot=year_juve%>%filter(DeathType== "S"&juvenile_survival=="1")%>%
  select(Code, BirthYear, Sex, MumCode, juvenile_survival)

table(juvenile_surv_shot$juvenile_survival)

## adding ids that survived to age 2 then were shot to the df
juvenile_surv=juvenile_surv%>%rbind(juvenile_surv_shot)
table(juvenile_surv$juvenile_survival)






## read in FROH values ##

setwd("H:/")

FROH_full<-read.table("PhD_4th_yr/2023_ROH_search/2021_sleuthed_052023.hom.indiv", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,KB) %>% dplyr::rename(Code=IID)%>%mutate(FROH=KB/2591865)%>%
  filter(nchar(Code)==5)%>% #removing IDs with non-sensical ID codes
  select(-KB)

FROH_mum<-read.table("PhD_4th_yr/2023_ROH_search/2021_sleuthed_052023.hom.indiv", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,KB) %>% dplyr::rename(MumCode=IID)%>%mutate(MumFROH=KB/2591865)%>%
  filter(nchar(MumCode)==5)%>% #removing IDs with non-sensical ID codes
  select(-KB)

mum_loc_calf_yr$BirthYear=as.numeric(as.character(mum_loc_calf_yr$BirthYear))


surv_loc_df=juvenile_surv%>%
    left_join(mum_stat)%>%left_join(birth_wt)%>%
    left_join(mum_birthyr)%>%
    mutate(mum_age=BirthYear-mum_birthyear)%>%
    mutate(mum_age_sq=mum_age^2)%>%
    left_join(mum_loc_calf_yr)%>%
    left_join(FROH_full)%>%
    left_join(FROH_mum)%>%
    rename(N=N_calf_yr, E=E_calf_yr)%>%
    left_join(DOBs)%>%
      select(-mum_birthyear)#dont need mum birth yr in df anymore
  



LocToReg6 <- function(E, N) {
  ifelse(N < 8019, "SG", #south glen
         ifelse(E < 1361, "LA", # laundry greens
                ifelse(E < 1366,
                       ifelse(N > 8033, "NG", "MG"), # north glen, mid  glen 
                       ifelse(E < 1373 , "IM", "SI")))) }


surv_loc_df$Reg <- with(surv_loc_df, LocToReg6(E, N))


write.table(surv_loc_df,
            file = "PhD_4th_yr/Spatial_var_inbreeding/survival_loc.txt",
            row.names = F, quote = F, sep = ",",na = "NA") #saving tables as txt file 



  