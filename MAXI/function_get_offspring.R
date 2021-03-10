library(dplyr)
library(purrr)

### This script is a function to get all the decendents of an individual or individuals in a list ###
### Compiling all offspring of an idiv then moving to the next generation to compile all their offfspring to gte the decendents ##
## Help from Martin Stoffel ##


setwd("H:/PHD 1st YR")

ped<-read.csv("Deer data/Database/RumPedigree_birth_07_2020.csv", header = TRUE, stringsAsFactors = FALSE)%>%
  select(-SireProbability)

ped<-ped[!((ped$MumCode=="") & ped$Sire==""), ] ##getting rid of ids with no mum or dad records 
head(ped)
str(ped)

#codes <- Maxipro$Code

#get_offspring <- function(codes) {
 # maxiF <- ped %>%  
  #  filter((MumCode %in% codes) | (Sire %in% codes))
#}


get_offspring <- function(sire_name) {
  # starting conditions
  current_offspring <- ped %>% 
    filter(Sire == sire_name)
  
  count_f <- 1
  all_offspring <- list()
  all_offspring[[count_f]] <- current_offspring
  
  if (nrow(current_offspring) == 0){
    offspring_left <- FALSE
  } else {
    offspring_left <- TRUE
  }
  
  while (offspring_left) {
    # go to next generation
    count_f <- count_f + 1
    # get offspring of next generation
    current_offspring <- ped %>%  
      filter((MumCode %in% current_offspring$Code) | (Sire %in% current_offspring$Code))
    # add offspring of next generation to all offspring
    all_offspring[[count_f]] <- current_offspring
    # if the next generation does not have offspring, leave the loop
    if (nrow(current_offspring) == 0) offspring_left <- FALSE
  }
  all_offspring <- all_offspring %>% bind_rows() %>% unique() 
  all_offspring
}

# all_individuals <- unique(ped$Code)
#all_individuals <- c("MAXIS", "ROGER")
all_individuals_2<-read.table("MAXI_gencont_stuff/IDsSires_numoffspring_1965_1975.txt", header = TRUE)%>%
  select(Code) #all individuals born between 1965-1975 loaded in 
all_ids<-as.vector(unlist(all_individuals_2$Code))## to convert df to list by code



offspring <- map(all_ids, get_offspring) ## offspring is list of 55 sires which had offspring in the 10 yr period and all of their decendents 


names(offspring) <- all_ids
offspring[["MAXIS"]]  ##can see all offspring of specific sire e.g. MAXI

# count number of offspring


##################################################################################
########~~~ now need to get number of offspring per year for each list element ######
##################################################################################



full_off_peryr <- list()

for(i in 1:55){
  full_off_peryr[[i]]<-offspring[[i]] %>% dplyr::select(Code, BirthYear)%>%
    dplyr::count(BirthYear)
  
} ## 

names(full_off_peryr) <- all_ids # replaces list numbers with names of sires 

full_off_peryr[["MAXIS"]] 



###############################################################################
########### now work out contribution per year compared to full cohort ########
###############################################################################

library(RODBC)
db<-"C:\\Users\\s1881212\\Documents\\Deer database 07_2020\\RedDeer1.93.accdb"
con<-odbcConnectAccess2007(db)
full_ped<-sqlFetch(con, "tbllife") %>% select(Code, BirthYear)
odbcClose(con)

yr_fullped<-full_ped%>%
  dplyr::count(BirthYear) 
names(yr_fullped)[2]<-c("full_cohort")


contribution_per_yr<-list()
library(plyr)
for (i in 1:55){
  contribution_per_yr[[i]]<-join(full_off_peryr[[i]], yr_fullped)
}    


join_func<-function(x){
  join(x, yr_fullped)
}

full_ped_joined<-map(full_off_peryr, join_func) ### this actually worked, more efficient that above 

###########################################################################
#### now to work out contribution #####################################

contribution<- function(x){
  x$cont<-NA
for (i in 1:nrow(x)){
  x$cont[i]<-(x$n[i]/x$full_cohort[i]) ## this is null
  
}
}
  
contribution<-function(sire){
  cont<-(sire)$n/(sire)$full_cohort 
  
  ## this gives correct numbers but not in a list
  
}

map(full_ped_joined, contribution)




fill<-n/full_cohort
was<-Map(cbind, full_ped_joined, contriibution=fill)

contr_func<-function(x){
  x$conttrbution<-x$n/x$full_cohort 
  
  
  ## this gives correct numbers but not in a list
  
}

contribution<- function(x){
    for (i in 1:nrow(x)){
    x$cont[i]<-(x$n[i]/x$full_cohort[i]) ## this is null
    
  }
}

who<-map(full_ped_joined, contribution)


maybe<-map(full_ped_joined, contr_func)


#### need to fix this ::
#contribution_per_yr[[i]]$contribution<-NA
#for (j in 1:nrow(contribution_per_yr[[i]]))
 # full_off_peryr[[i]]$contribution<-full_off_peryr[[i]]$n[j]/full_off_peryr[[i]]$full_cohort[j]







###### unlisting #####




un<-bind_rows(full_ped_joined, .id = "Name") ## unlists and bind rows together using Name as element names



un$contribution<-NA
for (i in 1:nrow(un)){
  
  un$contribution[i]<-un$n[i]/un$full_cohort[i]
  
}


library(ggplot2)
ggplot(un, aes(x = BirthYear, y = contribution, colour = "Name")) + geom_line() +
  labs(x = "Cohort year", y = "Genealogical contribution", title = "Proportion of cohorts genetically decendent from MAXI (1974 - 2019)") +
  theme_bw()


ggplot(un, aes(x = BirthYear, y = contribution)) + 
  geom_line(aes(linetype = Name))+
  theme(legend.position = "none")


ggplot(un, aes(x=BirthYear, y=contribution, group = Name, colour = Name)) +
  geom_line() +
  theme(legend.position = "none")

  geom_point( size=4, shape=21, fill="white")