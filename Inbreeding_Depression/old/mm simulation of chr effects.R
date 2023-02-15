library(MCMCglmm)
library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)


M<-matrix(runif(1000*33), 1000, 33)
colnames(M)<-paste0("Chr", 1:33)

beta_c<--0.13 #average variance (effect) of all chromosomes
sigma2_c<-0.02 #average variance of each individual chromosome 
sigma2_e<-1.5 #average residual var

ef<-rnorm(33, beta_c, sqrt(sigma2_c))#generating distributions based on average effects specified

view_ef<-as.data.frame(ef)

tot_ef<-matrix(0, 1000,1)

for(i in 1:33){
  tot_ef<-tot_ef+ef[i]*M[,i]
} #total effect 

#tot_ef<-M%*%ef## doing the same thing as above 

y_var<-rnorm(1000, tot_ef, sqrt(sigma2_e))# creating response variable 
y_mean<-rnorm(1000, 6.5, sqrt(1.6))
y<-y_mean + y_var ## 6.63 decided based on mean birth weight of 6.5 + estimated residuual variance estimated

dat<-data.frame(M, y=y, sumM=rowSums(M)) ## df with M = effect of each chr, y= response variable and sum M = sum of each chr effect

mean(dat$y)

m1<-MCMCglmm(y~sumM, random=~idv(Chr1+Chr2+Chr3+Chr4+Chr5+Chr6+Chr7+Chr8+
                                   Chr9+Chr10+Chr11+Chr12+Chr13+Chr14+Chr15+Chr16+
                                   Chr17+Chr18+Chr19+Chr20+Chr21+Chr22+Chr23+Chr24+
                                   Chr25+Chr26+Chr27+Chr28+Chr29+Chr30+Chr31+Chr32+
                                   Chr33), data=dat,pr=TRUE, nitt = 50000)

## can see from summary that the estimates are pretty close to what we had simulated
# for average effect etc 


plot(m1)

summary(m1)

names <- apply(m1$Sol,2,mean) %>% names
sols<-apply(m1$Sol,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(m1$Sol,2,quantile,probs = c(0.95)) #gets upper confidence interval for all solutions 
CI_lower<-apply(m1$Sol,2,quantile,probs = c(0.05)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols<-Random_table%>%filter(model_variable %like% "Chr")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH

FROH_sols$CHR<-as.factor(FROH_sols$CHR)

ggplot(data=FROH_sols, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange(colour="grey50") + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Linkage group", y="solution + CI", title = "Effect size of Linkage group FROH on Birthweight (kg)") +
  theme_bw()+  # use a white background                 
  theme(legend.position = "none")



birth_wt_stats<-birth_wt%>%na.omit()

mean(birth_wt_stats$BirthWt)
var(birth_wt_stats$BirthWt)
