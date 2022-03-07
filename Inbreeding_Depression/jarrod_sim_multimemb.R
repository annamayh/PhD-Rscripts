M<-matrix(runif(1000*33), 1000, 33)
colnames(M)<-paste0("C", 1:33)
  
beta_c<--0.1 #average variance (effect) of all chromosomes
sigma2_c<-0.1 #average variance of each individual chromosome 
sigma2_e<-0.5 #average residual var

ef<-rnorm(33, beta_c, sqrt(sigma2_c))#generating distributions based on average effects specified

#tot_ef<-matrix(0, 1000,1)

for(i in 1:33){
  tot_ef<-tot_ef+ef[i]*M[,i]
} #total effect 

#tot_ef<-M%*%ef## doing the same thing as above 

y<-rnorm(1000, tot_ef, sqrt(sigma_e))

dat<-data.frame(y=y, sumM=rowSums(M))
dat$M<-M #where M is the matrix of all chromosomal effects

m1<-MCMCglmm(y~sumM, random=~idv(M), data=dat)

## can see from summary that the estimates are pretty close to what we had simulated
# for average effect etc 
