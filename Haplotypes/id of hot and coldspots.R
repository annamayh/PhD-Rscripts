### Finding windows containing hotspots ###

library(ggplot2)

setwd("H:/PHD_2ndYR/Haplotype_diversity")
pall <- read.table("old_files/Pall_trueCM_filtered.txt", header = T, stringsAsFactors = F)
head(pall)

quantile(pall$pall, c(.01, .99)) 


## hot/cold spots ###

pall26 <- subset(pall, CHR == "2")
ggplot(pall26, aes(BP, pall)) +
  geom_point() +
  labs(x = "SNPs", y = "P all", title = "Frequecy of ids with ROH at SNPs across chromosome 18" ) +
  geom_hline(yintercept = 0.1447017, linetype="dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0.05881768, colour = "red") 

coldspots<-subset(pall, pall == 0)
unique(coldspots$CHR)
head(coldspots)
coldspots3<-subset(coldspots, CHR == 30)
coldspots3<-subset(coldspots3, BP < 50365260 )
arrange(coldspots3, BP)
head(coldspots3)
tail(coldspots3)


load("Haplotype_div_chr30.RData")


which(colnames(phased_lab)=="cela1_red_12_56063030" )
which(colnames(phased_lab)=="cela1_red_12_56723354" )

which(colnames(phased_lab)=="cela1_red_12_77089550" )
which(colnames(phased_lab)=="cela1_red_12_77872818" )

which(colnames(phased_lab)=="cela1_red_12_83667124" )
which(colnames(phased_lab)=="cela1_red_12_84962518" )



################# hotspots ########################
hotspots<-subset(pall, pall > 0.1496576) ##0.1447017 = genetic distance top 1%  hotspots are still same when using this threshold
head(hotspots)
unique(hotspots$CHR)
hotspots3<-subset(hotspots, CHR == 30)
hotspots3<-subset(hotspots3, BP > 25110457 )
arrange(hotspots3, BP)
head(hotspots3)
tail(hotspots3)
max(hotspots3$pall)
mean(hotspots3$pall)
sd(hotspots3$pall)
