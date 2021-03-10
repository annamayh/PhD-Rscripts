##########################################
### using plink to search for ROH in R ###
##########################################

setwd("H:/") # make sure directory is the same location as the plink application #


## Plink_ROH function is a function to use plink to search for ROH with different input files(vcf/bfile)and output file names 
# 
## file path can be changed directly from the plink input file path or the input file name later

plink_ROH<-function(Input_file, Output_file){
  
  system(paste0("plink --vcf PHD_2ndYR/Slim/Output/Neutral_model_outputs/",Input_file," ",
                "--autosome --autosome-num 33 --out PHD_2ndYR/Slim/Output/ROH_output/",Output_file," ",
                "--homozyg --homozyg-window-snp 35 --homozyg-snp 40 --homozyg-kb 2500 ",
                "--homozyg-density 70 --homozyg-window-missing 4 ",
                "--homozyg-het 0 ",
                "--maf 0.01 --freq --missing"))
  
  
}


## loop function over multiple files ##

for (i in 1:10){
  
  vcf<-paste0("test_neutral_loop",i) ## all files in model output are called test_neutral_loop1 etc 
  output<-paste0("testing_output_loop",i) ## want output files to be numbered too 
  
  plink_ROH(vcf,output) #run function for all files labeled 1:10
  
  
  
  
}
