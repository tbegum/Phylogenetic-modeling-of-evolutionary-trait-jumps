## This script is written to calculate tissue-specificity and average gene expression levels for 15 vertebrates using processing expression data 

library(hash)
library(gridExtra)

##Setting working directory to use functions
setwd("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Codes_Manuscipt3/")
folder1<-c("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Data/Processed_expression_data/")
folder2<-c("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Data/Tau_and_Average_expression/")

### Required functions
## Computing average gene expression level and tissue specificity for our data
## For this we used 3 functions created by Nadezda (KMRR 2016)

#Mean value per gene is calculated
fmean <- function(x)
{
  if(!all(is.na(x))) {
    res <- mean(x, na.rm=TRUE)
  } else {
    res <- NA
  }
  return(res)
}
###***###***###	


##Function require a vector with expression of one gene in different tissues.
#Max value per gene is calculated.
fmax <- function(x)
{
  if(!all(is.na(x))) {
    res <- max(x, na.rm=TRUE)
  } else {
    res <- NA
  }
  return(res)
}


#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
#Minimum 2 tissues
fTau <- function(x)
{
  if(all(!is.na(x))) {
    if(min(x, na.rm=TRUE) >= 0) {
      if(max(x)!=0) {
        x <- (1-(x/max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res/(length(x)-1)
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    } 
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  } 
  return(res)
}

######################### Computing tau and average gene expression levels from processed TPM data for 15 species #######################

species_of_interest<-c("northern_pike","cavefish","medaka","zebrafish","tilapia",
                       "human","mouse","rat","cow","macaque","chicken","dog","rabbit","ferret","opossum")
species_hash_key=hash(keys=c("northern_pike","cavefish","medaka","zebrafish","tilapia",
                             "human","mouse","cow","rat","macaque","chicken","dog","rabbit","ferret","opossum"),
                           values=c("northern_pike.processed.tpm.txt","cavefish.processed.tpm.txt",
                                    "medaka.processed.tpm.txt","TPM_data_zebrafish_SRP044781.txt","tilapia.processed.tpm.txt",
                                    "TPM_data_human_GSE30611.txt","TPM_data_mouse_GSE41637.txt","TPM_data_cow_GSE41637.txt","TPM_data_rat_GSE41637.txt","TPM_data_macaque_GSE41637.txt",
                                    "TPM_data_chicken_GSE41637.txt","dog.tpm.data.txt","rabbit.tpm.data.txt","ferret.tpm.data.txt",
                                    "opossum.tpm.data.txt"))

common_organs<-c("TPM.brain","TPM.heart","TPM.kidney","TPM.liver","TPM.muscle","TPM.testis")

for(i in 1:length(species_of_interest))
{
  print (i)
  ## Declaring variables
  sspecies_of_interest<-NULL
  gene_file<-NULL
  raw<-NULL
  y<-NULL
  organ_data<-NULL
  
  sspecies_of_interest<-species_of_interest[i]
  
  
   gene_file<-species_hash_key[[sspecies_of_interest]]
    ## Reading file and log transforming
    raw<-read.table(paste0(folder1,gene_file,sep=""),header=T, sep="\t")
    log.normalized.data<-raw
    x <- log.normalized.data[,c(-1)]
    x[x < 1] <- 1    ## for log transformation we changed values less than 1 to 1
    log.normalized.data[,c(-1)] <- log2(x)
    #summary(log.normalized.data)
  
    ## removing non-expressed genes based on sum of expression > 0 across tissue
    y<-log.normalized.data[,2:ncol(log.normalized.data)]
    log.normalized.data$sum<-apply(y, 1, sum) ## replaced mean with sum
    log.normalized.data2<-log.normalized.data[which(log.normalized.data$sum>0),] 
    #summary(log.normalized.data2)
    
    organ_data<-data.frame(Ens.ID=log.normalized.data2[,1])
    organ_data<-cbind(organ_data,log.normalized.data2[common_organs])
    
    ## calculating Tau, Mean, and Max
    organ_data$Tau<-apply(organ_data[,c(-1)], 1, fTau)
    organ_data$Mean<-apply(organ_data[,c(-1)], 1, fmean)
    organ_data$Max<-apply(organ_data[,c(-1)], 1, fmax)
    
   
  ## Writing the result in an output file
  outfile=paste0(folder2,"/",sspecies_of_interest,"_tau.txt",sep="")
  write.table(organ_data,outfile, row.names = F, sep="\t",quote = F)
} 
  
  
