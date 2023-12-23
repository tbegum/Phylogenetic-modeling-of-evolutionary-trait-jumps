## This script is written to process expression data collected from Chen et al. 2019

library(hash)
library(gridExtra)


##Setting working directory to use functions
setwd("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Codes_Manuscipt3/")
folder1<-c("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Data/Collected_preprocessed_expression_data/Chen_et_al_2019/")
folder2<-c("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Data/Processed_expression_data/")



######################### Processing output of Aviv Regev for Dog, Rabbit, and Ferret #######################################
## The processed data was collected using the perl program "sra2sh.pl" 
## Each directory of the species contains directories with number of tissue names and inside each directory there is one tsv file containing the tpm data for the corresponding tissue

## Reading the name of tetrapod species as a vector
tetrapod_species<-c("dog","rabbit","ferret")

## Collecting all the tissue names and their corresponding data
for(i in 1:length(tetrapod_species))
{
  ## Declaring variables
  species_of_interest<-NULL
  tissue_names<- NULL
  
  ## Collecting data
  species_of_interest<-tetrapod_species[i] 
  tissue_name_files<-as.vector(list.files(paste0(folder1,"/",species_of_interest,"/", sep="")))
  
  dataframe<-NULL
  
  ## For each tissue collecting tissue wise expression data the corresponding abundance.tsv file
  for(x in 1:length(tissue_name_files))
  {
    tissue<-NULL
    tissue_abundance<-NULL
    tissue<-tolower(tissue_name_files[x])
    original_tissue_name<-tissue_name_files[x]
    filename<-NULL
    filename<-list.files(paste0(folder1,species_of_interest,"/",original_tissue_name,"/",sep=""))
    tissue_abundance<-read.table(paste0(folder1,species_of_interest,"/",original_tissue_name,"/",filename,sep=""),header = T,sep="\t")
    tissue_abundance<-tissue_abundance[c(1,6)]
    colnames(tissue_abundance)<-c("Ens.ID",paste("TPM.",tissue,sep=""))
    
    ## Collecting the output files
    tissue_abundance$Ens.ID<-as.character(tissue_abundance$Ens.ID)
    if(x==1){dataframe<-data.frame(Ens.ID=tissue_abundance$Ens.ID)}
    dataframe<-merge(dataframe,tissue_abundance,by=c("Ens.ID"))
  }
  
  ## Writing in files
  filename=paste0(folder2,"/",species_of_interest,".tpm.data.txt",sep="")
  write.table(dataframe,filename,sep = "\t",row.names = F,quote = F)
}


  
