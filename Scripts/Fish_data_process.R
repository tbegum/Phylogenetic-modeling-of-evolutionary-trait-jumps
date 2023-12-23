
library(hash)
library(gridExtra)


##Setting working directory to use functions
setwd("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Codes_Manuscipt3/")
folder1<-c("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Data/Collected_preprocessed_expression_data/PhyloFish_data/")
folder2<-c("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Data/Processed_expression_data/")
folder3<-c("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Data/")
folder4<-c("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Data/Collected_preprocessed_expression_data/NCBI/")
folder5<-c("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/")


######################### Step-1: Processing kallisto output data for fish #######################################

## Accessing the phylofish kallisto data to format and to use for our study
## Declaring the name of fish species data we need
## Each directory of fish species contains directories with number of tissue names and inside each directory there is one tsv file containing the tpm data for the corresponding tissue

## Reading the name of fish species as a vector
fish_species<-c("medaka","cavefish","northern_pike","tilapia")

## Collecting all the tissue names and their corresponding data
for(i in 1:length(fish_species))
{
  ## Declaring variables
  species_of_interest<-NULL
  tissue_names<- NULL
  
  ## Collecting the data
  species_of_interest<-fish_species[i] 
  if(species_of_interest %in% c("medaka","cavefish","northern_pike"))
  {
    tissue_name_files<-as.vector(list.files(paste0(folder1,"/",species_of_interest,"/", sep="")))
  }
  if(species_of_interest %in% "tilapia")
  {
    tissue_name_files<-as.vector(list.files(paste0(folder4,"/",species_of_interest,"/", sep="")))
  }
  dataframe<-NULL
  
  ## For each tissue, we collect the corresponding abundance.tsv file
  for(x in 1:length(tissue_name_files))
  {
    tissue<-NULL
    tissue_abundance<-NULL
    tissue<-tolower(tissue_name_files[x])
    original_tissue_name<-tissue_name_files[x]
    if(species_of_interest %in% c("medaka","cavefish","northern_pike"))
    {
      tissue_abundance<-read.table(paste0(folder1,"/",species_of_interest,"/",original_tissue_name,"/",original_tissue_name,"/","abundance.tsv",sep=""),header = T,sep="\t")
    }
    if(species_of_interest %in% "tilapia")
    {
      tissue_abundance<-read.table(paste0(folder4,"/",species_of_interest,"/",original_tissue_name,"/","abundance.tsv",sep=""),header = T,sep="\t")
    }
    tissue_abundance<-tissue_abundance[c(1,5)]
    
    ## formatting tissue_names for the species
    tissue[which(tissue %in% "head_kidney")]<- "kidney"
    tissue[which(tissue %in% "headkidney_1")]<- "kidney_1"
    tissue[which(tissue %in% "headkidney_2")]<- "kidney_2"
    tissue[which(tissue %in% "gills")]<- "lung"
    tissue[which(tissue %in% "niovary")]<- "ovary"
    tissue[which(tissue %in% "nitestis")]<- "testis"
    colnames(tissue_abundance)<-c("transcript",paste("TPM.",tissue,sep=""))
    
    ## Collecting the kallisto output files
    tissue_abundance$transcript<-as.character(tissue_abundance$transcript)
    if(x==1){dataframe<-data.frame(transcript=tissue_abundance$transcript)}
    dataframe<-merge(dataframe,tissue_abundance,by=c("transcript"))
  }
  
  ## Writing the result in a file
  filename=paste0(folder2,"/",species_of_interest,"_tpm_data.txt",sep="")
  write.table(dataframe,filename,sep = "\t",row.names = F,quote = F)
}


################### Step-2: Processing the kallisto.tpm output of fish ##############################
###Here we aim to obtain corresponding Ensembl gene ID for the kallisto proceessed species data 
### We take mean values between replicates of tissues for a species

## Reading the name of fish species as a vector
fish_species<-c("medaka","cavefish","northern_pike","tilapia")

## To match the corresponding ensembl gene name, we created a hash
fish_species_ensembl=hash(keys=c("medaka","cavefish","northern_pike","tilapia"),
                          values=c("medaka_ens.txt","cavefish_ens.txt","northern_pike_ens.txt","tilapia_ens.txt"))

## Running loop to collect all the tissue names and their corresponding data
for(i in 1:length(fish_species))
{
  ## Declaring variables
  species_of_interest<-NULL
  ensembl_gene_file<-NULL
  gene_file<-NULL
  
  ## Collecting the species ensemble gene name file and reading it
  species_of_interest<-fish_species[i]
  ensembl_file<-fish_species_ensembl[[species_of_interest]]
  gene_file<-read.table(paste0(folder3,ensembl_file,sep=""),header=T, sep="\t")
  colnames(gene_file)<-c("Ens.ID","transcript.new")
  
  ## Reading the kallisto processed file for the species
  kallisto_data<-read.table(paste0(folder5,species_of_interest,"_tpm_data.txt",sep=""),header=T, sep="\t")
  
  ## To remove the "." in the transcript column
  kallisto_data$transcript<-as.character(kallisto_data$transcript)
  kallisto_data$transcript.new<-lapply(kallisto_data$transcript,function(x){unlist(strsplit(x, "\\."))[[1]]})
  
  ## Rearranging the data
  Tissue_number<-length(kallisto_data)
  kallisto_data<-kallisto_data[c(Tissue_number,2:(Tissue_number-1))]
  
  ## Merging the datafiles to obtain corresponding ensembl gene id
  gene_file<-merge(gene_file,kallisto_data,by=c("transcript.new"))
  nTissue<-length(gene_file)
  gene_file<-gene_file[c(2:nTissue)]
  gene_file.new<-gene_file[order(gene_file$Ens.ID),] ## Sorting files based on the Ensembl gene ID
  
  ## Next we aim to sum the TPM values of tissues for the same Ensembl gene ID
  gene_file.mod<-aggregate(.~Ens.ID, gene_file.new, FUN= sum)
  
  ## We also need to calculate the rowmeans for tissues with multiple entries
  ## For this we need to find out the columns with multiple entries
  column_name<-colnames(gene_file.mod)
  multiple.tissue<-column_name[grep("_",column_name)]
  
  if(length(multiple.tissue)>0)
  {
    new.colname <- lapply(multiple.tissue,function(x){unlist(strsplit(toString(x), split='_', fixed=TRUE))[1]})
    tissuesNames<-unique(new.colname)
    
    ##Collecting mean.tissue.tpm data for replicates for each gene
    for(x in 1:length(tissuesNames))
    {
      tissue<-NULL
      tissue<-tissuesNames[[x]]
      gene_file.mod$tissue.name<-rowMeans(gene_file.mod[,regexpr(tissue,colnames(gene_file.mod))>0],na.rm = T)
      names(gene_file.mod)[names(gene_file.mod)=="tissue.name"]<- tissue
    }
    
    ## dropping the multiple entries column
    columns_to_remove<-grep("_",column_name)
    if(length(columns_to_remove)>0)
    {
      gene_file.mod<-gene_file.mod[-columns_to_remove]
    }
  }
  
  ## Writing the result in an output file
  outfile=paste0(folder2,"/",species_of_interest,".processed.tpm.txt",sep="")
  write.table(gene_file.mod,outfile, row.names = F, sep="\t",quote = F)  
}



