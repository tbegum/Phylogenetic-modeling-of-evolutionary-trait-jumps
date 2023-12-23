######################### Step-1: Processing data from BGEE for vertebrates #######################################
## Accessing BGEE with R package to obtain datasets for multiple species
## To remove biasness in comparison we took tissues of same developmental stages (post juvenile adult stage)
library(BgeeDB)
library(hash)
library(preprocessCore)

#browseVignettes("BgeeDB")

setwd("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Codes_Manuscipt3/")
#dir.create("Expression_data_BGEE")
folder1<-c("~/Desktop/nwork/Empirical/protein_trees.nhx/Expression_data_BGEE/")
folder2<-c("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Data/Processed_expression_data/")

## To keep only adult stage data, we used the devstage ID provided by Juliyen
devstageID<-read.table("devStageIds_BGEE.txt",header = T)

## Creating hashes to generate expression data tables of 7 vertebrate species

speciess<- c("mouse","zebrafish","chicken","rat","macaque","cow","human","opossum")
species_name=hash(keys=c("mouse","zebrafish","chicken","rat","macaque","cow","human","opossum"),
                  values=c("10090","7955","9031","10116","9544","9913","9606","13616"))


species_experimental_id=hash(keys=c("mouse","zebrafish","chicken","rat","macaque","cow","human","opossum"),
                             values=c("GSE41637","SRP044781","GSE41637","GSE41637","GSE41637","GSE41637","GSE30611","GSE30352"))

## Creating data table of TPM expression data from BGEE v.14
for(i in 1:length(species_name))
{
  ## Reading the species name
  species<-NULL
  species<- speciess[i]
  
  ## Accessing bgee RNA-seq data for the corresponding species 
  bgee.RNAseq.data<-Bgee$new(species=species_name[[species]], dataType="rna_seq")
  annotation.bgee<-getAnnotation(bgee.RNAseq.data)
  
  #lapply(annotation.bgee, head)
  ## Obtaining experimentalID from the using the 
  #experimentId<-annotation.bgee$experiment.annotation$Experiment.ID[which(annotation.bgee$experiment.annotation$Organ.count==max(annotation.bgee$experiment.annotation$Organ.count))]
  experimentId<-species_experimental_id[[species]]
  #stage<-annotation.bgee$experiment.annotation$Organ.stage.count[which(annotation.bgee$experiment.annotation$Organ.count==max(annotation.bgee$experiment.annotation$Organ.count))]
 
  ## Obtaining TPM data for the corresponding experimentalID
  bgee_data <- getData(bgee.RNAseq.data, experimentId = experimentId)
  
  ## Checking for adult stages (omit if we wants data for other stages like embryonic stage)
  bgee_data_new<-bgee_data[which(bgee_data$Stage.ID %in% devstageID$stageId),]
  
  ## Choosing only the specific columns we need for our work
  bgee_data_neww <- bgee_data_new[c(4,6,12)]

  ## Extracting the TPM data for the corresponding tissues/organs for each genes
  expression_dataset<- reshape(bgee_data_neww, idvar = "Gene.ID", timevar = "Anatomical.entity.name", v.names = "TPM", direction = "wide")
  names(expression_dataset)[names(expression_dataset) =="TPM.skeletal muscle tissue"]<-"TPM.muscle"
  names(expression_dataset)[names(expression_dataset) =="TPM.adult mammalian kidney"]<-"TPM.kidney"
  names(expression_dataset)[names(expression_dataset) =="TPM.muscle tissue"]<-"TPM.muscle"
  names(expression_dataset)[names(expression_dataset) =="TPM.pharyngeal gill"]<-"TPM.gill"
  names(expression_dataset)[names(expression_dataset) =="TPM.female gonad"]<-"TPM.ovary"
  names(expression_dataset)[names(expression_dataset) =="TPM.mature ovarian follicle"]<-"TPM.unfertilized.egg"
  names(expression_dataset)[names(expression_dataset) =="TPM.bone tissue"]<-"TPM.bones"
  names(expression_dataset)[names(expression_dataset) =="Gene.ID"]<-"Ens.ID"
  
  if((species=="zebrafish"))
  {
     expression_dataset<-expression_dataset[c(1,12,7,5,3,10,6,2,9,11,8,4)]
  }
  
  ## Writing the datatables inside of a folder
  filename<-paste0(folder1,"TPM_data_",species,"_",experimentId,".txt",sep="")
  filename2<-paste0(folder2,"/TPM_data_",species,"_",experimentId,".txt",sep="")
  write.table(expression_dataset,filename, sep = "\t",quote = F,row.names = F)
  write.table(expression_dataset,filename2, sep = "\t",quote = F,row.names = F)
}




