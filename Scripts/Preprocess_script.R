
########## Main: This code is written to analyze expression evolution of vertebrates ################

########## Chunk1: Loaded required libraries, setting path and calling required functions ################
library(ape)
library(geiger)
library(phytools)
library(treeio)
library(stringr)
library(ggrepel)
library(dplyr)
library(ggtree)
library(treeio)
library(magrittr) 
library(parallel) ##Needed to run read.nhx from treeio 
library(digest)
library(ggplot2)
library(gtools)
library(caper)
#library(cowplot)
#devtools::install_github( "kassambara/easyGgplot2", quiet=TRUE )
#library(easyGgplot2)
library(ggplot2)
library(gridExtra)
library(R.utils)
library(hash)
#library(edgeR)

##Setting working directory to use functions
setwd("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/")
folder1 <- c("~/Desktop/nwork/Empirical/Manuscript_3/Data/Normalized_expression_data/")
folder2 <- c("~/Desktop/nwork/Empirical/Manuscript_3/Plot_Manuscipt3/")
folder3 <- c("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Data/Tau_and_Average_expression/")

## The functions written to test empirical data
source("functions_TPCM.R")

## Setting seed number for reproducibility 
set.seed(786)

## To perform parallel computation we need to detect core numbers
cores <- detectCores()

########## Chunk2: Processing Empirical trees from Genomicus v.95.01 for further analyses ############

## unzipping and parsing the phylogenetic tree from Genomicus v. 95.01 database
R.utils::bunzip2("protein_trees.nhx.bz2","protein_trees_95.01.nhx",remove = FALSE, skip = TRUE)
proteintree <- "protein_trees_95.01.nhx"
protree <- treeio::read.tree(proteintree) ##also treeio::read.tree

## No of protein trees total protein trees present in the database
Initial.trees.database <- length(protree)

## Reading each line of the tree file i.e. one tree per line
tree.lines <- readLines("protein_trees_95.01.nhx", warn = FALSE)

## Parsing each tree data with added events according to the node numbers
## This step takes some time
trees.all<-mclapply(tree.lines, Parse.trees.Genomicus.95, mc.cores=cores) ## started with 43491 trees
## save.image("Event_annotated_trees.rda")

########## Chunk3: Initial annotation of trees with Tau and with average expression levels ############
## The  5 fish species and 10 outgroup tetrapoda we kept are - (1) zebrafish (^ENSDARG), (2) northern pike (^ENSELUG), 
## (3) cavefish (^ENSAMXG), (4) medaka (^ENSORLG), (5) tilapia (^ENSONIG),(6) human (^ENSG), (7) mouse (^ENSMUSG),
## (8) rat (^ENSRNOG), (9) cow (^ENSBTAG), (10) macaque (^ENSMMUG), (11) chicken (^ENSGALG), (12) dog (^ENSCAFG),
## (13) rabbit (^ENSOCUG), (14) ferret (^ENSMPUG), and (15) opossum (^ENSMODG)

# Files with tau and average gene expression data of 15 species in 6 different tissues ("brain", "heart","kidney,"liver","muscle","testis")
# We used "/Users/admin/Desktop/nwork/Empirical/Manuscript_3/kallisto_output/tau.R" to generate tau and average gene expression on log2 normalized expression data of all species 
raw.expression.datafiles <- c(paste(folder3,"mouse_tau.txt",sep=""), 
                              paste(folder3,"human_tau.txt",sep=""),
                              paste(folder3,"rat_tau.txt",sep=""),
                              paste(folder3,"cow_tau.txt",sep=""),
                              paste(folder3,"macaque_tau.txt",sep=""),
                              paste(folder3,"chicken_tau.txt",sep=""),
                              paste(folder3,"dog_tau.txt",sep=""),
                              paste(folder3,"ferret_tau.txt",sep=""),
                              paste(folder3,"opossum_tau.txt",sep=""),
                              paste(folder3,"rabbit_tau.txt",sep=""),
                              paste(folder3,"zebrafish_tau.txt",sep=""),
                              paste(folder3,"medaka_tau.txt", sep=""),
                              paste(folder3,"cavefish_tau.txt",sep=""),
                              paste(folder3,"northern_pike_tau.txt",sep=""),
                              paste(folder3,"tilapia_tau.txt",sep=""))

## Reading trait data of all species in a single data frame
raw.exp.data<- bind_rows(lapply(raw.expression.datafiles, read.table, stringsAsFactors=FALSE, sep="\t", header=TRUE))
raw.exp.data<-raw.exp.data[c(1,8,9)]
colnames(raw.exp.data)[c(1,3)]<-c("label","Raw.avg")

## Adding exp data to our trees
exp.annotated.trees<- mclapply(trees.all,
                               function(tree){
                                 tree@data %<>% left_join(raw.exp.data, by = c("label")) 
                                 return(tree)}, 
                               mc.cores=cores)

################ Chunk4: Identifying phylogenies with atleast one speciation and one duplication events ################
## Now we need to prune trees with tips matching species name for which expression data are available
## From PhyloFish we have only 5 species common between Genomicus & PhyloFish & NCBI (We excluded atlantic cod because of bad Tau distribution plot)
## So we need to identify trees having exp data at the tips
## Annotating tree tips with trait data
pruned.annotated.trees <- mclapply(exp.annotated.trees, drop.leaf, mc.cores = cores)
pruned.annotated.trees <- pruned.annotated.trees[!is.na(pruned.annotated.trees)] ## Removing pruned trees with NA data; passed 17103 trees
initial.index.info2<-bind_rows(lapply(pruned.annotated.trees, Info_collection_tree))
initial.index.info2$tree_num<-rownames(initial.index.info2)

## Trees with all duplication events
trees.all.duplication<-mclapply(pruned.annotated.trees, all.duplication, mc.cores = cores)
trees.all.duplication <- trees.all.duplication[!is.na(trees.all.duplication)] ## only 148 trees passed

## We need to proceed with trees only having atleast one speciation and duplication events
## Identifying trees with atleast one duplication and one speciation events
trees.all.pruned<-mclapply(pruned.annotated.trees, filter.trees, mc.cores = cores)
trees.all.pruned <- trees.all.pruned[!is.na(trees.all.pruned)] ## 9039 trees passed the criteria

## Trees with all speciation events
trees.all.speciation<-length(pruned.annotated.trees) - (length(trees.all.duplication) + length(trees.all.pruned))
trees.all.spec<-mclapply(pruned.annotated.trees, all.speciation, mc.cores = cores)
trees.all.speciation.new <- trees.all.spec[!is.na(trees.all.spec)]## 7916 trees with tip number >= 4
speciation.index.info<-bind_rows(lapply(trees.all.speciation.new, Info_collection_tree))
speciation.index.info$tree_num<-rownames(speciation.index.info)
speciation.15species.treeno<-as.numeric(speciation.index.info$tree_num[which(speciation.index.info$tip_num==15)]) 
trees.all15.speciation<-trees.all.speciation[speciation.15species.treeno] ## 1077 trees with 15 species (1:1 ortholog)
rm(trees.all.spec)

################ Chunk5: Using the species trees to generate 1:1 orthologs to normalize the expression data ################
freq<-0
ortholog.set<-bind_rows(lapply(trees.all.speciation.new, ortholog.finder)) ## it's now 1077 genes with 1:1 ortholog

################ Chunk6a: Generating expression matrix for 6 common tissues and expression normalization following the method of Brawand et al. (2011) #############################
## To call the function, we need to provide:
# 1. dataframe of orthologs or orthogroups,
# 2. type = "ortholog" or "orthogroup"

pca.ortholog<-expression.all.sample(ortholog.set,"ortholog")
exp.normalized.ortholog<-NULL
exp.normalized.ortholog<-normalize.explevel.new(pca.ortholog)


################ Chunk7: PCA on orthologs  #############################
## To call the PCA_plot function, we need to provide:
# 1. normalized data matrix,
# 2. type = "ortholog" or "common.WGD", 
# 3. group.by = "tissue"/"species"
# 4. method = "normal.2D"

## Generating PCA plot to assess the result of the analysis
#par(mfrow=c(1,2))
PCA.plot(pca.ortholog,type="ortholog",group.by="tissue",method="normal.2D") ## before normalization
PCA.plot(exp.normalized.ortholog,type="ortholog",group.by="tissue",method="normal.2D") ## after normalization

################ Chunk 8: Generating normalized expression data for each species and annotating trees with normalized expression levels #############################
normalized.exp.data<- normalized.species.data(paste0(folder1,"/OG.scaling_factor.txt",sep=""))
normalized.exp.data$nTau <- apply(normalized.exp.data[c(2:7)],1,fTau)
normalized.exp.data$avg.exp <- apply(normalized.exp.data[c(2:7)],1,fmean)
normalized.exp.data$max.exp <- apply(normalized.exp.data[c(2:7)],1,fmax)
normalized.exp.data$sum.exp <- apply(normalized.exp.data[c(2:7)],1,sum)
colnames(normalized.exp.data)[1]<-c("label")

## Considering all genes with sum of expression level across all tissues are equal to or greater than 0
## This process filters genes which are no expressed in any of the six tissues
normalized.exp.data.mod <- normalized.exp.data[normalized.exp.data$sum.exp >= 0,]

## Adding normalized exp data our pruned trees
norm.exp.annotated.pruned.trees<- mclapply(trees.all.pruned,
                                           function(tree){
                                             tree@data %<>% left_join(normalized.exp.data.mod, by = c("label")) 
                                             return(tree)}, 
                                           mc.cores=cores)

norm.exp.annotated.species.trees<- mclapply(trees.all.speciation.new,
                                            function(tree){
                                              tree@data %<>% left_join(normalized.exp.data.mod, by = c("label")) 
                                              return(tree)}, 
                                            mc.cores=cores)

################ Chunk 9: Adding node height and node depth to the data, also modifying label of duplicates #############################
## This function aims to add branch length, node depth and node height required for building calibration matrix
trees.pruned.modified<-mclapply(norm.exp.annotated.pruned.trees, tree.nodedepth, mc.cores = cores) 
trees.pruned.modified<-mclapply(trees.pruned.modified, tree.height, mc.cores = cores) #9039 trees

trees.speciation.modified<-mclapply(norm.exp.annotated.species.trees, tree.nodedepth, mc.cores = cores)
trees.speciation.modified<-mclapply(trees.speciation.modified, tree.height, mc.cores = cores) #7916 trees

## Modifying gene tree clade label to distinguish speciation and duplication nodes
## This modification is necessary for time calibration of gene trees
trees.pruned2.modified<-mclapply(trees.pruned.modified, modify.label, mc.cores = cores) # 9039 trees with 4 tips

## Clearing memory
rm(raw.exp.data)
rm(exp.annotated.trees)
rm(trees.all)
rm(trees.fifteen.tips)

################ Chunk 10: Time calibration of trees and adding node age #############################

## Speciation time point
Calibration.time <- data.frame(
  label = c("Otophysi","Clupeocephala","Neopterygii","Acanthomorphata","Bilateria", "Vertebrata",
            "Chordata","Euteleostomi","Opisthokonta","Atherinomorphae","Ovalentaria","Percomorphaceae","Euteleosteomorpha","Pseudocrenilabrinae",
            "Osteoglossocephalai","Characoidei","Characiphysae","Euarchontoglires","Catarrhini","Murinae",
            "Mammalia","Amniota","Boreoeutheria","Primates","Homininae","Hominoidea","Rodentia","Glires","Theria","Caniformia","Lausasiatheria"),
  Mya = c(152.6, 229.9, 314.7, 148, 796.6, 615,
          676.4,435.3,1105,93.2,119.4,128,206.3,58.2,
          267.1,110,150,92,29,25,
          176.9, 311.9, 96.5, 73.8,9.1,20.2,70.5,82.1,162,39.9,78.5))
Calibration.time <- Calibration.time[order(-Calibration.time$Mya),]

## Time Calibrating all trees according to the our calibration time points
Time.calibrated.trees.pruned <- lapply(trees.pruned2.modified, 
                                       tree.calibrate,
                                       timeframe=Calibration.time,
                                       model="correlated")
rm(trees.pruned2.modified)

## Collecting all time calibrated trees 
Calibrated.all.pruned <-Time.calibrated.trees.pruned [!is.na(Time.calibrated.trees.pruned)] 
Calibrated.all.pruned2<-Calibrated.all.pruned[!sapply(Calibrated.all.pruned, is.null)]## 9022 trees passed with 4 tips
rm(Time.calibrated.trees.pruned)

## Time Calibrating 1634 speciation trees according to the our calibration time points
Time.calibrated.speciation.trees <- lapply(trees.speciation.modified, 
                                           tree.calibrate,
                                           timeframe=Calibration.time,
                                           model="correlated") #otherwise use "discrete" model
## Collecting all time calibrated trees 
Calibrated.speciation.trees <-Time.calibrated.speciation.trees[!is.na(Time.calibrated.speciation.trees)]
Calibrated.speciation.trees2<-Calibrated.speciation.trees[!sapply(Calibrated.speciation.trees, is.null)]## passed 4871 trees with 10 tips and 6515 trees with 4 tips minimum
rm(Time.calibrated.speciation.trees)
rm(Calibrated.speciation.trees)

## Remodifying gene tree clade label to initial label
## This modification is necessary for pic calculation as well as node age storage
Calibrated.all.modify<-mclapply(Calibrated.all.pruned2, remodify.label, mc.cores = cores)
Calibrated.speciation.trees.modify<-mclapply(Calibrated.speciation.trees2, remodify.label, mc.cores = cores)

## Adding nodeage to our data
Calibrated.all.modify<-mclapply(Calibrated.all.modify, node.age, mc.cores = cores)
Calibrated.speciation.trees.modify<-mclapply(Calibrated.speciation.trees.modify, node.age, mc.cores = cores)

#save.image("norm_exp_15spe_all_tree.rda")

######### Chunk 11: Filtering trees before analysis #################

## Dropping leaf with no average exp
trees.all.interestt<-mclapply(Calibrated.all.modify, drop.leaf.no.exp, mc.cores = cores)
trees.all.interestt<-trees.all.interestt[!is.na(trees.all.interestt)] ## 9022 trees passed with 4 tips trees since we considered average expression >=0

## However many trees can just contain outgroup species, not fish species of interest
## So, we need to keep trees for which at least 3 tips with fish data available, and atleast 1 tip should contain outgroup tetrapod species
trees.all.interest1<-mclapply(trees.all.interestt, filter1, mc.cores = cores)
trees.all.interest<-trees.all.interest1[!is.na(trees.all.interest1)]## 6923 trees

##Dropping leaf with no average exp for speciation tree
trees.speciation.interest<-mclapply(Calibrated.speciation.trees.modify, drop.leaf.no.exp, mc.cores = cores)
trees.speciation.interest<-trees.speciation.interest[!is.na(trees.speciation.interest)] ## 6515 trees passsed for 4 tips

#save.image("more_outgroup_added_exp_15spe_4tips.rda")

################### Chunk 12 indexing trees by names ####################
## Here we aim to write a function that provides an index to all trees of interest
## For this, we considered processed tree list of 6923 trees
Index<-0
for(Index in 1:length(trees.all.interest))
{
  trees.all.interest[[Index]]@phylo$index.tree<-Index ## Now all tree will have index with them
}
rm(Index)

####################### Chunk 13: Analysis on contrasts standardized trees for traits ###############################
# Diagnostic tests for Tau
count<-0
Calibrated.all.modify.standardized.tau1n<-lapply(trees.all.interest, diagnostic.plot.test, 0)
Calibrated.all.modify.standardized.tau2n<-Calibrated.all.modify.standardized.tau1n[!is.na(Calibrated.all.modify.standardized.tau1n)] 
Calibrated.all.standardized.tau<-Calibrated.all.modify.standardized.tau2n[!sapply(Calibrated.all.modify.standardized.tau2n, is.null)]## finally we obtained 4247 tree data passing the diagnostic test
rm(Calibrated.all.modify.standardized.tau1n)
rm(Calibrated.all.modify.standardized.tau2n)

# Diagnostic tests for Exp level
count<-0
Calibrated.all.modify.standardized.exp1<-lapply(trees.all.interest, diagnostic.plot.test.exp, 0)
Calibrated.all.modify.standardized.exp2<-Calibrated.all.modify.standardized.exp1[!is.na(Calibrated.all.modify.standardized.exp1)] 
Calibrated.all.standardized.exp<-Calibrated.all.modify.standardized.exp2[!sapply(Calibrated.all.modify.standardized.exp2, is.null)] ## 4603 trees passed the diagnostic tests
rm(Calibrated.all.modify.standardized.exp1)
rm(Calibrated.all.modify.standardized.exp2)

#### Analyses for speciation trees

# Diagnostic tests for Tau
count<-0
Calibrated.speciation.trees.modify.standardized.tau1<-lapply(trees.speciation.interest, diagnostic.plot.test, 0)
Calibrated.speciation.trees.modify.standardized.tau2<-Calibrated.all.modify.standardized.tau1[!is.na(Calibrated.all.modify.standardized.tau1)] 
Calibrated.speciation.trees.modify.standardized.tau<-Calibrated.all.modify.standardized.tau2[!sapply(Calibrated.all.modify.standardized.tau2, is.null)]##  3796 trees passed the diagnostic tests
rm(Calibrated.speciation.trees.modify.standardized.tau1)
rm(Calibrated.speciation.trees.modify.standardized.tau2)

# Diagnostic tests for Exp level
count<-0
Calibrated.speciation.trees.modify.standardized.exp1<-lapply(trees.speciation.interest, diagnostic.plot.test.exp, 0)
Calibrated.speciation.trees.modify.standardized.exp2<-Calibrated.speciation.trees.modify.standardized.exp1[!is.na(Calibrated.speciation.trees.modify.standardized.exp1)] 
Calibrated.speciation.trees.modify.standardized.exp<-Calibrated.speciation.trees.modify.standardized.exp2[!sapply(Calibrated.speciation.trees.modify.standardized.exp2, is.null)]## 5176 trees passed the diagnostic tests
rm(Calibrated.speciation.trees.modify.standardized.exp1)
rm(Calibrated.speciation.trees.modify.standardized.exp2)

#save.image("norm_exp_15spe_4tips.rda") 

################### Chunk- 14: synteny validated sure 3RWGD trees identification and fish WGD nodes annotation #############

## To download synteny validated ohnolog pairs from ohnologs database for 2 (zebrafish, medaka) out of 5 (zebrafish, medaka), cavefish, pike, and tilapia) teleosts used for our study
## Since we can validate our phylogeny based identification by ohnologs database, we relaxed our phylogeny based criteria for these two species, where we allow SSDs on recent clades following fishWGD duplication
## Zebrafish 3R ohnolog downloaded from Ohnologs database (http://ohnologs.curie.fr/cgi-bin/BrowsePage.cgi)

### Step1: Extracting synteny validated ohnolog pairs for zebrafish and medaka from ohnologs database 

##We used relaxed criteria (Relaxed : q-score(outgroups) < 0.05 AND q-score(self comparison) < 0.3) of ohnolog database to obtain maximal matches
### For zebrafish 
WGDpair.ohnolog.zfish<-read.table("drerio.Pairs.Relaxed.3R.txt",header = T, sep="\t") ## 3311 pairs of 3RWGD

## Only duplication time "FishWGD" and "protein_coding" genes were considered
WGDpair.ohnolog.zfish.new<-WGDpair.ohnolog.zfish[WGDpair.ohnolog.zfish$Duplication.time=="FishWGD" & WGDpair.ohnolog.zfish$Gene.type=="protein_coding",] ## Total 3080 pairs found

rm(WGDpair.ohnolog.zfish)

### For medaka
WGDpair.ohnolog.medaka<-read.table("olatipes.Pairs.Relaxed.3R.txt",header = T, sep="\t") ## 2160 pairs of 3RWGD

## Only duplication time "FishWGD" was considered
WGDpair.ohnolog.medaka.new<-WGDpair.ohnolog.medaka[WGDpair.ohnolog.medaka$Duplication.time=="FishWGD"& WGDpair.ohnolog.medaka$Gene.type=="protein_coding",] ## Total 2099 pairs found

rm(WGDpair.ohnolog.medaka)

## Ohnolog data for both the zebrafish and medaka
WGDpair.ohnolog <- NULL
WGDpair.ohnolog <- bind_rows(WGDpair.ohnolog.zfish.new,WGDpair.ohnolog.medaka.new)

## Identifying "sure" 3RWGD trees for zebrafish and medaka where ohnologs are verified by ohnologs database
Count<-0
Validated.WGD<-NULL
Validated.WGD<-lapply(trees.all.interest,validate.WGD.2species,WGDpair.ohnolog)
Validated.WGD.new1<-Validated.WGD[!is.na(Validated.WGD)]  
Validated.WGD.new<-Validated.WGD.new1[!sapply(Validated.WGD.new1, is.null)] ##2050 trees
rm(Validated.WGD)
rm(Validated.WGD.new1)

##Identifying "whole genome duplication" node for zebrafish and medaka considering ohnolog database 
count<-0
Validated.WGD.2species1<-lapply(Validated.WGD.new,mrca.WGD,WGDpair.ohnolog)
Validated.WGD.2species<-Validated.WGD.2species1[!is.na(Validated.WGD.2species1)]
Validated.WGD.2species<-Validated.WGD.2species[!sapply(Validated.WGD.2species, is.null)] ##1788 trees (sure + unsure)
rm(Validated.WGD.2species1)

### Validated "sure" WGD trees for these 2 species 
## Next we sought to find all validated WGD trees between zebrafish and medaka validated from ohnolog data base
summary.validated.trees<-summary.tree(Validated.WGD.2species)

## Identifying validated WGD trees for 2 species
list.index<-as.vector(unique(summary.validated.trees$index.tree)) ##1788

## Identifying sure teleosts specific WGD trees (considering only "Clupeocephala" and "Osteoglossocephalai" clade duplications) considering ohnolog database 
count<-0
Validated.WGD.2species.sure<-lapply(Validated.WGD.new,mrca.WGD.types,WGDpair.ohnolog,type="sure")
Sure.WGD.trees<-Validated.WGD.2species.sure[!is.na(Validated.WGD.2species.sure)]
Sure.WGD.trees<-Sure.WGD.trees[!sapply(Sure.WGD.trees, is.null)] ##1159 trees
rm(Validated.WGD.2species.sure)

summary.sure.WGD.trees<-summary.tree(Sure.WGD.trees)
list.sure.3RWGD<-as.vector(unique(summary.sure.WGD.trees$index.tree))

#save.image("norm_exp_15spe_4tips_new.RData")

################### Chunk- 15: Synteny validated unsure 3RWGDs trees identification #############

## Unsure 3RWGD trees (i.e., we identified sure 3RWGD trees using ohnolog database with Zebarish and Medaka data)
## These paralog pairs are present in the ohnolog database but their duplication time does not match with "FishWGD"
## This includes phylogenies with "FishWGD" at duplication nodes other than "Clupeocephala" and "Osteoglossocephalai" 
## + phylogenies with paralog pairs present in ohnolog database but do not show annotated duplications as "FishWGD" 


### For zebrafish 
WGDpair.ohnolog.zfish<-read.table("drerio.Pairs.Relaxed.3R.txt",header = T, sep="\t") ## 3311 pairs of 3RWGD

##only protein coding paralogs were considered
WGDpair.ohnolog.zfish.type2<-WGDpair.ohnolog.zfish[(WGDpair.ohnolog.zfish$Gene.type=="protein_coding"),] ## Total 201 pairs found
rm(WGDpair.ohnolog.zfish)

### For medaka 
WGDpair.ohnolog.medaka<-read.table("olatipes.Pairs.Relaxed.3R.txt",header = T, sep="\t") ## 2160 pairs of 3RWGD

##Only protein coding paralogs were considered
WGDpair.ohnolog.medaka.type2<-WGDpair.ohnolog.medaka[(WGDpair.ohnolog.medaka$Gene.type=="protein_coding"),] ## Total 44 pairs found
rm(WGDpair.ohnolog.medaka)

## Cadidate ohnologs for both the zebrafish and medaka
WGDpair.cadidate.ohnolog <- NULL
WGDpair.candidate.ohnolog <- bind_rows(WGDpair.ohnolog.zfish.type2,WGDpair.ohnolog.medaka.type2)

## Identifying 3RWGD nodes in phylogenies 
count<-0
Validated.WGD.2species2<-lapply(Validated.WGD.new,mrca.WGD,WGDpair.candidate.ohnolog) 
Validated.WGD.2species.new<-Validated.WGD.2species2[!is.na(Validated.WGD.2species2)]
Validated.WGD.2species.new<-Validated.WGD.2species.new[!sapply(Validated.WGD.2species.new, is.null)] ##1795  trees (sure + unsure)
rm(Validated.WGD.2species2)

## Identifying unsure teleosts specific WGD trees (excluding "Clupeocephala" and "Osteoglossocephalai" duplication) following ohnolog database 
count<-0
Validated.WGD.2species.unsure<-lapply(Validated.WGD.2species.new,mrca.WGD.types,WGDpair.ohnolog,type="unsure")
Unsure.WGD.trees<-Validated.WGD.2species.unsure[!is.na(Validated.WGD.2species.unsure)]
Unsure.WGD.trees<-Unsure.WGD.trees[!sapply(Unsure.WGD.trees, is.null)] ##628 trees
rm(Validated.WGD.2species.unsure)

summary.unsure.WGD.trees<-summary.tree(Unsure.WGD.trees)
list.unsure.3RWGD<-as.vector(unique(summary.unsure.WGD.trees$index.tree))

#save.image("norm_exp_15spe_4tips_new.RData")


################### Chunk- 16: Identifying sure and unsure SSD trees #############

## We could have only excluded the possible WGD trees to identify SSD trees
## Since ohnolog database provides FishWGD data only for zebrafish and medaka among the five teleosts we considered, 
## we should exclude trees for which "Clupeocephala" or "Osteoglossocephalai" duplication is present 
## This will help to identify sure SSD trees for this study
## For trees where "Clupeocephala" or "Osteoglossocephalai" duplication is present, they are considered as "unsure" SSD trees 

##Step1: Extracting tree indices which are not possible candidates of 3RWGDs
all.index<-seq(1:6923)
index.to.exclude<-c(list.sure.3RWGD,list.unsure.3RWGD)
possible.ssd.index<-all.index[!(all.index %in% index.to.exclude)]
Other.group.trees<-trees.all.interest[possible.ssd.index] #5136 trees

##Modify duplication node labels
Other.group.trees.new<-mclapply(Other.group.trees,modify.label, mc.cores = cores)

## Identifying sure SSD trees 
exclude.trees<-NULL
exclude.trees<-mclapply(Other.group.trees.new,candidate.unsureSSD.phylogeny, mc.cores = cores)
sure.SSD.trees<-exclude.trees[!is.na(exclude.trees)] ##4139 gene trees passed
sure.SSD.trees<-mclapply(sure.SSD.trees,remodify.label, mc.cores = cores) ##remodified duplication labels
rm(Other.group.trees)
rm(Other.group.trees.new)

## Identifying index of sure SSD trees
summary.sure.SSD.trees<-summary.tree(sure.SSD.trees)
list.sure.SSD<-as.vector(unique(summary.sure.SSD.trees$index.tree))

## Identifying unsure SSD trees
unsure.SSD.index<-possible.ssd.index[!(possible.ssd.index %in% list.sure.SSD)] ##997 gene trees passed 
unsure.SSD.trees<-trees.all.interest[unsure.SSD.index] #997 trees
rm(exclude.trees)

#save.image("norm_exp_15spe_4tips_new.RData")

################ Chunk 17: Using the WGD trees to generate orthogroup to normalize the expression data and checking PCA plot#############################

freq<-0 
orthogroup.set<-na.omit(bind_rows(lapply(Sure.WGD.trees, orthogroup.finder)))
rm(freq)

################ Chunk 18: Clustering of normalized extression data for WGDs in several ways #############################
pca.WGD<-expression.all.sample(orthogroup.set,"orthogroup")

############### cluster type 1 - Random classification of ohnologs (considering any of the duplicates as duplicate1 and the other as duplicate2) ######################
## Generating PCA plot to assess the result of the analysis
PCA.plot(pca.WGD,type="common.WGD",group.by="tissue",method="normal.2D")

############## cluster type 2 - Highly expressed duplicates are called as "duplicate1" and lowly expressed duplicates as "duplicate2" ######################

## Identifying highly and lowly expressed duplicate based on normalized average expression value
orthogroup.set1<-duplicate.classfication1(orthogroup.set,pca.WGD)
colnames(orthogroup.set1)<-colnames(orthogroup.set)

## Preparing dataset before pca 
pca.WGD1<-expression.all.sample(orthogroup.set1,"orthogroup")

## Separating highly expressed and lowly expressed duplicates
OG.ID<-colnames(pca.WGD1)
HE.OG.ID<-OG.ID[grepl(".1",OG.ID)]
LE.OG.ID<-OG.ID[grepl(".2",OG.ID)]

pca.WGD1.high.exp<-pca.WGD1[ , !(colnames(pca.WGD1) %in% LE.OG.ID)]
pca.WGD1.low.exp<-pca.WGD1[ , !(colnames(pca.WGD1) %in% HE.OG.ID)]

PCA.plot(pca.WGD1,type="common.WGD",group.by="tissue",method="normal.2D")

####################### Chunk 19: Analysis on contrasts standardized 3R sure WGDs and sure SSD trees for traits ###############################

## Cluepeocephala.d 3rd round whole genome duplication time point
## '@both  studies  place  it  before  the teleost radiation, which is consistent with a 
## genome duplication at the base of ray-finned fishes' - ref 1 cited 2 references (2017 pdf)
## ref 1: Whole-Genome Duplication in Teleost Fishes and Its Evolutionary Consequences
## ref 2: From 2R to 3R: evidence for a fish-specific genome duplication (FSGD)
## Time taken from Ensembl database + duplication time from published paper
## ref:" doi:10.1038/ng.3645" Nature Genetics 2016

## Diagnostic tests for Tau for sure 3R WGD trees 
count<-0
Calibrated.3R.tau1<-lapply(Sure.WGD.trees, diagnostic.plot.test, 0)
Calibrated.3R.tau2<-Calibrated.3R.tau1[!is.na(Calibrated.3R.tau1)] 
Calibrated.3R.tau<-Calibrated.3R.tau2[ ! sapply(Calibrated.3R.tau2, is.null) ]## finally we obtained 746/1159 tree data passing the diagnostic tests
rm(Calibrated.3R.tau1)
rm(Calibrated.3R.tau2)

##Diagnostic tests for Tau for sure SSD trees 
count<-0
Calibrated.SSD.tau1<-lapply(sure.SSD.trees, diagnostic.plot.test, 0)
Calibrated.SSD.tau2<-Calibrated.SSD.tau1[!is.na(Calibrated.SSD.tau1)] 
Calibrated.SSD.tau<-Calibrated.SSD.tau2[ ! sapply(Calibrated.SSD.tau2, is.null) ]## finally we obtained 2482/4139 trees passing the diagnostic tests
rm(Calibrated.SSD.tau1)
rm(Calibrated.SSD.tau2)

## Diagnostic tests for expression levels of sure 3R WGDs
count<-0
Calibrated.3R.exp1<-lapply(Sure.WGD.trees, diagnostic.plot.test.exp, 0)
Calibrated.3R.exp2<-Calibrated.3R.exp1[!is.na(Calibrated.3R.exp1)] 
Calibrated.3R.exp<-Calibrated.3R.exp2[ ! sapply(Calibrated.3R.exp2, is.null) ]## Finally we obtained 774 trees 
rm(Calibrated.3R.exp1)
rm(Calibrated.3R.exp2)

## Diagnostic tests for expression levels of sure SSD trees
count<-0
Calibrated.SSD.exp1<-lapply(sure.SSD.trees, diagnostic.plot.test.exp, 0)
Calibrated.SSD.exp2<-Calibrated.SSD.exp1[!is.na(Calibrated.SSD.exp1)] 
Calibrated.SSD.exp<-Calibrated.SSD.exp2[!sapply(Calibrated.SSD.exp2, is.null) ]## Finally we obtained 2717 out of 4139 trees 
rm(Calibrated.SSD.exp1)
rm(Calibrated.SSD.exp2)

## Summarizing contrasts for Tau for sure 3R WGDs
summary.tau.3RWGDs<-summary.tree(Calibrated.3R.tau)
summary.tau.3RWGDs$Tau.abs<-abs(summary.tau.3RWGDs$pic_Tau)
summary.tau.3RWGDs<-summary.tau.3RWGDs[,!(names(summary.tau.3RWGDs) %in% c("pic_Tau"))]
nodes.contrast.tau.3RWGDs <- summary.tau.3RWGDs[which(!is.na(summary.tau.3RWGDs$Tau.abs)),]
nodes.contrast.tau.3RWGDs$Event <- factor(nodes.contrast.tau.3RWGDs$events, levels=c("speciation", "duplication","FishWGD"))
rm(summary.tau.3RWGDs)

## Summarizing contrasts for Tau for sure SSDs
summary.tau.SSDs<-summary.tree(Calibrated.SSD.tau)
summary.tau.SSDs$Tau.abs<-abs(summary.tau.SSDs$pic_Tau)
summary.tau.SSDs<-summary.tau.SSDs[,!(names(summary.tau.SSDs) %in% c("pic_Tau"))]
nodes.contrast.tau.SSDs <- summary.tau.SSDs[which(!is.na(summary.tau.SSDs$Tau.abs)),]
nodes.contrast.tau.SSDs$Event <- factor(nodes.contrast.tau.SSDs$events, levels=c("speciation", "duplication"))
rm(summary.tau.SSDs)

## Summarizing contrasts for normalized average gene expression levels for sure 3R WGDs
summary.exp.3RWGDs<-summary.tree(Calibrated.3R.exp)
summary.exp.3RWGDs$Exp.abs <- abs(summary.exp.3RWGDs$pic_Exp)
summary.exp.3RWGDs<-summary.exp.3RWGDs[,!(names(summary.exp.3RWGDs) %in% c("pic_Exp"))]
nodes.contrast.exp.3RWGDs <- summary.exp.3RWGDs[which(!is.na(summary.exp.3RWGDs$Exp.abs)),]
nodes.contrast.exp.3RWGDs$Event <- factor(nodes.contrast.exp.3RWGDs$events, levels=c("speciation", "duplication","FishWGD"))
rm(summary.exp.3RWGDs)

## Summarizing contrasts for normalized average gene expression levels for sure SSDs
summary.exp.SSDs<-summary.tree(Calibrated.SSD.exp)
summary.exp.SSDs$Exp.abs <- abs(summary.exp.SSDs$pic_Exp)
summary.exp.SSDss<-summary.exp.SSDs[,!(names(summary.exp.SSDs) %in% c("pic_Exp"))]
nodes.contrast.exp.SSDs <- summary.exp.SSDs[which(!is.na(summary.exp.SSDs$Exp.abs)),]
nodes.contrast.exp.SSDs$Event <- factor(nodes.contrast.exp.SSDs$events, levels=c("speciation", "duplication"))
rm(summary.exp.SSDs)

### For species trees summary of Tau
summary.tau.species<-summary.tree(Calibrated.speciation.trees.modify.standardized.tau)
nodes.contrast.tau.species <- summary.tau.species[which(!is.na(summary.tau.species$pic_Tau)),]
nodes.contrast.tau.species$Tau.abs <-abs(nodes.contrast.tau.species$pic_Tau)

### For species trees summary of Exp
summary.exp.species<-summary.tree(Calibrated.speciation.trees.modify.standardized.exp)
nodes.contrast.exp.species <- summary.exp.species[which(!is.na(summary.exp.species$pic_Exp)),]
nodes.contrast.exp.species$Exp.abs <-abs(nodes.contrast.exp.species$pic_Exp)
#save.image("norm_exp_15spe_4tips_new.RData")

######################### Chunk20: Contrast analyses by our new method considering "speciation" as root event #############################
## Since our method needs a tree to start with a speciation event, many trees fail to pass such criteria due to vertebrata or chordata duplication at the base

###### Part1: Considering oldest duplication nodes#####

#### For Tau on contrasts standardized trees
count<-0
All.duplicates.interest<-lapply(Calibrated.all.standardized.tau,special.speciation.criteria.new)
Duplicates.our.interest1<-All.duplicates.interest[!is.na(All.duplicates.interest)] 
Duplicates.our.interest<-Duplicates.our.interest1[!sapply(Duplicates.our.interest1, is.null)]## 1375 trees passed the criteria 
rm(Duplicates.our.interest1)

## Checking effects of duplications for Tau on all contrasts standardized trees
count<-0
duplicate.effect1<-NULL
duplicate.effect<-NULL
duplicate.effect1<-bind_rows(lapply(Duplicates.our.interest,pic.spe.dup.pair.old,"pic_Tau"))
duplicate.effect<-duplicate.effect1[complete.cases(duplicate.effect1),]
rm(duplicate.effect1)
duplicate.effect$lab<-"all_tau"

#### For exp levels on contrasts standardized trees
count<-0
All.duplicates.interest.exp<-lapply(Calibrated.all.standardized.exp,special.speciation.criteria.new)
Duplicates.our.interest1.exp<-All.duplicates.interest.exp[!is.na(All.duplicates.interest.exp)] 
Duplicates.our.interest.exp<-Duplicates.our.interest1.exp[!sapply(Duplicates.our.interest1.exp, is.null)]## 1493 trees passed the criteria
rm(Duplicates.our.interest1.exp)

## Checking effects of duplications for exp levels on all contrasts standardized trees
count<-0
duplicate.effect1.exp<-NULL
duplicate.effect.exp<-NULL
duplicate.effect1.exp<-bind_rows(lapply(Duplicates.our.interest.exp,pic.spe.dup.pair.old,"pic_Exp"))
duplicate.effect.exp<-duplicate.effect1.exp[complete.cases(duplicate.effect1.exp),]
rm(duplicate.effect1.exp)
duplicate.effect.exp$lab<-"all_exp"

#### For Tau on contrasts standardized trees sure 3RWGD trees
count<-0
All.duplicates.interest.3R<-lapply(Calibrated.3R.tau,special.speciation.criteria.3RWGD.old)
Duplicates.our.interest1.3R<-All.duplicates.interest.3R[!is.na(All.duplicates.interest.3R)] 
Duplicates.our.interest.3R<-Duplicates.our.interest1.3R[!sapply(Duplicates.our.interest1.3R, is.null)]##581/746 trees passed the criteria
rm(Duplicates.our.interest1.3R)

## Checking effects of duplications for Tau on contrasts standardized trees 3RWGD trees
count<-0
duplicate.effect13R<-NULL
duplicate.effect3R<-NULL
duplicate.effect13R<-bind_rows(lapply(Duplicates.our.interest.3R,pic.speciation.pair.3RWGD,"pic_Tau"))
duplicate.effect3R1<-duplicate.effect13R[complete.cases(duplicate.effect13R),]
duplicate.effect3R<-duplicate.effect3R1[which(duplicate.effect3R1$duplication.label=="Clupeocephala" | duplicate.effect3R1$duplication.label=="Osteoglossocephalai"),]
duplicate.effect3R$lab<-"sure.3RWGD.Tau"
rm(duplicate.effect13R)
rm(duplicate.effect3R1)

#### For exp levels on contrasts standardized trees sure 3RWGD trees
count<-0
All.duplicates.interest.3R.exp<-lapply(Calibrated.3R.exp,special.speciation.criteria.3RWGD.old)
Duplicates.our.interest1.3R.exp<-All.duplicates.interest.3R.exp[!is.na(All.duplicates.interest.3R.exp)] 
Duplicates.our.interest.3R.exp<-Duplicates.our.interest1.3R.exp[!sapply(Duplicates.our.interest1.3R.exp, is.null)]##610/774 trees passed the criteria
rm(Duplicates.our.interest1.3R.exp)

## Checking effects of duplications for exp levels on contrasts standardized trees 3RWGD trees
count<-0
duplicate.effect13R.exp<-NULL
duplicate.effect3R.exp<-NULL
duplicate.effect13R.exp<-bind_rows(lapply(Duplicates.our.interest.3R.exp,pic.speciation.pair.3RWGD,"pic_Exp"))
duplicate.effect3R1.exp<-duplicate.effect13R.exp[complete.cases(duplicate.effect13R.exp),]
duplicate.effect3R.exp<-duplicate.effect3R1.exp[which(duplicate.effect3R1.exp$duplication.label=="Clupeocephala" | duplicate.effect3R1.exp$duplication.label=="Osteoglossocephalai"),]
duplicate.effect3R.exp$lab<-"sure.3RWGD.exp"
rm(duplicate.effect13R.exp)
rm(duplicate.effect3R1.exp)

#### For Tau on contrasts standardized trees sure SSD trees
count<-0
All.duplicates.interest.SSD<-lapply(Calibrated.SSD.tau,special.speciation.criteria.new)
All.duplicates.interest1.SSD<-All.duplicates.interest.SSD[!is.na(All.duplicates.interest.SSD)] 
All.duplicates.interest.SSD<-All.duplicates.interest1.SSD[!sapply(All.duplicates.interest1.SSD, is.null)]##only 339 trees passed the criteria
rm(All.duplicates.interest1.SSD)

## Checking effects of duplications for Tau on contrasts standardized trees SSD trees
count<-0
duplicate.effect.SSDn<-NULL
duplicate.effect.SSD<-NULL
duplicate.effect.SSDn<-bind_rows(lapply(All.duplicates.interest.SSD,pic.spe.dup.pair.old,"pic_Tau"))
duplicate.effect.SSD<-duplicate.effect.SSDn[complete.cases(duplicate.effect.SSDn),]
duplicate.effect.SSD$lab<-"sure.SSD.Tau"
rm(duplicate.effect.SSDn)

#### For exp levels on contrasts standardized trees sure SSD trees
count<-0
All.duplicates.interest.SSD.exp<-lapply(Calibrated.SSD.exp,special.speciation.criteria.new)
Duplicates.our.interest1.SSD.exp<-All.duplicates.interest.SSD.exp[!is.na(All.duplicates.interest.SSD.exp)] 
Duplicates.our.interest.SSD.exp<-Duplicates.our.interest1.SSD.exp[!sapply(Duplicates.our.interest1.SSD.exp, is.null)]##only 377 trees passed the criteria
rm(Duplicates.our.interest1.SSD.exp)

## Checking effects of duplications for exp levels on contrasts standardized trees SSD trees
count<-0
duplicate.effect.SSDn.exp<-NULL
duplicate.effect.SSD.exp<-NULL
duplicate.effect.SSDn.exp<-bind_rows(lapply(Duplicates.our.interest.SSD.exp,pic.spe.dup.pair.old,"pic_Exp"))
duplicate.effect.SSD.exp<-duplicate.effect.SSDn.exp[complete.cases(duplicate.effect.SSDn.exp),]
duplicate.effect.SSD.exp$lab<-"sure.SSD.Exp"
rm(duplicate.effect.SSDn.exp)

#save.image("norm_exp_15spe_4tips_new.RData")


####### Part2: Considering youngest duplication nodes ########

## ## Checking effects of duplications for Tau on all contrasts standardized trees
count<-0
duplicate.effect1.young<-NULL
duplicate.effect.young<-NULL
duplicate.effect1.young<-bind_rows(lapply(Duplicates.our.interest,pic.spe.dup.pair.young,"pic_Tau"))
duplicate.effect.young<-duplicate.effect1.young[complete.cases(duplicate.effect1.young),]
rm(duplicate.effect1.young)
duplicate.effect.young$lab<-"all_tau"

## Checking effects of duplications for exp levels on all contrasts standardized trees
count<-0
duplicate.effect1.young.exp<-NULL
duplicate.effect.young.exp<-NULL
duplicate.effect1.young.exp<-bind_rows(lapply(Duplicates.our.interest.exp,pic.spe.dup.pair.young,"pic_Exp"))
duplicate.effect.young.exp<-duplicate.effect1.young.exp[complete.cases(duplicate.effect1.young.exp),]
rm(duplicate.effect1.young.exp)
duplicate.effect.young.exp$lab<-"all_exp"

#### For Tau on contrasts standardized trees sure 3RWGD trees
count<-0
All.duplicates.interest.3R.young<-lapply(Calibrated.3R.tau,special.speciation.criteria.3RWGD.young)
Duplicates.our.interest1.3R.young<-All.duplicates.interest.3R.young[!is.na(All.duplicates.interest.3R.young)] 
Duplicates.our.interest.3R.young<-Duplicates.our.interest1.3R.young[!sapply(Duplicates.our.interest1.3R.young, is.null)]##581/746 trees passed the criteria
rm(Duplicates.our.interest1.3R.young)

## Checking effects of duplications for Tau on contrasts standardized trees 3RWGD trees
count<-0
duplicate.effect13R.young<-NULL
duplicate.effect3R.young<-NULL
duplicate.effect13R.young<-bind_rows(lapply(Duplicates.our.interest.3R.young,pic.speciation.pair.3RWGD.young,"pic_Tau"))
duplicate.effect3R1.young<-duplicate.effect13R.young[complete.cases(duplicate.effect13R.young),]
duplicate.effect3R.young<-duplicate.effect3R1.young[which(duplicate.effect3R1.young$duplication.label=="Clupeocephala" | duplicate.effect3R1.young$duplication.label=="Osteoglossocephalai"),]
duplicate.effect3R.young$lab<-"sure.3RWGD.Tau"
rm(duplicate.effect13R.young)

## Checking effects of duplications for exp levels on contrasts standardized trees 3RWGD trees
count<-0
duplicate.effect13R.young.exp<-NULL
duplicate.effect3R.young.exp<-NULL
duplicate.effect13R.young.exp<-bind_rows(lapply(Duplicates.our.interest.3R.exp,pic.speciation.pair.3RWGD.young,"pic_Exp"))
duplicate.effect3R1.young.exp<-duplicate.effect13R.young.exp[complete.cases(duplicate.effect13R.young.exp),]
duplicate.effect3R.young.exp<-duplicate.effect3R1.young.exp[which(duplicate.effect3R1.young.exp$duplication.label=="Clupeocephala"| duplicate.effect3R1.young.exp$duplication.label=="Osteoglossocephalai"),]
duplicate.effect3R.young.exp$lab<-"sure.3RWGD.exp"
rm(duplicate.effect13R.young.exp)
rm(duplicate.effect3R1.young.exp)

## Checking effects of duplications for Tau on contrasts standardized trees SSD trees
count<-0
duplicate.effect.SSDn.young<-NULL
duplicate.effect.SSD.young<-NULL
duplicate.effect.SSDn.young<-bind_rows(lapply(All.duplicates.interest.SSD,pic.spe.dup.pair.young,"pic_Tau"))
duplicate.effect.SSD.young<-duplicate.effect.SSDn.young[complete.cases(duplicate.effect.SSDn.young),]
duplicate.effect.SSD.young$lab<-"sure.SSD.Tau"
rm(duplicate.effect.SSDn.young)

## Checking effects of duplications for exp levels on contrasts standardized trees SSD trees
count<-0
duplicate.effect.SSDn.young.exp<-NULL
duplicate.effect.SSD.young.exp<-NULL
duplicate.effect.SSDn.young.exp<-bind_rows(lapply(Duplicates.our.interest.SSD.exp,pic.spe.dup.pair.young,"pic_Exp"))
duplicate.effect.SSD.young.exp<-duplicate.effect.SSDn.young.exp[complete.cases(duplicate.effect.SSDn.young.exp),]
duplicate.effect.SSD.young.exp$lab<-"sure.SSD.Exp"
rm(duplicate.effect.SSDn.young.exp)

#save.image("norm_exp_15spe_4tips_new.RData")

########################## Chunk 21: Trait randomization test for sure SSD and for sure 3RWGD using our method ################################

## In this test, we used tip shuffling for our trees, so that we use the same nodes for comparison
## After randomization, the trees no longer follow the Brownian model of trait evolution

## For calibrated sure 3RWGD trees 
## Permuted actual Tau of the tips 
Tau.randomized.tree.3R <- mclapply(Calibrated.3R.tau,shuffling.tau,mc.cores = cores)

Count<-0
Tree.pic.randomized.3R <- lapply(Tau.randomized.tree.3R,contrast.calc,trait="Tau")
rm(Count)
rm(Tau.randomized.tree.3R)

## Now we need to check for trees with contrast standardization 
count<-0
Randomized.tau.tree.3R1<-mclapply(Tree.pic.randomized.3R,special.speciation.criteria.3RWGD.old,mc.cores = cores)
Randomized.tau.tree.3R2<-Randomized.tau.tree.3R1[!is.na(Randomized.tau.tree.3R1)] 
Randomized.tau.tree.3R<-Randomized.tau.tree.3R2[!sapply(Randomized.tau.tree.3R2, is.null)] ##only 581 trees passed the criteria
rm(Randomized.tau.tree.3R1)
rm(Randomized.tau.tree.3R2)
rm(count)

## Effect of randomization on Tau
count<-0
duplicate.effect1.3R.random<-NULL
duplicate.effect.3R.random<-NULL
duplicate.effect1.3R.random<-bind_rows(lapply(Randomized.tau.tree.3R,pic.speciation.pair.3RWGD,"pic_Tau"))
duplicate.effect.3R.random<-duplicate.effect1.3R.random[complete.cases(duplicate.effect1.3R.random),]
duplicate.effect.3R.random<-duplicate.effect.3R.random[which(duplicate.effect.3R.random$duplication.label=="Clupeocephala" | duplicate.effect.3R.random$duplication.label=="Osteoglossocephalai"),]
rm(duplicate.effect1.3R.random)

## For calibrated sure SSD trees 
## Permuted actual Tau of the tips 
Tau.randomized.tree.SSD <- mclapply(Calibrated.SSD.tau,shuffling.tau,mc.cores = cores)

Count<-0
Tree.pic.randomized.SSD <- lapply(Tau.randomized.tree.SSD,contrast.calc,trait="Tau")
rm(Count)
rm(Tau.randomized.tree.SSD)

## Now we need to check for trees with contrast standardization 
count<-0
Randomized.tau.tree.SSD1<-mclapply(Tree.pic.randomized.SSD,special.speciation.criteria.new,mc.cores = cores)
Randomized.tau.tree.SSD2<-Randomized.tau.tree.SSD1[!is.na(Randomized.tau.tree.SSD1)] 
Randomized.tau.tree.SSD<-Randomized.tau.tree.SSD2[!sapply(Randomized.tau.tree.SSD2, is.null)] ##only 339 trees passed the criteria
rm(Randomized.tau.tree.SSD1)
rm(Randomized.tau.tree.SSD2)
rm(count)

## Effect of randomization on Tau
count<-0
duplicate.effect1.SSD.random<-NULL
duplicate.effect.SSD.random<-NULL
duplicate.effect1.SSD.random<-bind_rows(lapply(Randomized.tau.tree.SSD,pic.spe.dup.pair.old,"pic_Tau"))
duplicate.effect.SSD.random<-duplicate.effect1.SSD.random[complete.cases(duplicate.effect1.SSD.random),]
rm(duplicate.effect1.SSD.random)

## For calibrated sure WGD trees 
## Permuted actual Exp levels of the tips 
Exp.randomized.tree.3R <- mclapply(Calibrated.3R.exp,shuffling.exp,mc.cores = cores)

Count<-0
Tree.pic.randomized.3R.Exp <- lapply(Exp.randomized.tree.3R,contrast.calc,trait="Exp")
rm(Count)
rm(Exp.randomized.tree.3R)

## Now we need to check for trees with contrast standardization 
count<-0
Randomized.exp.tree.3R1<-mclapply(Tree.pic.randomized.3R.Exp,special.speciation.criteria.3RWGD.old,mc.cores = cores)
Randomized.exp.tree.3R2<-Randomized.exp.tree.3R1[!is.na(Randomized.exp.tree.3R1)] 
Randomized.exp.tree.3R<-Randomized.exp.tree.3R2[!sapply(Randomized.exp.tree.3R2, is.null)] ##only 610 trees passed the criteria
rm(Randomized.exp.tree.3R1)
rm(Randomized.exp.tree.3R2)
rm(count)

## Effect of randomization on Exp
count<-0
duplicate.effect1.3R.exp.random<-NULL
duplicate.effect.3R.exp.random<-NULL
duplicate.effect1.3R.exp.random<-bind_rows(lapply(Randomized.exp.tree.3R,pic.speciation.pair.3RWGD,"pic_Exp"))
duplicate.effect.3R.exp.random<-duplicate.effect1.3R.exp.random[complete.cases(duplicate.effect1.3R.exp.random),]
duplicate.effect.3R.exp.random<-duplicate.effect.3R.exp.random[which(duplicate.effect.3R.exp.random$duplication.label=="Clupeocephala" | duplicate.effect.3R.exp.random$duplication.label=="Osteoglossocephalai"),]
rm(duplicate.effect1.3R.exp.random)

## For calibrated sure SSD trees 
## Permuted actual Exp levels of the tips 
Exp.randomized.tree.SSD <- mclapply(Calibrated.SSD.exp,shuffling.exp,mc.cores = cores)

Count<-0
Tree.pic.randomized.SSD.Exp <- lapply(Exp.randomized.tree.SSD,contrast.calc,trait="Exp")
rm(Count)
rm(Exp.randomized.tree.SSD)

## Now we need to check for trees with contrast standardization 
count<-0
Randomized.exp.tree.SSD1<-mclapply(Tree.pic.randomized.SSD.Exp,special.speciation.criteria.new,mc.cores = cores)
Randomized.exp.tree.SSD2<-Randomized.exp.tree.SSD1[!is.na(Randomized.exp.tree.SSD1)] 
Randomized.exp.tree.SSD<-Randomized.exp.tree.SSD2[!sapply(Randomized.exp.tree.SSD2, is.null)] ##only 377 trees passed the criteria
rm(Randomized.exp.tree.SSD1)
rm(Randomized.exp.tree.SSD2)
rm(count)

## Effect of randomization on Exp
count<-0
duplicate.effect1.SSD.exp.random<-NULL
duplicate.effect.SSD.exp.random<-NULL
duplicate.effect1.SSD.exp.random<-bind_rows(lapply(Randomized.exp.tree.SSD,pic.spe.dup.pair.old,"pic_Exp"))
duplicate.effect.SSD.exp.random<-duplicate.effect1.SSD.exp.random[complete.cases(duplicate.effect1.SSD.exp.random),]
rm(duplicate.effect1.SSD.exp.random)

#save.image("norm_exp_15spe_4tips_new.RData")

########################## Chunk 22: Trait randomization test for contrasts standardized trees using our method ################################

## In this test, we used tip shuffling for our trees, so that we use the same nodes for comparison
## After randomization, the trees no longer follow the Brownian model of trait evolution

## For calibrated trees of Tau  
## Permuted actual Tau of the tips 
Tau.randomized.tree <- mclapply(Calibrated.all.standardized.tau,shuffling.tau,mc.cores = cores)

Count<-0
Tree.pic.randomized <- lapply(Tau.randomized.tree,contrast.calc,trait="Tau")
rm(Count)
rm(Tau.randomized.tree)

## Now we need to check for trees with contrast standardization 
count<-0
Randomized.tau.tree1<-mclapply(Tree.pic.randomized,special.speciation.criteria.new,mc.cores = cores)
Randomized.tau.tree2<-Randomized.tau.tree1[!is.na(Randomized.tau.tree1)] 
Randomized.tau.tree<-Randomized.tau.tree2[!sapply(Randomized.tau.tree2, is.null)] ##only 1375 trees passed the criteria
rm(Randomized.tau.tree1)
rm(Randomized.tau.tree2)
rm(count)

## Effect of randomization on Tau
count<-0
duplicate.effect1.tau.random<-NULL
duplicate.effect.tau.random<-NULL
duplicate.effect1.tau.random<-bind_rows(lapply(Randomized.tau.tree,pic.spe.dup.pair.old,"pic_Tau"))
duplicate.effect.tau.random<-duplicate.effect1.tau.random[complete.cases(duplicate.effect1.tau.random),]
rm(duplicate.effect1.tau.random)

## For calibrated trees of exp
## Permuted actual Exp levels of the tips 
Exp.randomized.tree <- mclapply(Calibrated.all.standardized.exp,shuffling.exp,mc.cores = cores)

Count<-0
Tree.pic.randomized.Exp <- lapply(Exp.randomized.tree,contrast.calc,trait="Exp")
rm(Count)
rm(Exp.randomized.tree)

## Now we need to check for trees with contrast standardization 
count<-0
Randomized.exp.tree1<-mclapply(Tree.pic.randomized.Exp,special.speciation.criteria.new,mc.cores = cores)
Randomized.exp.tree2<-Randomized.exp.tree1[!is.na(Randomized.exp.tree1)] 
Randomized.exp.tree<-Randomized.exp.tree2[!sapply(Randomized.exp.tree2, is.null)] ##only 1493 trees passed the criteria
rm(Randomized.exp.tree1)
rm(Randomized.exp.tree2)
rm(count)

## Effect of randomization on Exp
count<-0
duplicate.effect1.exp.random<-NULL
duplicate.effect.exp.random<-NULL
duplicate.effect1.exp.random<-bind_rows(lapply(Randomized.exp.tree,pic.spe.dup.pair.old,"pic_Exp"))
duplicate.effect.exp.random<-duplicate.effect1.exp.random[complete.cases(duplicate.effect1.exp.random),]
rm(duplicate.effect1.exp.random)

#save.image("Preprocess.rda")
save.image("norm_exp_15spe_4tips_new.RData")



