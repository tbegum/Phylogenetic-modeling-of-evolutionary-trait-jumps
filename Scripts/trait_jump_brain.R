

### This code is written to work with "trait jump model" for brain

##########Loaded required libraries, setting path and calling required functions 
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
library(tidyverse)
library(digest)
library(ggplot2)
library(gtools)
library(caper)
library(cowplot)
library(easyGgplot2)
library(ggplot2)
library(gridExtra)
library(R.utils)
library(hash)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")
library(edgeR)

##Setting working directory to use functions
setwd("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/")
folder1 <- c("~/Desktop/nwork/Empirical/Manuscript_3/processed_data_files/")
source("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/functions_TPCM.R")

## Setting seed number for reproducibility 
set.seed(786)

## To perform parallel computation we need to detect core numbers
cores <- detectCores()

## loading previously saved data
load("norm_exp_15spe_4tips_new.RData")  ## latest modified exp 

######################### Chunk1: Data given to pablo for applying "Levolution" to find trait jump #############################
## trees.all.interest with 6923 gene trees given to pablo

######################### Chunk2: reading output of trait jump #############################

test.output.brain<-read.csv(paste0(folder1,"To_pablo/Brain/results_Brain.csv", sep=""),sep=",",header=T) ##168740 obs

## Only considering trees with posterior probability of jump more than 70%
test.output.brain1<-test.output.brain[which(test.output.brain$JumpProbability>=0.7),] ##33010 obs


########################## Chunk 3: Analysis on all contrasts standardized trees for traits ###############################

### Step 1a: Considering (trees.all.interest) 6923 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for brain tissue for all trees 
count<-0
Calibrated.brain.exp1<-lapply(trees.all.interest, diagnostic.plot.test.exp.all.tissue, "TPM.brain")
Calibrated.brain.exp2<-Calibrated.brain.exp1[!is.na(Calibrated.brain.exp1)] 
Calibrated.all.standardized.brain<-Calibrated.brain.exp2[ ! sapply(Calibrated.brain.exp2, is.null) ]## finally we obtained 4662/6923 tree data passing the diagnostic tests
rm(Calibrated.brain.exp1)
rm(Calibrated.brain.exp2)

### Step 1b: Considering (Calibrated.all.standardized.brain) 4662 tree data passing the diagnostic test for brain expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.all.standardized.brain,tree.index.collect))
rm(count)

test.output.brain.all<-NULL
test.output.brain.all<-merge(index.info,test.output.brain1,by="Gene") ## 19138 observations
index.trees.jump<-unique(test.output.brain.all$tree.num) ##2854 trees
Calibrated.standardized.brain.jump<-Calibrated.all.standardized.brain[c(index.trees.jump)] ##2854/4662=61.21% Brownian trees supported trait jump model for brain expressions 
Calibrated.standardized.brain.nojump<-Calibrated.all.standardized.brain[-c(index.trees.jump)] ## 1808/4662=38.78% trees with no support for jump in brain expression
rm(index.info)
rm(index.trees.jump)

### Step2: Calculating PICs of required trait for trees to match the index number of the trees on which levolution was run
Count<-0
trees.all.interest<-lapply(trees.all.interest,contrast.calc,trait="TPM.brain")
# Count<-0
# Calibrated.standardized.brain.jump<-lapply(Calibrated.standardized.brain.jump,contrast.calc,trait="TPM.brain") 
# rm(Count)

### Step3: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.brain.all))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.brain.all)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 19318 obs
}

test.output.brain.processed<-bind_cols(test.output.brain.all,processed.fout) ## 19318 obs for all the 15 species
rm(i)
rm(test.output.brain.all)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

## Step 4: Analysis on all 15 vertebrates clades
test.output.brain.new.mod<-unique(test.output.brain.processed[c(6,3,7,8,10,2,1)])  #14511 obs

########### Proportion of jump events in 15 vertebrates species in the 2854 trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.trees.brain<-NULL
count<-0
info.jump.trees.brain<-bind_rows(lapply(Calibrated.standardized.brain.jump, tree.data.collection))
info.jump.trees.brain$spe.num<-info.jump.trees.brain$internal.events-info.jump.trees.brain$dup.num
info.jump.trees.brain<-info.jump.trees.brain[c(11,12,3:5)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.brain<-
  test.output.brain.processed  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.brain)<-c("index.old","jump.dup")

test.spe.brain<-
  test.output.brain.processed%>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.brain)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.brain<-full_join(test.spe.brain,test.dup.brain)
info.jump.trees.brain.final<-merge(info.jump.trees.brain,merged.brain, by = c("index.old"))
info.jump.trees.brain.final[is.na(info.jump.trees.brain.final)] <- 0

## Have to consider trees where atleast one duplication and speciation exists
info.jump.trees.brain.final1<-info.jump.trees.brain.final[which(info.jump.trees.brain.final$dup.num>0 & info.jump.trees.brain.final$spe.num>0),]
rm(info.jump.trees.brain.final)

## Proportion of jumps
info.jump.trees.brain.final1$dup.jump.prob<-(info.jump.trees.brain.final1$jump.dup/(info.jump.trees.brain.final1$dup.num*2))
info.jump.trees.brain.final1$spe.jump.prob<-(info.jump.trees.brain.final1$jump.spe/(info.jump.trees.brain.final1$spe.num*2))

## Stats
median(na.omit(info.jump.trees.brain.final1$dup.jump.prob)) #0.25
median(na.omit(info.jump.trees.brain.final1$spe.jump.prob)) #0.0769
mean(na.omit(info.jump.trees.brain.final1$dup.jump.prob)) #0.313
mean(na.omit(info.jump.trees.brain.final1$spe.jump.prob)) #0.1604

pval.jump.brain <- paired.wilcox(info.jump.trees.brain.final1$spe.jump.prob,info.jump.trees.brain.final1$dup.jump.prob) #p-value = 2.2e-16

## For median data
dtfr.brain<-NULL
dtfr.brain<-data.frame(prop.jump=info.jump.trees.brain.final1$spe.jump.prob,events="speciation")
dtfr.brain<-rbind(dtfr.brain,data.frame(prop.jump=info.jump.trees.brain.final1$dup.jump.prob,events="duplication"))
dtfr.brain$Exp.abs <-abs(dtfr.brain$prop.jump)
dtfr.brain$Event<-factor(dtfr.brain$events, levels=c("speciation","duplication"))

library(tidyverse)
median.data.jump.brain<-NULL
median.data.jump.brain<- dtfr.brain%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.brain$pic.round<-paste0("Median = ",round(median.data.jump.brain$Exp.abs,4), sep="")
median.data.jump.brain$count<- paste0(median.data.jump.brain$Event,"\n\n","(n = ",median.data.jump.brain$freq,")")

plot4A<-boxplot.new2(dtfr.brain,pval.jump.brain,"Exp",median.data.jump.brain) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.brain, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.brain$count) ## check

# Create Data
data.jump.brain <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(2854,1808)
)

# Compute the position of labels
data.jump.brain <- data.jump.brain %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.brain$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.brain$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.brain<-ggplot(data.jump.brain, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 2854 trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.brain.jump<-summary.tree(Calibrated.standardized.brain.jump)
summary.brain.jump$Exp.abs <- abs(summary.brain.jump$pic_TPM.brain)
summary.brain.jumps<-summary.brain.jump[,!(names(summary.brain.jump) %in% c("pic_TPM.brain"))]
nodes.contrast.brain.jumps.15spe <- summary.brain.jumps[which(!is.na(summary.brain.jumps$Exp.abs)),]
nodes.contrast.brain.jumps.15spe$Event <- factor(nodes.contrast.brain.jumps.15spe$events, levels=c("speciation", "duplication"))
rm(summary.brain.jump)
rm(summary.brain.jumps)

############ Testing OC using all 15 species for 2854 trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.brain.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.brain.jumps.15spe$index<-nodes.contrast.brain.jumps.15spe$index.tree
test.output.brain.processed$node<-test.output.brain.processed$From
nodes.contrast.brain.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.brain.jumps.15spe,test.output.brain.processed,by=c("index","node"))
nodes.contrast.brain.jumps.15spe.removing.jump.node <- nodes.contrast.brain.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.brain.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.brain.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.brain.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.brain.jumps.15spe.removing.jump.node$Exp.abs<-nodes.contrast.brain.jumps.15spe.removing.jump.node$pic


########### Testing OC for all nodes of 1808 trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.brain.nojump<-summary.tree(Calibrated.standardized.brain.nojump)
summary.brain.nojump$Exp.abs <- abs(summary.brain.nojump$pic_TPM.brain)
summary.brain.nojumps<-summary.brain.nojump[,!(names(summary.brain.nojump) %in% c("pic_TPM.brain"))]
nodes.contrast.brain.nojumps.15spe <- summary.brain.nojump[which(!is.na(summary.brain.nojump$Exp.abs)),]
nodes.contrast.brain.nojumps.15spe$Event <- factor(nodes.contrast.brain.nojumps.15spe$events, levels=c("speciation", "duplication"))
rm(summary.brain.nojump)
rm(summary.brain.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.brain.new.fish<-test.output.brain.processed[which(test.output.brain.processed$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##7591 obs

## Is OC supported for exp for fish clades supporting jump?
pic.spe.fish<-test.output.brain.new.fish$pic[which(test.output.brain.new.fish$node.event=="speciation")]
pic.dup.fish<-test.output.brain.new.fish$pic[which(test.output.brain.new.fish$node.event=="duplication")]
median(pic.spe.fish) ##0.057
median(pic.dup.fish) ##0.1118
Pval.brain.jump.fish<-two.tailed.wilcox(pic.spe.fish,pic.dup.fish) # < 2.2 e-16

## Step 6: Analysis on fish specific clades
test.output.brain.new.fish.mod<-unique(test.output.brain.new.fish[c(6,3,7,8,10,2,1)])  #5595 obs

########### Proportion of jump events in teleosts in the 2854 trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.trees.brain.fish<-NULL
count<-0
info.jump.trees.brain.fish<-bind_rows(lapply(Calibrated.standardized.brain.jump, tree.data.collection.fish))
info.jump.trees.brain.fish$spe.num<-info.jump.trees.brain.fish$internal.events.fish-info.jump.trees.brain.fish$dup.num
info.jump.trees.brain.fish<-info.jump.trees.brain.fish[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.brain.fish<-
  test.output.brain.new.fish  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.brain.fish)<-c("index.old","jump.dup")

test.spe.brain.fish<-
  test.output.brain.new.fish %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.brain.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.brain.fish<-full_join(test.spe.brain.fish,test.dup.brain.fish)
info.jump.trees.brain.final.fish<-merge(info.jump.trees.brain.fish,merged.brain.fish, by = c("index.old"))
info.jump.trees.brain.final.fish[is.na(info.jump.trees.brain.final.fish)] <- 0

## Have to consider trees where atleast one fish specific duplication and speciation exists
info.jump.trees.brain.final1.fish<-info.jump.trees.brain.final.fish[which(info.jump.trees.brain.final.fish$dup.num>0 & info.jump.trees.brain.final.fish$spe.num>0),]

## Proportion of jumps
info.jump.trees.brain.final1.fish$dup.jump.prob<-(info.jump.trees.brain.final1.fish$jump.dup/(info.jump.trees.brain.final1.fish$dup.num*2))
info.jump.trees.brain.final1.fish$spe.jump.prob<-(info.jump.trees.brain.final1.fish$jump.spe/(info.jump.trees.brain.final1.fish$spe.num*2))


median(na.omit(info.jump.trees.brain.final1.fish$dup.jump.prob)) #0.25
median(na.omit(info.jump.trees.brain.final1.fish$spe.jump.prob)) #0.167
mean(na.omit(info.jump.trees.brain.final1.fish$dup.jump.prob)) #0.3077
mean(na.omit(info.jump.trees.brain.final1.fish$spe.jump.prob)) #0.2628

pval.jump.brain.fish <- paired.wilcox(info.jump.trees.brain.final1.fish$spe.jump.prob,info.jump.trees.brain.final1.fish$dup.jump.prob) #p-value = 3.87e-08


## For median data
dtfr.brain.fish<-NULL
dtfr.brain.fish<-data.frame(prop.jump=info.jump.trees.brain.final1.fish$spe.jump.prob,events="speciation")
dtfr.brain.fish<-rbind(dtfr.brain.fish,data.frame(prop.jump=info.jump.trees.brain.final1.fish$dup.jump.prob,events="duplication"))
dtfr.brain.fish$Exp.abs <-abs(dtfr.brain.fish$prop.jump)
dtfr.brain.fish$Event<-factor(dtfr.brain.fish$events, levels=c("speciation","duplication"))


library(tidyverse)
median.data.jump.brain.brain.fish<-NULL
median.data.jump.brain.brain.fish<- dtfr.brain.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.brain.brain.fish$pic.round<-paste0("Median = ",round(median.data.jump.brain.brain.fish$Exp.abs,4), sep="")
median.data.jump.brain.brain.fish$count<- paste0(median.data.jump.brain.brain.fish$Event,"\n\n","(n = ",median.data.jump.brain.brain.fish$freq,")")

plot4B<-boxplot.new2(dtfr.brain.fish,pval.jump.brain.fish,"Exp",median.data.jump.brain.brain.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 5 teleosts"))))+
  #ylab(expression(bold(paste("Proportions of jump events for teleosts ")))) +
  ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.brain.fish, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.brain.brain.fish$count) ## check

########### Testing OC for all fish specific nodes of 2854 trees supporting jump##################
nodes.contrast.brain.jumps.fish<-nodes.contrast.brain.jumps.15spe[which(nodes.contrast.brain.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 2854 trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.brain.jumps.fish.removing.jump.node <- NULL
nodes.contrast.brain.jumps.fish$index<-nodes.contrast.brain.jumps.fish$index.tree
test.output.brain.new.fish$node<-test.output.brain.new.fish$From
nodes.contrast.brain.jumps.fish.removing.jump.node <-anti_join(nodes.contrast.brain.jumps.fish,test.output.brain.new.fish,by=c("index","node"))
nodes.contrast.brain.jumps.fish.removing.jump.node <- nodes.contrast.brain.jumps.fish.removing.jump.node[which(!is.na(nodes.contrast.brain.jumps.fish.removing.jump.node$pic)),]
nodes.contrast.brain.jumps.fish.removing.jump.node$Event <- factor(nodes.contrast.brain.jumps.fish.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.brain.jumps.fish.removing.jump.node$Exp.abs<-nodes.contrast.brain.jumps.fish.removing.jump.node$pic

pic.spe.removing.jump.fish.2854trees<-speciation.contrast(nodes.contrast.brain.jumps.fish.removing.jump.node,"pic") #(n=10105)
pic.dup.removing.jump.fish.2854trees<-duplication.contrast(nodes.contrast.brain.jumps.fish.removing.jump.node,"pic") #(n=2398)

rm(test.output.brain.new.fish)

########### Testing OC for all fish specific nodes of 1808 trees that do not support trait jump for teleosts ##################
nodes.contrast.brain.nojumps.fish<-nodes.contrast.brain.nojumps.15spe[which(nodes.contrast.brain.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")

################## Chunk 4: Analysis on contrasts standardized sure FishWGD trees for traits ##########################

### Step 1a: Considering (Sure.WGD.trees) 1159 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for brain tissue for sure WGD trees 
count<-0
Calibrated.brain.3Rexp1<-lapply(Sure.WGD.trees, diagnostic.plot.test.exp.all.tissue, "TPM.brain")
Calibrated.brain.3Rexp2<-Calibrated.brain.3Rexp1[!is.na(Calibrated.brain.3Rexp1)] 
Calibrated.3R.standardized.brain<-Calibrated.brain.3Rexp2[ ! sapply(Calibrated.brain.3Rexp2, is.null) ]## finally we obtained 800/1159 trees passing the diagnostic tests
rm(Calibrated.brain.3Rexp1)
rm(Calibrated.brain.3Rexp2)

### Step 1b: Considering (Calibrated.3R.standardized.brain) 800 trees passing the diagnostic test for brain expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.3R.standardized.brain,tree.index.collect))
rm(count)

test.output.3Rexp.brain<-NULL
test.output.3Rexp.brain<-merge(index.info,test.output.brain1,by="Gene") ## 2716 observations
index.trees.jump<-unique(test.output.3Rexp.brain$tree.num) ##445 trees
Calibrated.3RWGD.brain.jump<-Calibrated.3R.standardized.brain[c(index.trees.jump)] ##445/800=55.62% Brownian trees supported trait jump model for brain expressions
Calibrated.3RWGD.brain.nojump<-Calibrated.3R.standardized.brain[-c(index.trees.jump)] ## 355/800=44.37% trees with no support for jump in trait
rm(index.info)
rm(index.trees.jump)

### Step2: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.3Rexp.brain))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.3Rexp.brain)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 2716 obs
}

test.output.3RWGD.brain.processed<-bind_cols(test.output.3Rexp.brain,processed.fout) ## 2716 obs for all the 15 species
rm(i)
rm(test.output.3Rexp.brain)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

##Since our calibrated tree shows both the FishWGD and duplication nodes, 
##we want to modify the node events of the jump output to separate out the effect of FishWGD from the SSDs in those trees.

jump.3RWGD.brain.trees.summary<-summary.tree(Calibrated.3RWGD.brain.jump)
jump.3RWGD.brain.trees.summary$Gene<-jump.3RWGD.brain.trees.summary$index.tree
jump.3RWGD.brain.trees.summary$From<-jump.3RWGD.brain.trees.summary$node
test.output.3RWGD.brain.processed.mod<-merge(test.output.3RWGD.brain.processed,jump.3RWGD.brain.trees.summary,by=c("Gene","From"))
test.output.3RWGD.brain.processed.mod<-test.output.3RWGD.brain.processed.mod[c(1,3,2,4:11,16)]
test.output.3RWGD.brain.processed.mod$node.event<-test.output.3RWGD.brain.processed.mod$events
colnames(test.output.3RWGD.brain.processed.mod)[10]<-"pic"
test.output.3RWGD.brain.processed.mod<-test.output.3RWGD.brain.processed.mod[c(1:11)] ## 2716 obs

## Step 3: Analysis on all 15 vertebrates clades
test.output.3RWGD.brain.processed.new<-unique(test.output.3RWGD.brain.processed.mod[c(6,3,7,8,10,2,1)])  #2139 obs

########### Proportion of jump events in 15 vertebrates species in the 445 3RWGD trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.3RWGD.brain<-NULL
info.jump.3RWGD.brain<-bind_rows(lapply(Calibrated.3RWGD.brain.jump, tree.data.collection.3RWGD))
info.jump.3RWGD.brain$spe.num<-info.jump.3RWGD.brain$internal.events-(info.jump.3RWGD.brain$dup.num + info.jump.3RWGD.brain$WGD.num)
info.jump.3RWGD.brain<-info.jump.3RWGD.brain[c(11,14,3,12,5,4,13)]


## Counting frquency of jump for "speciation" and "duplication" events
test.dup.3RWGD.brain<-
  test.output.3RWGD.brain.processed.mod %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.3RWGD.brain)<-c("index.old","jump.dup")

test.WGD.3RWGD.brain<-
  test.output.3RWGD.brain.processed.mod %>%
  filter(node.event=="FishWGD" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.WGD.3RWGD.brain)<-c("index.old","jump.WGD")

test.spe.3RWGD.brain<-
  test.output.3RWGD.brain.processed.mod %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.3RWGD.brain)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.3RWGD.brain1<-full_join(test.spe.3RWGD.brain,test.dup.3RWGD.brain)
merged.3RWGD.brain<-full_join(merged.3RWGD.brain1,test.WGD.3RWGD.brain)
info.jump.3RWGD.brain.final1<-merge(info.jump.3RWGD.brain,merged.3RWGD.brain, by = c("index.old"))
info.jump.3RWGD.brain.final1[is.na(info.jump.3RWGD.brain.final1)] <- 0
rm(merged.3RWGD.brain1)

## Have to consider trees where atleast one WGD and one speciation exists
info.jump.3RWGD.brain.final<-info.jump.3RWGD.brain.final1[which(info.jump.3RWGD.brain.final1$spe.num>0 & info.jump.3RWGD.brain.final1$WGD.num>0),]
rm(info.jump.3RWGD.brain.final1)

## Proportion of jumps
info.jump.3RWGD.brain.final$dup.jump.prob<-(info.jump.3RWGD.brain.final$jump.dup/(info.jump.3RWGD.brain.final$dup.num*2))
info.jump.3RWGD.brain.final$spe.jump.prob<-(info.jump.3RWGD.brain.final$jump.spe/(info.jump.3RWGD.brain.final$spe.num*2))
info.jump.3RWGD.brain.final$WGD.jump.prob<-(info.jump.3RWGD.brain.final$jump.WGD/(info.jump.3RWGD.brain.final$WGD.num*2))

## Since,many WGD trees didnot have a duplication node, it produces NA for dup.jump.prob
## To get rid of this
info.jump.3RWGD.brain.final2<-info.jump.3RWGD.brain.final[complete.cases(info.jump.3RWGD.brain.final),] ##227 obs

## Stats
median(info.jump.3RWGD.brain.final2$dup.jump.prob) #0.25
median(info.jump.3RWGD.brain.final2$spe.jump.prob) #0.067
median(info.jump.3RWGD.brain.final2$WGD.jump.prob) #0
mean(info.jump.3RWGD.brain.final2$dup.jump.prob) #0.3478
mean(info.jump.3RWGD.brain.final2$spe.jump.prob) #0.1573
mean(info.jump.3RWGD.brain.final2$WGD.jump.prob) #0.1552

pval.jump.3RWGD.brain.spe2dup <- paired.wilcox(info.jump.3RWGD.brain.final2$spe.jump.prob,info.jump.3RWGD.brain.final2$dup.jump.prob) #p-value = 9.78e-16
pval.jump.3RWGD.brain.spe2WGD <- paired.wilcox(info.jump.3RWGD.brain.final2$spe.jump.prob,info.jump.3RWGD.brain.final2$WGD.jump.prob) #p-value = 1.31e-3

## For median data
dtfr.3RWGD.brain.15spe<-NULL
dtfr.3RWGD.brain.15spe<-data.frame(prop.jump=info.jump.3RWGD.brain.final2$spe.jump.prob,events="speciation")
dtfr.3RWGD.brain.15spe<-rbind(dtfr.3RWGD.brain.15spe,data.frame(prop.jump=info.jump.3RWGD.brain.final2$dup.jump.prob,events="duplication"))
dtfr.3RWGD.brain.15spe<-rbind(dtfr.3RWGD.brain.15spe,data.frame(prop.jump=info.jump.3RWGD.brain.final2$WGD.jump.prob,events="FishWGD"))
dtfr.3RWGD.brain.15spe$brain.abs <-abs(dtfr.3RWGD.brain.15spe$prop.jump)
dtfr.3RWGD.brain.15spe$Event<-factor(dtfr.3RWGD.brain.15spe$events, levels=c("speciation","duplication","FishWGD"))


median.data.jump.3RWGD.brain.15spe<-NULL
median.data.jump.3RWGD.brain.15spe<- dtfr.3RWGD.brain.15spe%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.3RWGD.brain.15spe$pic.round<-paste0("Median = ",round(median.data.jump.3RWGD.brain.15spe$Exp.abs,4), sep="")
median.data.jump.3RWGD.brain.15spe$count<- paste0(median.data.jump.3RWGD.brain.15spe$Event,"\n\n","(n = ",median.data.jump.3RWGD.brain.15spe$freq,")")

plot3R1<-boxplot.new.3RWGD2(dtfr.3RWGD.brain.15spe,pval.jump.3RWGD.brain.spe2dup,pval.jump.3RWGD.brain.spe2WGD,"Exp.3RWGD",median.data.jump.3RWGD.brain.15spe) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.05, xmax = 1.7, ymin = 0.67, ymax =0.67, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.4, y = 0.75, label= pval.jump.3RWGD.brain.spe2dup, fontface = 4)+
  annotate("rect", xmin = 1.05, xmax = 2.8, ymin = 0.88, ymax =0.88, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.7, y = 0.95, label= pval.jump.3RWGD.brain.spe2WGD, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.3RWGD.brain.15spe$count) ## check

# Create Data
data.jump.3RWGD.exp <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(445,355)
)

# Compute the position of labels
data.jump.3RWGD.brain <- data.jump.3RWGD.brain %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.3RWGD.brain$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.3RWGD.brain$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.3RWGD.brain<-ggplot(data.jump.3RWGD.brain, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 445 WGD trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.3RWGD.brain.jump<-summary.tree(Calibrated.3RWGD.brain.jump)
summary.3RWGD.brain.jump$Exp.abs <- abs(summary.3RWGD.brain.jump$pic)
summary.3RWGD.brain.jumps<-summary.3RWGD.brain.jump[,!(names(summary.3RWGD.brain.jump) %in% c("pic"))]
nodes.contrast.3RWGD.brain.jumps.15spe <- summary.3RWGD.brain.jumps[which(!is.na(summary.3RWGD.brain.jumps$Exp.abs)),]
nodes.contrast.3RWGD.brain.jumps.15spe$Event <- factor(nodes.contrast.3RWGD.brain.jumps.15spe$events, levels=c("speciation", "duplication", "FishWGD"))
rm(summary.3RWGD.brain.jump)
rm(summary.3RWGD.brain.jumps)

############ Testing OC using all 15 species for 445 WGD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.3RWGD.brain.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.3RWGD.brain.jumps.15spe$index<-nodes.contrast.3RWGD.brain.jumps.15spe$index.tree
test.output.3RWGD.brain.processed.mod$node<-test.output.3RWGD.brain.processed.mod$From
nodes.contrast.3RWGD.brain.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.3RWGD.brain.jumps.15spe,test.output.3RWGD.brain.processed.mod,by=c("index","node"))
nodes.contrast.3RWGD.brain.jumps.15spe.removing.jump.node <- nodes.contrast.3RWGD.brain.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.3RWGD.brain.jumps.15spe.removing.jump.node$Exp.abs)),]
nodes.contrast.3RWGD.brain.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.3RWGD.brain.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication","FishWGD"))
nodes.contrast.3RWGD.brain.jumps.15spe.removing.jump.node$pic<-abs(nodes.contrast.3RWGD.brain.jumps.15spe.removing.jump.node$Exp.abs)

########### Testing OC for all nodes of 355 3RWGD exp trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.3RWGD.brain.nojump<-summary.tree(Calibrated.3RWGD.brain.nojump)
summary.3RWGD.brain.nojump$Exp.abs <- abs(summary.3RWGD.brain.nojump$pic.brain)
summary.3RWGD.brain.nojumps<-summary.3RWGD.brain.nojump[,!(names(summary.3RWGD.brain.nojump) %in% c("pic.brain"))]
nodes.contrast.3RWGD.brain.nojumps.15spe <- summary.3RWGD.brain.nojumps[which(!is.na(summary.3RWGD.brain.nojumps$Exp.abs)),]
nodes.contrast.3RWGD.brain.nojumps.15spe$Event <- factor(nodes.contrast.3RWGD.brain.nojumps.15spe$events, levels=c("speciation", "duplication","FishWGD"))
rm(summary.3RWGD.brain.nojump)
rm(summary.3RWGD.brain.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.3RWGD.brain.new.fish<-test.output.3RWGD.brain.processed.mod[which(test.output.3RWGD.brain.processed.mod$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##7101 obs

## Step 6: Analysis on fish specific clades
test.output.3RWGD.brain.processed.new.fish<-unique(test.output.3RWGD.brain.new.fish[c(6,3,7,8,10,2,1)])  #1024obs
rm(test.output.3RWGD.brain.new.fish)

########### Proportion of jump events in teleosts in the 445 WGD trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.3RWGD.fish.brain<-NULL
info.jump.3RWGD.fish.brain<-bind_rows(lapply(Calibrated.3RWGD.brain.jump, tree.data.collection.3RWGD.fish))
info.jump.3RWGD.fish.brain$spe.num<-info.jump.3RWGD.fish.brain$internal.events.fish-(info.jump.3RWGD.fish.brain$dup.num + info.jump.3RWGD.fish.brain$WGD.num)
info.jump.3RWGD.fish.brain<-info.jump.3RWGD.fish.brain[c(11,14,3,12,5,4,13)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.3RWGD.fish.brain<-
  test.output.3RWGD.brain.new.fish %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.3RWGD.fish.brain)<-c("index.old","jump.dup")

test.WGD.3RWGD.fish.brain<-
  test.output.3RWGD.brain.new.fish %>%
  filter(node.event=="FishWGD" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.WGD.3RWGD.fish.brain)<-c("index.old","jump.WGD")

test.spe.3RWGD.fish.brain<-
  test.output.3RWGD.brain.new.fish %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.3RWGD.fish.brain)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.3RWGD.brain.fish1<-full_join(test.spe.3RWGD.fish.brain,test.dup.3RWGD.fish.brain)
merged.3RWGD.brain.fish<-full_join(merged.3RWGD.brain.fish1,test.WGD.3RWGD.fish.brain)
info.jump.brain.final.3RWGD.fish<-merge(info.jump.3RWGD.fish.brain,merged.3RWGD.brain.fish, by = c("index.old"))
info.jump.brain.final.3RWGD.fish[is.na(info.jump.brain.final.3RWGD.fish)] <- 0
rm(merged.3RWGD.brain.fish1)

## Have to consider trees where atleast one WGD and one speciation exists
info.jump.brain.final1.3RWGD.fish<-info.jump.brain.final.3RWGD.fish[which(info.jump.brain.final.3RWGD.fish$spe.num>0 & info.jump.brain.final.3RWGD.fish$WGD.num>0),]
rm(info.jump.brain.final.3RWGD.fish)

## Proportion of jumps
info.jump.brain.final1.3RWGD.fish$dup.jump.prob<-(info.jump.brain.final1.3RWGD.fish$jump.dup/(info.jump.brain.final1.3RWGD.fish$dup.num*2))
info.jump.brain.final1.3RWGD.fish$spe.jump.prob<-(info.jump.brain.final1.3RWGD.fish$jump.spe/(info.jump.brain.final1.3RWGD.fish$spe.num*2))
info.jump.brain.final1.3RWGD.fish$WGD.jump.prob<-(info.jump.brain.final1.3RWGD.fish$jump.WGD/(info.jump.brain.final1.3RWGD.fish$WGD.num*2))

## Since,many WGD trees didnot have a duplication node, it produces NA for dup.jump.prob
## To get rid of this
info.jump.brain.final2.3RWGD.fish<-info.jump.brain.final1.3RWGD.fish[complete.cases(info.jump.brain.final1.3RWGD.fish),] #97 obs

## Stats
median(na.omit(info.jump.brain.final2.3RWGD.fish$dup.jump.prob)) #0.25
median(na.omit(info.jump.brain.final2.3RWGD.fish$spe.jump.prob)) #0.167
median(na.omit(info.jump.brain.final2.3RWGD.fish$WGD.jump.prob)) #0
mean(na.omit(info.jump.brain.final2.3RWGD.fish$dup.jump.prob)) #0.321
mean(na.omit(info.jump.brain.final2.3RWGD.fish$spe.jump.prob)) #0.251
mean(na.omit(info.jump.brain.final2.3RWGD.fish$WGD.jump.prob)) #0.201

pval.jump.3RWGD.brain.fish.spe2dup <- paired.wilcox(info.jump.brain.final2.3RWGD.fish$spe.jump.prob,info.jump.brain.final2.3RWGD.fish$dup.jump.prob) #p-value = 2.33e-2
pval.jump.3RWGD.brain.fish.spe2WGD <- paired.wilcox(info.jump.brain.final2.3RWGD.fish$spe.jump.prob,info.jump.brain.final2.3RWGD.fish$WGD.jump.prob) #p-value = 2.14e-2

## For median data
dtfr.3RWGD.brain.fish<-NULL
dtfr.3RWGD.brain.fish<-data.frame(prop.jump=info.jump.brain.final2.3RWGD.fish$spe.jump.prob,events="speciation")
dtfr.3RWGD.brain.fish<-rbind(dtfr.3RWGD.brain.fish,data.frame(prop.jump=info.jump.brain.final2.3RWGD.fish$dup.jump.prob,events="duplication"))
dtfr.3RWGD.brain.fish<-rbind(dtfr.3RWGD.brain.fish,data.frame(prop.jump=info.jump.brain.final2.3RWGD.fish$WGD.jump.prob,events="FishWGD"))
dtfr.3RWGD.brain.fish$Exp.abs <-abs(dtfr.3RWGD.brain.fish$prop.jump)
dtfr.3RWGD.brain.fish$Event<-factor(dtfr.3RWGD.brain.fish$events, levels=c("speciation","duplication","FishWGD"))

library(tidyverse)
median.data.jump.3RWGD.brain.fish<-NULL
median.data.jump.3RWGD.brain.fish<- dtfr.3RWGD.brain.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.3RWGD.brain.fish$pic.round<-paste0("Median = ",round(median.data.jump.3RWGD.brain.fish$brain.abs,4), sep="")
median.data.jump.3RWGD.brain.fish$count<- paste0(median.data.jump.3RWGD.brain.fish$Event,"\n","(n = ",median.data.jump.3RWGD.brain.fish$freq,")")

plot3R2<-boxplot.new.3RWGD2(dtfr.3RWGD.brain.fish,pval.jump.3RWGD.brain.fish.spe2dup,pval.jump.3RWGD.brain.fish.spe2WGD,"Exp.3RWGD",median.data.jump.3RWGD.brain.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("Jump in exp levels for 5 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.05, xmax = 1.7, ymin = 0.67, ymax =0.67, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.4, y = 0.75, label= pval.jump.3RWGD.brain.fish.spe2dup, fontface = 4)+
  annotate("rect", xmin = 1.05, xmax = 2.8, ymin = 0.88, ymax =0.88, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.7, y = 0.95, label= pval.jump.3RWGD.brain.fish.spe2WGD, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.3RWGD.brain.fish$count) ## check

########### Testing OC for all fish specific nodes of 445 WGD trees supporting jump##################
nodes.contrast.3RWGD.brain.jumps.fish<-nodes.contrast.3RWGD.brain.jumps.15spe[which(nodes.contrast.3RWGD.brain.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 445 WGD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.3RWGD.brain.jumps.fish.removing.jump.node<-nodes.contrast.3RWGD.brain.jumps.15spe.removing.jump.node[which(nodes.contrast.3RWGD.brain.jumps.15spe.removing.jump.node$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

########### Testing OC for all fish specific nodes of 355 3RWGD trees that do not support trait jump for teleosts ##################
nodes.contrast.3RWGD.brain.nojumps.fish<-nodes.contrast.3RWGD.brain.nojumps.15spe[which(nodes.contrast.3RWGD.brain.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")

################## Chunk 5: Analysis on contrasts standardized sure SSD trees for traits ##########################

### Step 1a: Considering (sure.SSD.trees) 4139 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for brain tissue for sure SSD trees 
count<-0
Calibrated.brain.SSD1<-lapply(sure.SSD.trees, diagnostic.plot.test.exp.all.tissue, "TPM.brain")
Calibrated.brain.SSD2<-Calibrated.brain.SSD1[!is.na(Calibrated.brain.SSD1)] 
Calibrated.standardized.SSD.brain<-Calibrated.brain.SSD2[ ! sapply(Calibrated.brain.SSD2, is.null) ]## finally we obtained 2721/4139 trees passing the diagnostic tests
rm(Calibrated.brain.SSD1)
rm(Calibrated.brain.SSD2)

### Step 1b: Considering (Calibrated.standardized.SSD.brain) 2721 trees passing the diagnostic test for brain expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.standardized.SSD.brain,tree.index.collect))
rm(count)

test.output.SSD.brain<-NULL
test.output.SSD.brain<-merge(index.info,test.output.brain1,by="Gene") ## 11705 observations
index.trees.jump<-unique(test.output.SSD.brain$tree.num) ##1752 trees
Calibrated.SSD.brain.jump<-Calibrated.standardized.SSD.brain[c(index.trees.jump)] ##1752/2721=64.38% Brownian trees supported trait jump model for brain expressions
Calibrated.SSD.brain.nojump<-Calibrated.standardized.SSD.brain[-c(index.trees.jump)] ## 969/2721=35.61%  trees with no support for jump in tau
rm(index.info)
rm(index.trees.jump)

### Step2: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.SSD.brain))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.SSD.brain)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 11705 obs
}

test.output.SSD.brain.processed<-bind_cols(test.output.SSD.brain,processed.fout) ## 11705 obs for all the 15 species
test.output.SSD.brain.processed<-test.output.SSD.brain.processed[(which(!(test.output.SSD.brain.processed$clade.name %in% c("Clupeocephala","Osteoglossocephalai")))),] ##11179 obs
rm(i)
rm(test.output.SSD.brain)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

## Step 3: Analysis on all 15 vertebrates clades
test.output.SSD.brain.processed.mod<-unique(test.output.SSD.brain.processed[c(6,3,7,8,10,2,1)])  #8256 obs

########### Proportion of jump events in 15 vertebrates species in the 1752 SSD trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.SSD.brain<-NULL
count<-0
info.jump.SSD.brain<-bind_rows(lapply(Calibrated.SSD.brain.jump, tree.data.collection))
info.jump.SSD.brain$spe.num<-info.jump.SSD.brain$internal.events-info.jump.SSD.brain$dup.num
info.jump.SSD.brain<-info.jump.SSD.brain[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.SSD.brain<-
  test.output.SSD.brain.processed %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.SSD.brain)<-c("index.old","jump.dup")

test.spe.SSD.brain<-
  test.output.SSD.brain.processed %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.SSD.brain)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.SSD.brain<-full_join(test.spe.SSD.brain,test.dup.SSD.brain)
info.jump.brain.final.SSD.all<-merge(info.jump.SSD.brain,merged.SSD.brain, by = c("index.old"))
info.jump.brain.final.SSD.all[is.na(info.jump.brain.final.SSD.all)] <- 0

## Have to consider trees where atleast one duplication and speciation exists
info.jump.brain.final1.SSD.all<-info.jump.brain.final.SSD.all[which(info.jump.brain.final.SSD.all$dup.num>0 & info.jump.brain.final.SSD.all$spe.num>0),]
rm(info.jump.brain.final.SSD.all)

## Proportion of jumps
info.jump.brain.final1.SSD.all$dup.jump.prob<-(info.jump.brain.final1.SSD.all$jump.dup/(info.jump.brain.final1.SSD.all$dup.num*2))
info.jump.brain.final1.SSD.all$spe.jump.prob<-(info.jump.brain.final1.SSD.all$jump.spe/(info.jump.brain.final1.SSD.all$spe.num*2))

## Stats
median(na.omit(info.jump.brain.final1.SSD.all$dup.jump.prob)) #0.25
median(na.omit(info.jump.brain.final1.SSD.all$spe.jump.prob)) #0.071
mean(na.omit(info.jump.brain.final1.SSD.all$dup.jump.prob)) #0.3642
mean(na.omit(info.jump.brain.final1.SSD.all$spe.jump.prob)) #0.1486


########### Testing OC for all nodes of 1752 SSD trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.SSD.brain.jump<-summary.tree(Calibrated.SSD.brain.jump)
summary.SSD.brain.jump$Exp.abs <- abs(summary.SSD.brain.jump$pic.brain)
summary.SSD.brain.jumps<-summary.SSD.brain.jump[,!(names(summary.SSD.brain.jump) %in% c("pic.brain"))]
nodes.contrast.SSD.brain.jumps.15spe <- summary.SSD.brain.jumps[which(!is.na(summary.SSD.brain.jumps$Exp.abs)),]
nodes.contrast.SSD.brain.jumps.15spe$Event <- factor(nodes.contrast.SSD.brain.jumps.15spe$events, levels=c("speciation", "duplication"))
nodes.contrast.SSD.brain.jumps.15spe<-nodes.contrast.SSD.brain.jumps.15spe[(which(!(nodes.contrast.SSD.brain.jumps.15spe$label %in% c("Clupeocephala","Osteoglossocephalai")))),]
rm(summary.SSD.brain.jump)
rm(summary.SSD.brain.jumps)

############ Testing OC using all 15 species for 1752 SSD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.SSD.brain.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.SSD.brain.jumps.15spe$index<-nodes.contrast.SSD.brain.jumps.15spe$index.tree
test.output.SSD.brain.processed$node<-test.output.SSD.brain.processed$From
nodes.contrast.SSD.brain.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.SSD.brain.jumps.15spe,test.output.SSD.brain.processed,by=c("index","node"))
nodes.contrast.SSD.brain.jumps.15spe.removing.jump.node <- nodes.contrast.SSD.brain.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.SSD.brain.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.SSD.brain.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.SSD.brain.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.SSD.brain.jumps.15spe.removing.jump.node$Exp.abs<-nodes.contrast.SSD.brain.jumps.15spe.removing.jump.node$pic
nodes.contrast.SSD.brain.jumps.15spe.removing.jump.node<-nodes.contrast.SSD.brain.jumps.15spe.removing.jump.node[(which(!(nodes.contrast.SSD.brain.jumps.15spe.removing.jump.node$label %in% c("Clupeocephala","Osteoglossocephalai")))),]

########### Testing OC for all nodes of 969 SSD exp trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.SSD.brain.nojump<-summary.tree(Calibrated.SSD.brain.nojump)
summary.SSD.brain.nojump$Exp.abs <- abs(summary.SSD.brain.nojump$pic)
summary.SSD.brain.nojumps<-summary.SSD.brain.nojump[,!(names(summary.SSD.brain.nojump) %in% c("pic"))]
nodes.contrast.SSD.brain.nojumps.15spe <- summary.SSD.brain.nojumps[which(!is.na(summary.SSD.brain.nojumps$Exp.abs)),]
nodes.contrast.SSD.brain.nojumps.15spe$Event <- factor(nodes.contrast.SSD.brain.nojumps.15spe$events, levels=c("speciation", "duplication"))
nodes.contrast.SSD.brain.nojumps.15spe<-nodes.contrast.SSD.brain.nojumps.15spe[(which(!(nodes.contrast.SSD.brain.nojumps.15spe$label %in% c("Clupeocephala","Osteoglossocephalai")))),] 
rm(summary.SSD.brain.nojump)
rm(summary.SSD.brain.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.SSD.brain.new.fish<-test.output.SSD.brain.processed[which(test.output.SSD.brain.processed$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##7101 obs

## Step 6: Analysis on fish specific clades
test.output.SSD.brain.new.fish.mod<-unique(test.output.SSD.brain.new.fish[c(6,3,7,8,10,2,1)])  #2254 obs


########### Proportion of jump events in teleosts in the 1752 SSD trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.SSD.brain.trees.fish<-NULL
count<-0
info.jump.SSD.brain.trees.fish<-bind_rows(lapply(Calibrated.SSD.brain.jump, tree.data.collection.fish))
info.jump.SSD.brain.trees.fish$spe.num<-info.jump.SSD.brain.trees.fish$internal.events.fish-info.jump.SSD.brain.trees.fish$dup.num
info.jump.SSD.brain.trees.fish<-info.jump.SSD.brain.trees.fish[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.SSD.brain.fish<-
  test.output.SSD.brain.new.fish  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.SSD.brain.fish)<-c("index.old","jump.dup")


test.spe.SSD.brain.fish<-
  test.output.SSD.brain.new.fish %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.SSD.brain.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.SSD.brain.fish<-full_join(test.spe.SSD.brain.fish,test.dup.SSD.brain.fish)
info.jump.trees.SSD.brain.final.fish<-merge(info.jump.SSD.brain.trees.fish,merged.SSD.brain.fish, by = c("index.old"))
info.jump.trees.SSD.brain.final.fish[is.na(info.jump.trees.SSD.brain.final.fish)] <- 0

## Have to consider trees where atleast one fish specific duplication and speciation exists
info.jump.trees.SSD.brain.final1.fish<-info.jump.trees.SSD.brain.final.fish[which(info.jump.trees.SSD.brain.final.fish$dup.num>0 & info.jump.trees.SSD.brain.final.fish$spe.num>0),]

## Proportion of jumps
info.jump.trees.SSD.brain.final1.fish$dup.jump.prob<-(info.jump.trees.SSD.brain.final1.fish$jump.dup/(info.jump.trees.SSD.brain.final1.fish$dup.num*2))
info.jump.trees.SSD.brain.final1.fish$spe.jump.prob<-(info.jump.trees.SSD.brain.final1.fish$jump.spe/(info.jump.trees.SSD.brain.final1.fish$spe.num*2))


median(na.omit(info.jump.trees.SSD.brain.final1.fish$dup.jump.prob)) #0.5
median(na.omit(info.jump.trees.SSD.brain.final1.fish$spe.jump.prob)) #0.125
mean(na.omit(info.jump.trees.SSD.brain.final1.fish$dup.jump.prob)) #0.3818
mean(na.omit(info.jump.trees.SSD.brain.final1.fish$spe.jump.prob)) #0.1993

########### Testing OC for all fish specific nodes of 1752 SSD trees supporting jump##################
nodes.contrast.SSD.brain.jumps.fish<-nodes.contrast.SSD.brain.jumps.15spe[which(nodes.contrast.SSD.brain.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 1752 SSD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 

nodes.contrast.SSD.brain.jumps.fish.removing.jump.node <- NULL
nodes.contrast.SSD.brain.jumps.fish$index<-nodes.contrast.SSD.brain.jumps.fish$index.tree
test.output.SSD.brain.new.fish$node<-test.output.SSD.brain.new.fish$From
nodes.contrast.SSD.brain.jumps.fish.removing.jump.node <-anti_join(nodes.contrast.SSD.brain.jumps.fish,test.output.SSD.brain.new.fish,by=c("index","node"))
nodes.contrast.SSD.brain.jumps.fish.removing.jump.node <- nodes.contrast.SSD.brain.jumps.fish.removing.jump.node[which(!is.na(nodes.contrast.SSD.brain.jumps.fish.removing.jump.node$pic)),]
nodes.contrast.SSD.brain.jumps.fish.removing.jump.node$Event <- factor(nodes.contrast.SSD.brain.jumps.fish.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.SSD.brain.jumps.fish.removing.jump.node$Exp.abs<-nodes.contrast.SSD.brain.jumps.fish.removing.jump.node$pic
nodes.contrast.SSD.brain.jumps.fish.removing.jump.node<-nodes.contrast.SSD.brain.jumps.fish.removing.jump.node[(which(!(nodes.contrast.SSD.brain.jumps.fish.removing.jump.node$label %in% c("Clupeocephala","Osteoglossocephalai")))),]

########### Testing OC for all fish specific nodes of 969 SSD exp trees that do not support trait jump for teleosts ##################
nodes.contrast.SSD.brain.nojumps.fish<-nodes.contrast.SSD.brain.nojumps.15spe[which(nodes.contrast.SSD.brain.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")



