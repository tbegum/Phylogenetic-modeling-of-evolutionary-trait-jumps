

### This code is written to work with "trait jump model" for tau
## This code has a few modifications since the tree sets were given to Pablo, somehow it's index was different than the later version
## So we needed to match index before processing the data

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
load("norm_exp_15spe_4tips_new.RData")  ## latest modified tau 

######################### Chunk1: Data given to pablo for applying "Levolution" to find trait jump #############################
## For Tau "trees.all.interestt.old" of 7711 gene trees were given to pablo
## So the output index was according to this tree

######################### Chunk2: reading output of trait jump #############################

test.output.tau<-read.csv(paste0(folder1,"To_pablo/Tau/trait_jump_output_tau.csv", sep=""),sep=",",header=T) ##201206 obs

## Only considering trees with posterior probability of jump more than 70%
test.output1.tau<-test.output.tau[which(test.output.tau$JumpProbability>=0.7),] ##36936 obs

## This code has a few modifications since the tree sets were given to Pablo, somehow it's index was different than the later version
## So we needed to match indices before processing the data
## We had indices for trees.all.interest (6923) trees
## So we need to match first the indices of trees.all.interestt.old (7711) trees to trees.all.interest (6923) trees
index.map<-data.frame(index.6923tree.tau=NA,index.7711tree.tau=NA)

for(j in 1:length(trees.all.interest))
{
  cat(j,"\t")
  for(k in 1:length(trees.all.interestt.old))
  {
    if(trees.all.interest[[j]]@phylo$tip.label %in% trees.all.interestt.old[[k]]@phylo$tip.label)
    {
      index.map<-rbind(index.map,data.frame(index.6923tree.tau=j,index.7711tree.tau=k))
      cat(k,"\n")
      break
    }
  }
}
rm(j)
rm(k)

index.map<-index.map[-1,] ## 6703 trees mapped
index.map$Gene<-index.map$index.7711tree.tau


########################## Chunk 3: Analysis on all contrasts standardized trees for traits ###############################

### Step 1: Considering (Calibrated.all.standardized.tau) 4247 tree data passing the diagnostic test for tau
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.all.standardized.tau,tree.index.collect))
rm(count)

## We need to match this indices with the previous one
colnames(index.info)[2]<-colnames(index.map)[1]
index.info<-merge(index.info,index.map,by=c("index.6923tree.tau"))


test.output.tau.all<-NULL
test.output.tau.all<-merge(index.info,test.output1.tau,by="Gene") ## 14771 observations
index.trees.jump<-unique(test.output.tau.all$tree.num) ##2636 trees
Calibrated.standardized.tau.jump<-Calibrated.all.standardized.tau[c(index.trees.jump)] ##2636/4247=62% Brownian trees supported trait jump model for tau
Calibrated.standardized.tau.nojump<-Calibrated.all.standardized.tau[-c(index.trees.jump)] ## 1611 trees with no support for jump in tau
#rm(index.info)
rm(index.trees.jump)

### Step2: Calculating PICs of required trait for trees to match the index number of the trees on which levolution was run
Count<-0
trees.all.interest<-lapply(trees.all.interest,contrast.calc,trait="Tau")

Count<-0
Calibrated.standardized.tau.jump<-lapply(Calibrated.standardized.tau.jump,contrast.calc,trait="Tau") 

rm(Count)

### Step3: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.tau.all))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output(test.output.tau.all)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 14771 obs
}

test.output.tau.processed<-bind_cols(test.output.tau.all,processed.fout) ## 14771 obs for all the 15 species
rm(i)
#rm(test.output.tau.all)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

## Is OC supported for tau for vertebrate clades supporting jump?
pic.spe.all<-test.output.tau.processed$pic[which(test.output.tau.processed$node.event=="speciation")]
pic.dup.all<-test.output.tau.processed$pic[which(test.output.tau.processed$node.event=="duplication")]
median(pic.spe.all) ##0.00998
median(pic.dup.all) ##0.0118
Pval.all.jump.all<-two.tailed.wilcox(pic.spe.all,pic.dup.all) # p= 4.65e-16

## Step 4: Analysis on all 15 vertebrates clades
test.output.tau.new.mod<-unique(test.output.tau.processed[c(8,5,9,10,12,4,3)])  #11306 obs

########### Proportion of jump events in 15 vertebrates species in the 2636 trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.trees.tau<-NULL
count<-0
info.jump.trees.tau<-bind_rows(lapply(Calibrated.standardized.tau.jump, tree.data.collection))
info.jump.trees.tau$spe.num<-info.jump.trees.tau$internal.events-info.jump.trees.tau$dup.num
info.jump.trees.tau<-info.jump.trees.tau[c(11,12,3:5)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.tau.all<-
  test.output.tau.processed  %>%
  filter(node.event=="duplication" ) %>%
  .$index %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.tau.all)<-c("index.old","jump.dup")

test.spe.tau.all<-
  test.output.tau.processed%>%
  filter(node.event=="speciation") %>%
  .$index %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.tau.all)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.tau.all<-full_join(test.spe.tau.all,test.dup.tau.all)
info.jump.trees.tau.final.all<-merge(info.jump.trees.tau,merged.tau.all, by = c("index.old"))
info.jump.trees.tau.final.all[is.na(info.jump.trees.tau.final.all)] <- 0

## Have to consider trees where atleast one duplication and speciation exists
info.jump.trees.tau.final1.all<-info.jump.trees.tau.final.all[which(info.jump.trees.tau.final.all$dup.num>0 & info.jump.trees.tau.final.all$spe.num>0),]
rm(info.jump.trees.tau.final.all)

## Proportion of jumps
info.jump.trees.tau.final1.all$dup.jump.prob<-(info.jump.trees.tau.final1.all$jump.dup/(info.jump.trees.tau.final1.all$dup.num*2))
info.jump.trees.tau.final1.all$spe.jump.prob<-(info.jump.trees.tau.final1.all$jump.spe/(info.jump.trees.tau.final1.all$spe.num*2))

## Stats
median(na.omit(info.jump.trees.tau.final1.all$dup.jump.prob)) #0.1667
median(na.omit(info.jump.trees.tau.final1.all$spe.jump.prob)) #0.0555
mean(na.omit(info.jump.trees.tau.final1.all$dup.jump.prob)) #0.2624
mean(na.omit(info.jump.trees.tau.final1.all$spe.jump.prob)) #0.1414

pval.jump.tau.all <- paired.wilcox(info.jump.trees.tau.final1.all$spe.jump.prob,info.jump.trees.tau.final1.all$dup.jump.prob) #p-value = 2.2e-16


plot4A<-boxplot.new2(dtfr.tau.all,pval.jump.tau.all,"Tau",median.data.jump.tau) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste(tau,", for 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.tau.all, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.tau$count) ## check

# Create Data
data.jump.tau <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(2636,1611)
)

# Compute the position of labels
data.jump.tau <- data.jump.tau %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.tau$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.tau$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.tau<-ggplot(data.jump.tau, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 2636 trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized tau levels 
summary.tau.jump<-summary.tree(Calibrated.standardized.tau.jump)
summary.tau.jump$Tau.abs <- abs(summary.tau.jump$pic_Tau)
summary.tau.jumps<-summary.tau.jump[,!(names(summary.tau.jump) %in% c("pic_Tau"))]
nodes.contrast.tau.jumps.15spe <- summary.tau.jumps[which(!is.na(summary.tau.jumps$Tau.abs)),]
nodes.contrast.tau.jumps.15spe$Event <- factor(nodes.contrast.tau.jumps.15spe$events, levels=c("speciation", "duplication")) ##48476 obs
rm(summary.tau.jump)
rm(summary.tau.jumps)

############ Testing OC using all 15 species for 2636 trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
library(dplyr)
nodes.contrast.tau.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.tau.jumps.15spe$index<-nodes.contrast.tau.jumps.15spe$index.tree
test.output.tau.processed$node<-test.output.tau.processed$From
nodes.contrast.tau.jumps.15spe.removing.jump.node <- anti_join(nodes.contrast.tau.jumps.15spe,test.output.tau.processed,by=c("index","node"))
nodes.contrast.tau.jumps.15spe.removing.jump.node <- nodes.contrast.tau.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.tau.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.tau.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.tau.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication")) #37170 obs
nodes.contrast.tau.jumps.15spe.removing.jump.node$Tau.abs<-nodes.contrast.tau.jumps.15spe.removing.jump.node$pic

pic.spe.removing.jump.tau.15species.2636trees<-speciation.contrast(nodes.contrast.tau.jumps.15spe.removing.jump.node,"pic") #(n=32894)
pic.dup.removing.jump.tau.15species.2636trees<-duplication.contrast(nodes.contrast.tau.jumps.15spe.removing.jump.node,"pic") #(n=4276)


########### Testing OC for all nodes of 1611 trees that do not support trait jump for tau in 15 species ##################
## Summarizing contrasts for normalized tau levels 
summary.tau.nojump<-summary.tree(Calibrated.standardized.tau.nojump)
summary.tau.nojump$Tau.abs <- abs(summary.tau.nojump$pic_Tau)
summary.tau.nojumps<-summary.tau.nojump[,!(names(summary.tau.nojump) %in% c("pic_Tau"))]
nodes.contrast.tau.nojumps.15spe <- summary.tau.nojump[which(!is.na(summary.tau.nojump$Tau.abs)),]
nodes.contrast.tau.nojumps.15spe$Event <- factor(nodes.contrast.tau.nojumps.15spe$events, levels=c("speciation", "duplication"))
rm(summary.tau.nojump)
rm(summary.tau.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.tau.new.fish<-test.output.tau.processed[which(test.output.tau.processed$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##5984 obs

## Is OC supported for tau for fish clades supporting jump?
pic.spe.fish<-test.output.tau.new.fish$pic[which(test.output.tau.new.fish$node.event=="speciation")]
pic.dup.fish<-test.output.tau.new.fish$pic[which(test.output.tau.new.fish$node.event=="duplication")]
median(pic.spe.fish) ##0.0088
median(pic.dup.fish) ##0.018
Pval.all.jump.fish<-two.tailed.wilcox(pic.spe.fish,pic.dup.fish) # < 2.21e-16

## Step 6: Analysis on fish specific clades
test.output.tau.new.fish.mod<-unique(test.output.tau.new.fish[c(8,5,9,10,12,4,3)])  #4516 obs
rm(test.output.tau.new.fish)

########### Proportion of jump events in teleosts in the 2636 trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.trees.tau.fish<-NULL
count<-0
info.jump.trees.tau.fish<-bind_rows(lapply(Calibrated.standardized.tau.jump, tree.data.collection.fish))
info.jump.trees.tau.fish$spe.num<-info.jump.trees.tau.fish$internal.events.fish-info.jump.trees.tau.fish$dup.num
info.jump.trees.tau.fish<-info.jump.trees.tau.fish[c(11,12,3:5)] #2636 obs
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.tau.fish<-
  test.output.tau.new.fish  %>%
  filter(node.event=="duplication" ) %>%
  .$index %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.tau.fish)<-c("index.old","jump.dup")

test.spe.tau.fish<-
  test.output.tau.new.fish%>%
  filter(node.event=="speciation") %>%
  .$index %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.tau.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
#test.spe.tau.fish$index.old<-as.character(test.spe.tau.fish$index.old)
#test.dup.tau.fish$index.old<-as.character(test.dup.tau.fish$index.old)
merged.tau.fish<-full_join(test.spe.tau.fish,test.dup.tau.fish)
info.jump.trees.tau.final.fish<-merge(info.jump.trees.tau.fish,merged.tau.fish, by = c("index.old"))
info.jump.trees.tau.final.fish[is.na(info.jump.trees.tau.final.fish)] <- 0


## Have to consider trees where atleast one fish specific duplication and speciation exists
info.jump.trees.tau.final1.fish<-info.jump.trees.tau.final.fish[which(info.jump.trees.tau.final.fish$dup.num>0 & info.jump.trees.tau.final.fish$spe.num>0),]

## Proportion of jumps
info.jump.trees.tau.final1.fish$dup.jump.prob<-(info.jump.trees.tau.final1.fish$jump.dup/(info.jump.trees.tau.final1.fish$dup.num*2))
info.jump.trees.tau.final1.fish$spe.jump.prob<-(info.jump.trees.tau.final1.fish$jump.spe/(info.jump.trees.tau.final1.fish$spe.num*2))


median(na.omit(info.jump.trees.tau.final1.fish$dup.jump.prob)) #0.1667
median(na.omit(info.jump.trees.tau.final1.fish$spe.jump.prob)) #0.125
mean(na.omit(info.jump.trees.tau.final1.fish$dup.jump.prob)) #0.2622
mean(na.omit(info.jump.trees.tau.final1.fish$spe.jump.prob)) #0.32210


plot4B<-boxplot.new2(dtfr.tau.fish,pval.jump.tau.fish,"Tau",median.data.jump.tau.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste(tau,", for 5 teleosts"))))+
  #ylab(expression(bold(paste("Proportions of jump events for teleosts ")))) +
  ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.tau.fish, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.tau.fish$count) ## check

########### Testing OC for all fish specific nodes of 2636 trees supporting jump##################
nodes.contrast.tau.jumps.fish<-nodes.contrast.tau.jumps.15spe[which(nodes.contrast.tau.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 2636 trees by removing nodes supporting trait jump ###############
nodes.contrast.tau.jumps.fish.removing.jump.node<-nodes.contrast.tau.jumps.15spe.removing.jump.node[which(nodes.contrast.tau.jumps.15spe.removing.jump.node$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]


##Removing nodes supporting trait jump 
library(dplyr)
nodes.contrast.tau.jumps.fish.removing.jump.node <- NULL
nodes.contrast.tau.jumps.fish$index<-nodes.contrast.tau.jumps.fish$index.tree
test.output.tau.new.fish$node<-test.output.tau.new.fish$From
nodes.contrast.tau.jumps.fish.removing.jump.node <-anti_join(nodes.contrast.tau.jumps.fish,test.output.tau.new.fish,by=c("index","node"))
nodes.contrast.tau.jumps.fish.removing.jump.node <- nodes.contrast.tau.jumps.fish.removing.jump.node[which(!is.na(nodes.contrast.tau.jumps.fish.removing.jump.node$pic)),]
nodes.contrast.tau.jumps.fish.removing.jump.node$Event <- factor(nodes.contrast.tau.jumps.fish.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.tau.jumps.fish.removing.jump.node$Tau.abs<-nodes.contrast.tau.jumps.fish.removing.jump.node$pic

pic.spe.removing.jump.tau.fish.2636trees<-speciation.contrast(nodes.contrast.tau.jumps.fish.removing.jump.node,"pic") #(n=9863)
pic.dup.removing.jump.tau.fish.2636trees<-duplication.contrast(nodes.contrast.tau.jumps.fish.removing.jump.node,"pic") #(n=2282)


########### Testing OC for all fish specific nodes of 1611 trees that do not support trait jump for teleosts ##################
nodes.contrast.tau.nojumps.fish<-nodes.contrast.tau.nojumps.15spe[which(nodes.contrast.tau.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")

################## Chunk 4: Analysis on contrasts standardized sure 3RWGD trees for traits ##########################

### Step 1: Considering (Calibrated.3R.tau) 746 3RWGD trees data passing the diagnostic test for average expression levels 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info1<-NULL
count<-0
index.info1<-bind_rows(lapply(Calibrated.3R.tau,tree.index.collect))
rm(count)

## We need to match this indices with the previous one
colnames(index.info1)[2]<-colnames(index.map)[1]
index.info1<-merge(index.info1,index.map,by=c("index.6923tree.tau"))


test.output.tau.all<-NULL
test.output.tau.all<-merge(index.info1,test.output1.tau,by="Gene") ## 2216 observations
index.trees.jump<-unique(test.output.tau.all$tree.num) ##435 trees
Calibrated.3RWGD.tau.jump<-Calibrated.3R.tau[c(index.trees.jump)] ##435/746=58.31% Brownian trees supported trait jump model for average expression levels
Calibrated.3RWGD.tau.nojump<-Calibrated.3R.tau[-c(index.trees.jump)] ## 311 trees with no support for jump in tau
rm(index.info1)
rm(index.trees.jump)

### Step2: Calculating PICs of required trait for trees to match the index number of the trees on which levolution was run
Count<-0
trees.all.interest<-lapply(trees.all.interest,contrast.calc,trait="Tau")

Count<-0
Calibrated.3RWGD.tau.jump<-lapply(Calibrated.3RWGD.tau.jump,contrast.calc,trait="Tau") 
rm(Count)

### Step3: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.tau.all))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output(test.output.tau.all)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 2216 obs
}

test.output.3RWGD.tau.processed<-bind_cols(test.output.tau.all,processed.fout) ## 2216 obs for all the 15 species
rm(i)
rm(test.output.tau.all)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

##Since our calibrated tree shows both the FishWGD and duplication nodes, 
##we want to modify the node events of the jump output to separate out the effect of FishWGD from the SSDs in those trees.

jump.3RWGD.tautrees.summary<-summary.tree(Calibrated.3RWGD.tau.jump)
jump.3RWGD.tautrees.summary$index<-jump.3RWGD.tautrees.summary$index.tree
jump.3RWGD.tautrees.summary$From<-jump.3RWGD.tautrees.summary$node
test.output.3RWGD.tau.processed$node<-test.output.3RWGD.tau.processed$From
test.output.3RWGD.tau.processed.mod<-merge(test.output.3RWGD.tau.processed,jump.3RWGD.tautrees.summary,by=c("index","From"))
test.output.3RWGD.tau.processed.mod<-test.output.3RWGD.tau.processed.mod[c(1,6,5,2,7:8,10,12,19)]
test.output.3RWGD.tau.processed.mod$node.event<-test.output.3RWGD.tau.processed.mod$events
test.output.3RWGD.tau.processed.mod$node<-test.output.3RWGD.tau.processed.mod$From ## 2216 obs

## Is OC supported for vertebrate clades of 3RWGD tau trees that are supporting jump?
## Answer is no
pic.spe.3RWGD.tau<-test.output.3RWGD.tau.processed.mod$pic[which(test.output.3RWGD.tau.processed.mod$node.event=="speciation")]
pic.WGD.3RWGD.tau<-test.output.3RWGD.tau.processed.mod$pic[which(test.output.3RWGD.tau.processed.mod$node.event=="FishWGD")]
pic.SSD.3RWGD.tau<-test.output.3RWGD.tau.processed.mod$pic[which(test.output.3RWGD.tau.processed.mod$node.event=="duplication")]

## Step 4: Analysis on all 15 vertebrates clades
test.output.3RWGD.tau.processed.tau.new.mod<-unique(test.output.3RWGD.tau.processed.mod[c(1,4,10,7,8,3,2)])  #1783 obs

########### Proportion of jump events in 15 vertebrates species in the 435 3RWGD trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.3RWGD.tau.trees<-NULL
info.jump.3RWGD.tau.trees<-bind_rows(lapply(Calibrated.3RWGD.tau.jump, tree.data.collection.3RWGD))
info.jump.3RWGD.tau.trees$spe.num<-info.jump.3RWGD.tau.trees$internal.events-(info.jump.3RWGD.tau.trees$dup.num + info.jump.3RWGD.tau.trees$WGD.num)
info.jump.3RWGD.tau.trees<-info.jump.3RWGD.tau.trees[c(11,14,3,12,5,4,13)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.3RWGD.tau.all<-
  test.output.3RWGD.tau.processed.mod %>%
  filter(node.event=="duplication" ) %>%
  .$index %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.3RWGD.tau.all)<-c("index.old","jump.dup")

test.WGD.3RWGD.tau.all<-
  test.output.3RWGD.tau.processed.mod %>%
  filter(node.event=="FishWGD" ) %>%
  .$index %>% 
  table() %>% 
  as.data.frame()
colnames(test.WGD.3RWGD.tau.all)<-c("index.old","jump.WGD")

test.spe.3RWGD.tau.all<-
  test.output.3RWGD.tau.processed.mod%>%
  filter(node.event=="speciation") %>%
  .$index %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.3RWGD.tau.all)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.3RWGD.tau.all1<-full_join(test.spe.3RWGD.tau.all,test.dup.3RWGD.tau.all)
merged.3RWGD.tau.all<-full_join(merged.3RWGD.tau.all1,test.WGD.3RWGD.tau.all)
info.jump.trees.final.3RWGD.tau.all<-merge(info.jump.3RWGD.tau.trees,merged.3RWGD.tau.all, by = c("index.old"))
info.jump.trees.final.3RWGD.tau.all[is.na(info.jump.trees.final.3RWGD.tau.all)] <- 0
rm(merged.3RWGD.tau.all1)

## Have to consider trees where atleast one WGD and one speciation exists
info.jump.trees.final1.3RWGD.tau.all<-info.jump.trees.final.3RWGD.tau.all[which(info.jump.trees.final.3RWGD.tau.all$spe.num>0 & info.jump.trees.final.3RWGD.tau.all$WGD.num>0),]
rm(info.jump.trees.final.3RWGD.tau.all)

## Proportion of jumps
info.jump.trees.final1.3RWGD.tau.all$dup.jump.prob<-(info.jump.trees.final1.3RWGD.tau.all$jump.dup/(info.jump.trees.final1.3RWGD.tau.all$dup.num*2))
info.jump.trees.final1.3RWGD.tau.all$spe.jump.prob<-(info.jump.trees.final1.3RWGD.tau.all$jump.spe/(info.jump.trees.final1.3RWGD.tau.all$spe.num*2))
info.jump.trees.final1.3RWGD.tau.all$WGD.jump.prob<-(info.jump.trees.final1.3RWGD.tau.all$jump.WGD/(info.jump.trees.final1.3RWGD.tau.all$WGD.num*2))

## Since,many WGD trees didnot have a duplication node, it produces NA for dup.jump.prob
## To get rid of this
info.jump.trees.final2.3RWGD.tau.all<-info.jump.trees.final1.3RWGD.tau.all[complete.cases(info.jump.trees.final1.3RWGD.tau.all),]

## Stats
median(na.omit(info.jump.trees.final2.3RWGD.tau.all$dup.jump.prob)) #0.125
median(na.omit(info.jump.trees.final2.3RWGD.tau.all$spe.jump.prob)) #0.042
median(na.omit(info.jump.trees.final2.3RWGD.tau.all$WGD.jump.prob)) #0
mean(na.omit(info.jump.trees.final2.3RWGD.tau.all$dup.jump.prob)) #0.280
mean(na.omit(info.jump.trees.final2.3RWGD.tau.all$spe.jump.prob)) #0.121
mean(na.omit(info.jump.trees.final2.3RWGD.tau.all$WGD.jump.prob)) #0.125

pval.jump.3RWGD.tau.all.spe2dup <- paired.wilcox(info.jump.trees.final2.3RWGD.tau.all$spe.jump.prob,info.jump.trees.final2.3RWGD.tau.all$dup.jump.prob) #p-value = 2.77e-11
pval.jump.3RWGD.tau.all.spe2WGD <- paired.wilcox(info.jump.trees.final2.3RWGD.tau.all$spe.jump.prob,info.jump.trees.final2.3RWGD.tau.all$WGD.jump.prob) #p-value = 5.74e-5

## For median data
dtfr.3RWGD.tau.15spe<-NULL
dtfr.3RWGD.tau.15spe<-data.frame(prop.jump=info.jump.trees.final2.3RWGD.tau.all$spe.jump.prob,events="speciation")
dtfr.3RWGD.tau.15spe<-rbind(dtfr.3RWGD.tau.15spe,data.frame(prop.jump=info.jump.trees.final2.3RWGD.tau.all$dup.jump.prob,events="duplication"))
dtfr.3RWGD.tau.15spe<-rbind(dtfr.3RWGD.tau.15spe,data.frame(prop.jump=info.jump.trees.final2.3RWGD.tau.all$WGD.jump.prob,events="FishWGD"))
dtfr.3RWGD.tau.15spe$Tau.abs <-abs(dtfr.3RWGD.tau.15spe$prop.jump)
dtfr.3RWGD.tau.15spe$Event<-factor(dtfr.3RWGD.tau.15spe$events, levels=c("speciation","duplication","FishWGD"))

library(tidyverse)
median.data.jump.3RWGD.tau.15spe<-NULL
median.data.jump.3RWGD.tau.15spe<- dtfr.3RWGD.tau.15spe%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Tau.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.3RWGD.tau.15spe$pic.round<-paste0("Median = ",round(median.data.jump.3RWGD.tau.15spe$Tau.abs,4), sep="")
median.data.jump.3RWGD.tau.15spe$count<- paste0(median.data.jump.3RWGD.tau.15spe$Event,"\n\n","(n = ",median.data.jump.3RWGD.tau.15spe$freq,")")

plot3R1<-boxplot.new.3RWGD2(dtfr.3RWGD.tau.15spe,pval.jump.3RWGD.tau.all.spe2dup,pval.jump.3RWGD.tau.all.spe2WGD,"Tau.3RWGD",median.data.jump.3RWGD.tau.15spe) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.05, xmax = 1.7, ymin = 0.67, ymax =0.67, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.4, y = 0.75, label= pval.jump.3RWGD.tau.all.spe2dup, fontface = 4)+
  annotate("rect", xmin = 1.05, xmax = 2.8, ymin = 0.88, ymax =0.88, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.7, y = 0.95, label= pval.jump.3RWGD.tau.all.spe2WGD, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.3RWGD.tau.15spe$count) ## check

# Create Data
data.jump.3RWGD.tau <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(435,311)
)

# Compute the position of labels
data.jump.3RWGD.tau <- data.jump.3RWGD.tau %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.3RWGD.tau$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.3RWGD.tau$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.3RWGD.tau<-ggplot(data.jump.3RWGD.tau, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 435 WGD trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized tau levels 
summary.3RWGD.tau.jump<-summary.tree(Calibrated.3RWGD.tau.jump)
summary.3RWGD.tau.jump$Tau.abs <- abs(summary.3RWGD.tau.jump$pic_Tau)
summary.3RWGD.tau.jumps<-summary.3RWGD.tau.jump[,!(names(summary.3RWGD.tau.jump) %in% c("pic_Tau"))]
nodes.contrast.3RWGD.tau.jumps.15spe <- summary.3RWGD.tau.jumps[which(!is.na(summary.3RWGD.tau.jumps$Tau.abs)),]
nodes.contrast.3RWGD.tau.jumps.15spe$Event <- factor(nodes.contrast.3RWGD.tau.jumps.15spe$events, levels=c("speciation", "duplication", "FishWGD"))
rm(summary.3RWGD.tau.jump)
rm(summary.3RWGD.tau.jumps)

############ Testing OC using all 15 species for 435 WGD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.3RWGD.tau.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.3RWGD.tau.jumps.15spe$index<-nodes.contrast.3RWGD.tau.jumps.15spe$index.tree
test.output.3RWGD.tau.processed.mod$node<-test.output.3RWGD.tau.processed.mod$From
nodes.contrast.3RWGD.tau.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.3RWGD.tau.jumps.15spe,test.output.3RWGD.tau.processed.mod,by=c("index","node"))
nodes.contrast.3RWGD.tau.jumps.15spe.removing.jump.node <- nodes.contrast.3RWGD.tau.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.3RWGD.tau.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.3RWGD.tau.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.3RWGD.tau.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication","FishWGD"))
nodes.contrast.3RWGD.tau.jumps.15spe.removing.jump.node$Tau.abs<-abs(nodes.contrast.3RWGD.tau.jumps.15spe.removing.jump.node$pic)

pic.spe.removing.jump.15species.3RWGD.tau.trees<-speciation.contrast(nodes.contrast.3RWGD.tau.jumps.15spe.removing.jump.node,"pic") 
pic.WGD.removing.jump.15species.3RWGD.tau.trees<-duplication.contrast.3R(nodes.contrast.3RWGD.tau.jumps.15spe.removing.jump.node,"pic")
pic.dup.removing.jump.15species.3RWGD.tau.trees<-duplication.contrast(nodes.contrast.3RWGD.tau.jumps.15spe.removing.jump.node,"pic")


########### Testing OC for all nodes of 311 3RWGD tau trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized tau levels 
summary.3RWGD.tau.nojump<-summary.tree(Calibrated.3RWGD.tau.nojump)
summary.3RWGD.tau.nojump$Tau.abs <- abs(summary.3RWGD.tau.nojump$pic_Tau)
summary.3RWGD.tau.nojumps<-summary.3RWGD.tau.nojump[,!(names(summary.3RWGD.tau.nojump) %in% c("pic_Tau"))]
nodes.contrast.3RWGD.tau.nojumps.15spe <- summary.3RWGD.tau.nojumps[which(!is.na(summary.3RWGD.tau.nojumps$Tau.abs)),]
nodes.contrast.3RWGD.tau.nojumps.15spe$Event <- factor(nodes.contrast.3RWGD.tau.nojumps.15spe$events, levels=c("speciation", "duplication","FishWGD"))
rm(summary.3RWGD.tau.nojump)
rm(summary.3RWGD.tau.nojumps)


### Step 5: Identifying teleosts fish specific clades
test.output.3RWGD.tau.new.fish<-test.output.3RWGD.tau.processed.mod[which(test.output.3RWGD.tau.processed.mod$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##7101 obs

## Step 6: Analysis on fish specific clades
test.output.3RWGD.tau.processed.fish.mod<-unique(test.output.3RWGD.tau.new.fish[c(1,4,10,7,8,3,2)])  #910 obs
rm(test.output.3RWGD.tau.new.fish)

########### Proportion of jump events in teleosts in the 435 WGD trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.3RWGD.tau.fish.trees<-NULL
info.jump.3RWGD.tau.fish.trees<-bind_rows(lapply(Calibrated.3RWGD.tau.jump, tree.data.collection.3RWGD.fish))
info.jump.3RWGD.tau.fish.trees$spe.num<-info.jump.3RWGD.tau.fish.trees$internal.events.fish-(info.jump.3RWGD.tau.fish.trees$dup.num + info.jump.3RWGD.tau.fish.trees$WGD.num)
info.jump.3RWGD.tau.fish.trees<-info.jump.3RWGD.tau.fish.trees[c(11,14,3,12,5,4,13)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.3RWGD.tau.fish<-
  test.output.3RWGD.tau.new.fish %>%
  filter(node.event=="duplication" ) %>%
  .$index %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.3RWGD.tau.fish)<-c("index.old","jump.dup")

test.WGD.3RWGD.tau.fish<-
  test.output.3RWGD.tau.new.fish %>%
  filter(node.event=="FishWGD" ) %>%
  .$index %>% 
  table() %>% 
  as.data.frame()
colnames(test.WGD.3RWGD.tau.fish)<-c("index.old","jump.WGD")

test.spe.3RWGD.tau.fish<-
  test.output.3RWGD.tau.new.fish%>%
  filter(node.event=="speciation") %>%
  .$index %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.3RWGD.tau.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.3RWGD.tau.fish1<-full_join(test.spe.3RWGD.tau.fish,test.dup.3RWGD.tau.fish)
merged.3RWGD.tau.fish<-full_join(merged.3RWGD.tau.fish1,test.WGD.3RWGD.tau.fish)
info.jump.trees.final.3RWGD.tau.fish<-merge(info.jump.3RWGD.tau.fish.trees,merged.3RWGD.tau.fish, by = c("index.old"))
info.jump.trees.final.3RWGD.tau.fish[is.na(info.jump.trees.final.3RWGD.tau.fish)] <- 0
rm(merged.3RWGD.tau.fish1)

## Have to consider trees where atleast one WGD and one speciation exists
info.jump.trees.final1.3RWGD.tau.fish<-info.jump.trees.final.3RWGD.tau.fish[which(info.jump.trees.final.3RWGD.tau.fish$spe.num>0 & info.jump.trees.final.3RWGD.tau.fish$WGD.num>0),]
rm(info.jump.trees.final.3RWGD.tau.fish)

## Proportion of jumps
info.jump.trees.final1.3RWGD.tau.fish$dup.jump.prob<-(info.jump.trees.final1.3RWGD.tau.fish$jump.dup/(info.jump.trees.final1.3RWGD.tau.fish$dup.num*2))
info.jump.trees.final1.3RWGD.tau.fish$spe.jump.prob<-(info.jump.trees.final1.3RWGD.tau.fish$jump.spe/(info.jump.trees.final1.3RWGD.tau.fish$spe.num*2))
info.jump.trees.final1.3RWGD.tau.fish$WGD.jump.prob<-(info.jump.trees.final1.3RWGD.tau.fish$jump.WGD/(info.jump.trees.final1.3RWGD.tau.fish$WGD.num*2))

## Since,many WGD trees didnot have a duplication node, it produces NA for dup.jump.prob
## To get rid of this
info.jump.trees.final2.3RWGD.tau.fish<-info.jump.trees.final1.3RWGD.tau.fish[complete.cases(info.jump.trees.final1.3RWGD.tau.fish),]

## Stats
median(na.omit(info.jump.trees.final2.3RWGD.tau.fish$dup.jump.prob)) #0.098
median(na.omit(info.jump.trees.final2.3RWGD.tau.fish$spe.jump.prob)) #0.083
median(na.omit(info.jump.trees.final2.3RWGD.tau.fish$WGD.jump.prob)) #0
mean(na.omit(info.jump.trees.final2.3RWGD.tau.fish$dup.jump.prob)) #0.211
mean(na.omit(info.jump.trees.final2.3RWGD.tau.fish$spe.jump.prob)) #0.1398
mean(na.omit(info.jump.trees.final2.3RWGD.tau.fish$WGD.jump.prob)) #0.083

pval.jump.3RWGD.tau.fish.spe2dup <- paired.wilcox(info.jump.trees.final2.3RWGD.tau.fish$spe.jump.prob,info.jump.trees.final2.3RWGD.tau.fish$dup.jump.prob) #p-value = 2.36e-2
pval.jump.3RWGD.tau.fish.spe2WGD <- paired.wilcox(info.jump.trees.final2.3RWGD.tau.fish$spe.jump.prob,info.jump.trees.final2.3RWGD.tau.fish$WGD.jump.prob) #p-value = 1.81e-4

## For median data
dtfr.3RWGD.tau.fish<-NULL
dtfr.3RWGD.tau.fish<-data.frame(prop.jump=info.jump.trees.final2.3RWGD.tau.fish$spe.jump.prob,events="speciation")
dtfr.3RWGD.tau.fish<-rbind(dtfr.3RWGD.tau.fish,data.frame(prop.jump=info.jump.trees.final2.3RWGD.tau.fish$dup.jump.prob,events="duplication"))
dtfr.3RWGD.tau.fish<-rbind(dtfr.3RWGD.tau.fish,data.frame(prop.jump=info.jump.trees.final2.3RWGD.tau.fish$WGD.jump.prob,events="FishWGD"))
dtfr.3RWGD.tau.fish$Exp.abs <-abs(dtfr.3RWGD.tau.fish$prop.jump)
dtfr.3RWGD.tau.fish$Event<-factor(dtfr.3RWGD.tau.fish$events, levels=c("speciation","duplication","FishWGD"))

library(tidyverse)
median.data.jump.3RWGD.tau.fish<-NULL
median.data.jump.3RWGD.tau.fish<- dtfr.3RWGD.tau.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Tau.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.3RWGD.tau.fish$pic.round<-paste0("Median = ",round(median.data.jump.3RWGD.tau.fish$Tau.abs,4), sep="")
median.data.jump.3RWGD.tau.fish$count<- paste0(median.data.jump.3RWGD.tau.fish$Event,"\n","(n = ",median.data.jump.3RWGD.tau.fish$freq,")")

plot3R2<-boxplot.new.3RWGD2(dtfr.3RWGD.tau.fish,pval.jump.3RWGD.tau.fish.spe2duppval.jump.3RWGD.tau.fish.spe2WGD,"Tau.3RWGD",median.data.jump.3RWGD.tau.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("Jump in ", tau, ", for 5 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.05, xmax = 1.7, ymin = 0.67, ymax =0.67, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.4, y = 0.75, label= pval.jump.3RWGD.tau.fish.spe2dup, fontface = 4)+
  annotate("rect", xmin = 1.05, xmax = 2.8, ymin = 0.88, ymax =0.88, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.7, y = 0.95, label= pval.jump.3RWGD.tau.fish.spe2WGD, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.3RWGD.tau.fish$count) ## check


########### Testing OC for all fish specific nodes of 435 WGD trees supporting jump##################
nodes.contrast.3RWGD.tau.jumps.fish<-nodes.contrast.3RWGD.tau.jumps.15spe[which(nodes.contrast.3RWGD.tau.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 435 WGD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.3RWGD.tau.jumps.fish.removing.jump.node<-nodes.contrast.3RWGD.tau.jumps.15spe.removing.jump.node[which(nodes.contrast.3RWGD.tau.jumps.15spe.removing.jump.node$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

########### Testing OC for all fish specific nodes of 311 3RWGD trees that do not support trait jump for teleosts ##################
nodes.contrast.3RWGD.tau.nojumps.fish<-nodes.contrast.3RWGD.tau.nojumps.15spe[which(nodes.contrast.3RWGD.tau.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")

################## Chunk 5: Analysis on contrasts standardized sure SSD trees for traits ##########################

### Step 1: Considering (Calibrated.SSD.tau) 2482 SSD tree data passing the diagnostic test for Tau
### and matching their index with the output to identify standardized trees with and without trait jump

index.info1<-NULL
count<-0
index.info1<-bind_rows(lapply(Calibrated.SSD.tau,tree.index.collect))
rm(count)

## We need to match this indices with the previous one
colnames(index.info1)[2]<-colnames(index.map)[1]
index.info1<-merge(index.info1,index.map,by=c("index.6923tree.tau"))

test.output.tau.all<-NULL
test.output.tau.all<-merge(index.info1,test.output1.tau,by="Gene") ## 8861 observations
index.trees.jump<-unique(test.output.tau.all$tree.num) ##1557 trees
Calibrated.SSD.tau.jump<-Calibrated.SSD.tau[c(index.trees.jump)] ##1557/2482=62.73% Brownian trees supported trait jump model for Tau
Calibrated.SSD.tau.nojump<-Calibrated.SSD.tau[-c(index.trees.jump)] ## 925 trees with no support for jump in tau
rm(index.info1)
rm(index.trees.jump)

### Step2: Calculating PICs of required trait for trees to match the index number of the trees on which levolution was run
Count<-0
trees.all.interest<-lapply(trees.all.interest,contrast.calc,trait="Tau")

Count<-0
Calibrated.SSD.tau.jump<-lapply(Calibrated.SSD.tau.jump,contrast.calc,trait="Tau") 
rm(Count)

### Step3: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.tau.all))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output(test.output.tau.all)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 8861 obs
}

test.output.SSD.tau.processed<-bind_cols(test.output.tau.all,processed.fout) ## 8861 obs for all the 15 species
test.output.SSD.tau.processed<-test.output.SSD.tau.processed[(which(!(test.output.SSD.tau.processed$clade.name %in% c("Clupeocephala","Osteoglossocephalai")))),] ##8429 obs
rm(i)
rm(test.output.tau.all)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

## Is OC supported for vertebrate clades of SSD tau trees that are supporting jump?
pic.spe.SSD.tau<-test.output.SSD.tau.processed$pic[which(test.output.SSD.tau.processed$node.event=="speciation")]
pic.dup.SSD.tau<-test.output.SSD.tau.processed$pic[which(test.output.SSD.tau.processed$node.event=="duplication")]
median(pic.spe.SSD.tau) ##0.0098
median(pic.dup.SSD.tau) ##0.0131
Pval.SSD.tau.jump.all<-two.tailed.wilcox(pic.spe.SSD.tau,pic.dup.SSD.tau) # < 2.2e-16

## Step 4: Analysis on all 15 vertebrates clades
test.output.SSD.tau.processed.tau.new.mod<-unique(test.output.SSD.tau.processed[c(8,5,9,10,12,4,3)])  #6294 obs

########### Proportion of jump events in 15 vertebrates species in the 1557 SSD trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.tau.SSD.trees<-NULL
count<-0
info.jump.tau.SSD.trees<-bind_rows(lapply(Calibrated.SSD.tau.jump, tree.data.collection))
info.jump.tau.SSD.trees$spe.num<-info.jump.tau.SSD.trees$internal.events-info.jump.tau.SSD.trees$dup.num
info.jump.tau.SSD.trees<-info.jump.tau.SSD.trees[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.SSD.tau.all<-
  test.output.SSD.tau.processed %>%
  filter(node.event=="duplication" ) %>%
  .$index %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.SSD.tau.all)<-c("index.old","jump.dup")

test.spe.SSD.tau.all<-
  test.output.SSD.tau.processed%>%
  filter(node.event=="speciation") %>%
  .$index %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.SSD.tau.all)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.SSD.tau.all<-full_join(test.spe.SSD.tau.all,test.dup.SSD.tau.all)
info.jump.trees.final.SSD.tau.all<-merge(info.jump.tau.SSD.trees,merged.SSD.tau.all, by = c("index.old"))
info.jump.trees.final.SSD.tau.all[is.na(info.jump.trees.final.SSD.tau.all)] <- 0

## Have to consider trees where atleast one duplication and speciation exists
info.jump.trees.final1.SSD.tau.all<-info.jump.trees.final.SSD.tau.all[which(info.jump.trees.final.SSD.tau.all$dup.num>0 & info.jump.trees.final.SSD.tau.all$spe.num>0),]
rm(info.jump.trees.final.SSD.tau.all)

## Proportion of jumps
info.jump.trees.final1.SSD.tau.all$dup.jump.prob<-(info.jump.trees.final1.SSD.tau.all$jump.dup/(info.jump.trees.final1.SSD.tau.all$dup.num*2))
info.jump.trees.final1.SSD.tau.all$spe.jump.prob<-(info.jump.trees.final1.SSD.tau.all$jump.spe/(info.jump.trees.final1.SSD.tau.all$spe.num*2))

## Stats
median(na.omit(info.jump.trees.final1.SSD.tau.all$dup.jump.prob)) #0.25
median(na.omit(info.jump.trees.final1.SSD.tau.all$spe.jump.prob)) #0.05
mean(na.omit(info.jump.trees.final1.SSD.tau.all$dup.jump.prob)) #0.30
mean(na.omit(info.jump.trees.final1.SSD.tau.all$spe.jump.prob)) #0.13

pval.jump.SSD.tau.all <- paired.wilcox(info.jump.trees.final1.SSD.tau.all$spe.jump.prob,info.jump.trees.final1.SSD.tau.all$dup.jump.prob) #p-value = 2.2e-16

## For median data
dtfr.SSD.tau.15spe<-NULL
dtfr.SSD.tau.15spe<-data.frame(prop.jump=info.jump.trees.final1.SSD.tau.all$spe.jump.prob,events="speciation")
dtfr.SSD.tau.15spe<-rbind(dtfr.SSD.tau.15spe,data.frame(prop.jump=info.jump.trees.final1.SSD.tau.all$dup.jump.prob,events="duplication"))
dtfr.SSD.tau.15spe$tau.abs <-abs(dtfr.SSD.tau.15spe$prop.jump)
dtfr.SSD.tau.15spe$Event<-factor(dtfr.SSD.tau.15spe$events, levels=c("speciation","duplication"))

library(tidyverse)
median.data.jump.SSD.tau<-NULL
median.data.jump.SSD.tau<- dtfr.SSD.tau.15spe%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            tau.abs=median(prop.jump),
            freq = n())
median.data.jump.SSD.tau$pic.round<-paste0("Median = ",round(median.data.jump.SSD.tau$Tau.abs,4), sep="")
median.data.jump.SSD.tau$count<- paste0(median.data.jump.SSD.tau$Event,"\n\n","(n = ",median.data.jump.SSD.tau$freq,")")

plot1<-boxplot.new2(dtfr.SSD.tau.15spe,pval.jump.SSD.tau.alll,"Tau",median.data.jump.SSD.tau) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.SSD.tau.all, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.SSD.tau$count) ## check

# Create Data
data.jump.SSD.tau <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(1557,925)
)

# Compute the position of labels
data.jump.SSD.tau <- data.jump.SSD.tau %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.SSD.tau$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.SSD.tau$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.SSD.tau<-ggplot(data.jump.SSD.tau, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 1557 SSD trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized tau levels 
summary.SSD.tau.jump<-summary.tree(Calibrated.SSD.tau.jump)
summary.SSD.tau.jump$Tau.abs <- abs(summary.SSD.tau.jump$pic_Tau)
summary.SSD.tau.jumps<-summary.SSD.tau.jump[,!(names(summary.SSD.tau.jump) %in% c("pic_Tau"))]
nodes.contrast.SSD.tau.jumps.15spe <- summary.SSD.tau.jumps[which(!is.na(summary.SSD.tau.jumps$Tau.abs)),]
nodes.contrast.SSD.tau.jumps.15spe$Event <- factor(nodes.contrast.SSD.tau.jumps.15spe$events, levels=c("speciation", "duplication"))
nodes.contrast.SSD.tau.jumps.15spe<-nodes.contrast.SSD.tau.jumps.15spe[(which(!(nodes.contrast.SSD.tau.jumps.15spe$label %in% c("Clupeocephala","Osteoglossocephalai")))),]
rm(summary.SSD.tau.jump)
rm(summary.SSD.tau.jumps)

############ Testing OC using all 15 species for 1557 SSD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
library(dplyr)
nodes.contrast.SSD.tau.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.SSD.tau.jumps.15spe$index<-nodes.contrast.SSD.tau.jumps.15spe$index.tree
test.output.SSD.tau.processed$node<-test.output.SSD.tau.processed$From
nodes.contrast.SSD.tau.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.SSD.tau.jumps.15spe,test.output.SSD.tau.processed,by=c("index","node"))
nodes.contrast.SSD.tau.jumps.15spe.removing.jump.node <- nodes.contrast.SSD.tau.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.SSD.tau.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.SSD.tau.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.SSD.tau.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.SSD.tau.jumps.15spe.removing.jump.node$Tau.abs<-nodes.contrast.SSD.tau.jumps.15spe.removing.jump.node$pic
nodes.contrast.SSD.tau.jumps.15spe.removing.jump.node<-nodes.contrast.SSD.tau.jumps.15spe.removing.jump.node[(which(!(nodes.contrast.SSD.tau.jumps.15spe.removing.jump.node$label %in% c("Clupeocephala","Osteoglossocephalai")))),]

pic.spe.removing.jump.15species.SSD.tau.trees<-speciation.contrast(nodes.contrast.SSD.tau.jumps.15spe.removing.jump.node,"pic") #(n=17064)
pic.dup.removing.jump.15species.SSD.tau.trees<-duplication.contrast(nodes.contrast.SSD.tau.jumps.15spe.removing.jump.node,"pic") #(n=2257)


########### Testing OC for all nodes of 925 SSD tau trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized tau levels 
summary.SSD.tau.nojump<-summary.tree(Calibrated.SSD.tau.nojump)
summary.SSD.tau.nojump$Tau.abs <- abs(summary.SSD.tau.nojump$pic_Tau)
summary.SSD.tau.nojumps<-summary.SSD.tau.nojump[,!(names(summary.SSD.tau.nojump) %in% c("pic_Tau"))]
nodes.contrast.SSD.tau.nojumps.15spe <- summary.SSD.tau.nojumps[which(!is.na(summary.SSD.tau.nojumps$Tau.abs)),]
nodes.contrast.SSD.tau.nojumps.15spe$Event <- factor(nodes.contrast.SSD.tau.nojumps.15spe$events, levels=c("speciation", "duplication"))
nodes.contrast.SSD.tau.nojumps.15spe<-nodes.contrast.SSD.tau.nojumps.15spe[(which(!(nodes.contrast.SSD.tau.nojumps.15spe$label %in% c("Clupeocephala","Osteoglossocephalai")))),] #13239 obs
rm(summary.SSD.tau.nojump)
rm(summary.SSD.tau.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.SSD.tau.new.fish<-test.output.SSD.tau.processed[which(test.output.SSD.tau.processed$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##2640 obs

## Step 6: Analysis on fish specific clades
test.output.SSD.tau.new.fish.mod<-unique(test.output.SSD.tau.new.fish[c(8,5,9,10,12,4,3)])  #1880 obs
rm(test.output.SSD.tau.new.fish)

########### Proportion of jump events in teleosts in the 1557 SSD trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.SSD.tau.trees.fish<-NULL
count<-0
info.jump.SSD.tau.trees.fish<-bind_rows(lapply(Calibrated.SSD.tau.jump, tree.data.collection.fish))
info.jump.SSD.tau.trees.fish$spe.num<-info.jump.SSD.tau.trees.fish$internal.events.fish-info.jump.SSD.tau.trees.fish$dup.num
info.jump.SSD.tau.trees.fish<-info.jump.SSD.tau.trees.fish[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.SSD.tau.fish<-
  test.output.SSD.tau.new.fish  %>%
  filter(node.event=="duplication") %>%
  .$index %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.SSD.tau.fish)<-c("index.old","jump.dup")

test.spe.SSD.tau.fish<-
  test.output.SSD.tau.new.fish %>%
  filter(node.event=="speciation") %>%
  .$index %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.SSD.tau.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.SSD.tau.fish<-full_join(test.spe.SSD.tau.fish,test.dup.SSD.tau.fish)
info.jump.trees.SSD.tau.final.fish<-merge(info.jump.SSD.tau.trees.fish,merged.SSD.tau.fish, by = c("index.old"))
info.jump.trees.SSD.tau.final.fish[is.na(info.jump.trees.SSD.tau.final.fish)] <- 0


## Have to consider trees where atleast one fish specific duplication and speciation exists
info.jump.trees.SSD.tau.final1.fish<-info.jump.trees.SSD.tau.final.fish[which(info.jump.trees.SSD.tau.final.fish$dup.num>0 & info.jump.trees.SSD.tau.final.fish$spe.num>0),]

## Proportion of jumps
info.jump.trees.SSD.tau.final1.fish$dup.jump.prob<-(info.jump.trees.SSD.tau.final1.fish$jump.dup/(info.jump.trees.SSD.tau.final1.fish$dup.num*2))
info.jump.trees.SSD.tau.final1.fish$spe.jump.prob<-(info.jump.trees.SSD.tau.final1.fish$jump.spe/(info.jump.trees.SSD.tau.final1.fish$spe.num*2))


median(na.omit(info.jump.trees.SSD.tau.final1.fish$dup.jump.prob)) #0.25
median(na.omit(info.jump.trees.SSD.tau.final1.fish$spe.jump.prob)) #0.125
mean(na.omit(info.jump.trees.SSD.tau.final1.fish$dup.jump.prob)) #0.34
mean(na.omit(info.jump.trees.SSD.tau.final1.fish$spe.jump.prob)) #0.18

pval.jump.SSD.tau.fish <- paired.wilcox(info.jump.trees.SSD.tau.final1.fish$spe.jump.prob,info.jump.trees.SSD.tau.final1.fish$dup.jump.prob) #p-value = 2.2e-16

## For median data
dtfr.SSD.tau.fish<-NULL
dtfr.SSD.tau.fish<-data.frame(prop.jump=info.jump.trees.SSD.tau.final1.fish$spe.jump.prob,events="speciation")
dtfr.SSD.tau.fish<-rbind(dtfr.SSD.tau.fish,data.frame(prop.jump=info.jump.trees.SSD.tau.final1.fish$dup.jump.prob,events="duplication"))
dtfr.SSD.tau.fish$Tau.abs <-abs(dtfr.SSD.tau.fish$prop.jump)
dtfr.SSD.tau.fish$Event<-factor(dtfr.SSD.tau.fish$events, levels=c("speciation","duplication"))


library(tidyverse)
median.data.jump.SSD.tau.fish<-NULL
median.data.jump.SSD.tau.fish<- dtfr.SSD.tau.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            tau.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.SSD.tau.fish$pic.round<-paste0("Median = ",round(median.data.jump.SSD.tau.fish$Tau.abs,4), sep="")
median.data.jump.SSD.tau.fish$count<- paste0(median.data.jump.SSD.tau.fish$Event,"\n\n","(n = ",median.data.jump.SSD.tau.fish$freq,")")

plot2<-boxplot.new2(dtfr.SSD.tau.fish,pval.jump.SSD.tau.fish,"tau",median.data.jump.SSD.tau.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 5 teleosts"))))+
  #ylab(expression(bold(paste("Proportions of jump events for teleosts ")))) +
  ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.SSD.tau.fish, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.SSD.tau.fish$count) ## check

########### Testing OC for all fish specific nodes of 1557 SSD trees supporting jump##################
nodes.contrast.SSD.tau.jumps.fish<-nodes.contrast.SSD.tau.jumps.15spe[which(nodes.contrast.SSD.tau.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 1557 SSD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
library(dplyr)
nodes.contrast.SSD.tau.jumps.fish.removing.jump.node <- NULL
nodes.contrast.SSD.tau.jumps.fish$index<-nodes.contrast.SSD.tau.jumps.fish$index.tree
test.output.SSD.tau.new.fish$node<-test.output.SSD.tau.new.fish$From
nodes.contrast.SSD.tau.jumps.fish.removing.jump.node <-anti_join(nodes.contrast.SSD.tau.jumps.fish,test.output.SSD.tau.new.fish,by=c("index","node"))
nodes.contrast.SSD.tau.jumps.fish.removing.jump.node <- nodes.contrast.SSD.tau.jumps.fish.removing.jump.node[which(!is.na(nodes.contrast.SSD.tau.jumps.fish.removing.jump.node$pic)),]
nodes.contrast.SSD.tau.jumps.fish.removing.jump.node$Event <- factor(nodes.contrast.SSD.tau.jumps.fish.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.SSD.tau.jumps.fish.removing.jump.node$Tau.abs<-nodes.contrast.SSD.tau.jumps.fish.removing.jump.node$pic
nodes.contrast.SSD.tau.jumps.fish.removing.jump.node<-nodes.contrast.SSD.tau.jumps.fish.removing.jump.node[(which(!(nodes.contrast.SSD.tau.jumps.fish.removing.jump.node$label %in% c("Clupeocephala","Osteoglossocephalai")))),]


########### Testing OC for all fish specific nodes of 925 SSD tau trees that do not support trait jump for teleosts ##################
nodes.contrast.SSD.tau.nojumps.fish<-nodes.contrast.SSD.tau.nojumps.15spe[which(nodes.contrast.SSD.tau.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] #3453 obs

#save.image("norm_exp_15spe_4tips_new.RData")



