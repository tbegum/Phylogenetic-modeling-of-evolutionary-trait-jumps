

### This code is written to work with "trait jump model" for muscle

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
load("norm_exp_15spe_4tips_new.RData")  ## latest modified exp 

######################### Chunk1: Data given to pablo for applying "Levolution" to find trait jump #############################
## trees.all.interest with 6923 gene trees given to pablo

######################### Chunk2: reading output of trait jump #############################

test.output.muscle<-read.csv(paste0(folder1,"To_pablo/muscle/results_muscle.csv", sep=""),sep=",",header=T) ##207192 obs

## Only considering trees with posterior probability of jump more than 70%
test.output.muscle1<-test.output.muscle[which(test.output.muscle$JumpProbability>=0.7),] ##24825 obs

########################## Chunk 3: Analysis on all contrasts standardized trees for traits ###############################
### Step 1a: Considering (trees.all.interest) 6923 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for muscle tissue for all trees 
count<-0
Calibrated.muscle.exp1<-lapply(trees.all.interest, diagnostic.plot.test.exp.all.tissue, "TPM.muscle")
Calibrated.muscle.exp2<-Calibrated.muscle.exp1[!is.na(Calibrated.muscle.exp1)] 
Calibrated.all.standardized.muscle<-Calibrated.muscle.exp2[ ! sapply(Calibrated.muscle.exp2, is.null) ]## finally we obtained 4600/6923 tree data passing the diagnostic tests
rm(Calibrated.muscle.exp1)
rm(Calibrated.muscle.exp2)

### Step 1b: Considering (Calibrated.all.standardized.muscle) 4600 tree data passing the diagnostic test for muscle expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.all.standardized.muscle,tree.index.collect))
rm(count)

test.output.muscle.all<-NULL
test.output.muscle.all<-merge(index.info,test.output.muscle1,by="Gene") ## 15378 observations
index.trees.jump<-unique(test.output.muscle.all$tree.num) ##3800 trees
Calibrated.standardized.muscle.jump<-Calibrated.all.standardized.muscle[c(index.trees.jump)] ##3800/4600=82.61% Brownian trees supported trait jump model for muscle expressions 
Calibrated.standardized.muscle.nojump<-Calibrated.all.standardized.muscle[-c(index.trees.jump)] ## 800/4600=17.39% trees with no support for jump in muscle expression
rm(index.info)
rm(index.trees.jump)

### Step2: Calculating PICs of required trait for trees to match the index number of the trees on which levolution was run
Count<-0
trees.all.interest<-lapply(trees.all.interest,contrast.calc,trait="TPM.muscle")
# Count<-0
# Calibrated.standardized.muscle.jump<-lapply(Calibrated.standardized.muscle.jump,contrast.calc,trait="TPM.muscle") 
# rm(Count)

### Step3: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.muscle.all))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.muscle.all)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 15378 obs
}
test.output.muscle.processed<-bind_cols(test.output.muscle.all,processed.fout) ## 15378 obs for all the 15 species
rm(i)
rm(test.output.muscle.all)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

## Step 4: Analysis on all 15 vertebrates clades
test.output.muscle.new.mod<-unique(test.output.muscle.processed[c(6,3,7,8,10,2,1)])  #12727 obs

########### Proportion of jump events in 15 vertebrates species in the 3495 trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.trees.muscle<-NULL
count<-0
info.jump.trees.muscle<-bind_rows(lapply(Calibrated.standardized.muscle.jump, tree.data.collection))
info.jump.trees.muscle$spe.num<-info.jump.trees.muscle$internal.events-info.jump.trees.muscle$dup.num
info.jump.trees.muscle<-info.jump.trees.muscle[c(11,12,3:5)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.muscle<-
  test.output.muscle.processed  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.muscle)<-c("index.old","jump.dup")

test.spe.muscle<-
  test.output.muscle.processed %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.muscle)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.muscle<-full_join(test.spe.muscle,test.dup.muscle)
info.jump.trees.muscle.final<-merge(info.jump.trees.muscle,merged.muscle, by = c("index.old"))
info.jump.trees.muscle.final[is.na(info.jump.trees.muscle.final)] <- 0

## Have to consider trees where atleast one duplication and speciation exists
info.jump.trees.muscle.final1<-info.jump.trees.muscle.final[which(info.jump.trees.muscle.final$dup.num>0 & info.jump.trees.muscle.final$spe.num>0),]
rm(info.jump.trees.muscle.final)

## Proportion of jumps
info.jump.trees.muscle.final1$dup.jump.prob<-(info.jump.trees.muscle.final1$jump.dup/(info.jump.trees.muscle.final1$dup.num*2))
info.jump.trees.muscle.final1$spe.jump.prob<-(info.jump.trees.muscle.final1$jump.spe/(info.jump.trees.muscle.final1$spe.num*2))

## Stats
median(na.omit(info.jump.trees.muscle.final1$dup.jump.prob)) #0
median(na.omit(info.jump.trees.muscle.final1$spe.jump.prob)) #0.045
mean(na.omit(info.jump.trees.muscle.final1$dup.jump.prob)) #0.1618
mean(na.omit(info.jump.trees.muscle.final1$spe.jump.prob)) #0.1027

pval.jump.muscle <- paired.wilcox(info.jump.trees.muscle.final1$spe.jump.prob,info.jump.trees.muscle.final1$dup.jump.prob) #p-value = 2.05e-1

## For median data
dtfr.muscle<-NULL
dtfr.muscle<-data.frame(prop.jump=info.jump.trees.muscle.final1$spe.jump.prob,events="speciation")
dtfr.muscle<-rbind(dtfr.muscle,data.frame(prop.jump=info.jump.trees.muscle.final1$dup.jump.prob,events="duplication"))
dtfr.muscle$Exp.abs <-abs(dtfr.muscle$prop.jump)
dtfr.muscle$Event<-factor(dtfr.muscle$events, levels=c("speciation","duplication"))

library(tidyverse)
median.data.jump.muscle<-NULL
median.data.jump.muscle<- dtfr.muscle%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.muscle$pic.round<-paste0("Median = ",round(median.data.jump.muscle$Exp.abs,4), sep="")
median.data.jump.muscle$count<- paste0(median.data.jump.muscle$Event,"\n\n","(n = ",median.data.jump.muscle$freq,")")

plot4A<-boxplot.new2(dtfr.muscle,pval.jump.muscle,"Exp",median.data.jump.muscle) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.muscle, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.muscle$count) ## check

# Create Data
data.jump.muscle <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(3800,800)
)

# Compute the position of labels
data.jump.muscle <- data.jump.muscle %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.muscle$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.muscle$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.muscle<-ggplot(data.jump.muscle, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 3495 trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.muscle.jump<-summary.tree(Calibrated.standardized.muscle.jump)
summary.muscle.jump$Exp.abs <- abs(summary.muscle.jump$pic.muscle)
summary.muscle.jumps<-summary.muscle.jump[,!(names(summary.muscle.jump) %in% c("pic.muscle"))]
nodes.contrast.muscle.jumps.15spe <- summary.muscle.jumps[which(!is.na(summary.muscle.jumps$Exp.abs)),]
nodes.contrast.muscle.jumps.15spe$Event <- factor(nodes.contrast.muscle.jumps.15spe$events, levels=c("speciation", "duplication"))
rm(summary.muscle.jump)
rm(summary.muscle.jumps)

############ Testing OC using all 15 species for 3495 trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.muscle.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.muscle.jumps.15spe$index<-nodes.contrast.muscle.jumps.15spe$index.tree
test.output.muscle.processed$node<-test.output.muscle.processed$From
nodes.contrast.muscle.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.muscle.jumps.15spe,test.output.muscle.processed,by=c("index","node"))
nodes.contrast.muscle.jumps.15spe.removing.jump.node <- nodes.contrast.muscle.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.muscle.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.muscle.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.muscle.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.muscle.jumps.15spe.removing.jump.node$Exp.abs<-nodes.contrast.muscle.jumps.15spe.removing.jump.node$pic

########### Testing OC for all nodes of 1427 trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.muscle.nojump<-summary.tree(Calibrated.standardized.muscle.nojump)
summary.muscle.nojump$Exp.abs <- abs(summary.muscle.nojump$pic.muscle)
summary.muscle.nojumps<-summary.muscle.nojump[,!(names(summary.muscle.nojump) %in% c("pic.muscle"))]
nodes.contrast.muscle.nojumps.15spe <- summary.muscle.nojump[which(!is.na(summary.muscle.nojump$Exp.abs)),]
nodes.contrast.muscle.nojumps.15spe$Event <- factor(nodes.contrast.muscle.nojumps.15spe$events, levels=c("speciation", "duplication"))
rm(summary.muscle.nojump)
rm(summary.muscle.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.muscle.new.fish<-test.output.muscle.processed[which(test.output.muscle.processed$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##5908 obs

## Step 6: Analysis on fish specific clades
test.output.muscle.new.fish.mod<-unique(test.output.muscle.new.fish[c(6,3,7,8,10,2,1)])  #4830 obs

########### Proportion of jump events in teleosts in the 3495 trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.trees.muscle.fish<-NULL
count<-0
info.jump.trees.muscle.fish<-bind_rows(lapply(Calibrated.standardized.muscle.jump, tree.data.collection.fish))
info.jump.trees.muscle.fish$spe.num<-info.jump.trees.muscle.fish$internal.events.fish-info.jump.trees.muscle.fish$dup.num
info.jump.trees.muscle.fish<-info.jump.trees.muscle.fish[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.muscle.fish<-
  test.output.muscle.new.fish  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.muscle.fish)<-c("index.old","jump.dup")

test.spe.muscle.fish<-
  test.output.muscle.new.fish %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.muscle.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.muscle.fish<-full_join(test.spe.muscle.fish,test.dup.muscle.fish)
info.jump.trees.muscle.final.fish<-merge(info.jump.trees.muscle.fish,merged.muscle.fish, by = c("index.old"))
info.jump.trees.muscle.final.fish[is.na(info.jump.trees.muscle.final.fish)] <- 0

## Have to consider trees where atleast one fish specific duplication and speciation exists
info.jump.trees.muscle.final1.fish<-info.jump.trees.muscle.final.fish[which(info.jump.trees.muscle.final.fish$dup.num>0 & info.jump.trees.muscle.final.fish$spe.num>0),]

## Proportion of jumps
info.jump.trees.muscle.final1.fish$dup.jump.prob<-(info.jump.trees.muscle.final1.fish$jump.dup/(info.jump.trees.muscle.final1.fish$dup.num*2))
info.jump.trees.muscle.final1.fish$spe.jump.prob<-(info.jump.trees.muscle.final1.fish$jump.spe/(info.jump.trees.muscle.final1.fish$spe.num*2))


median(na.omit(info.jump.trees.muscle.final1.fish$dup.jump.prob)) #0.0
median(na.omit(info.jump.trees.muscle.final1.fish$spe.jump.prob)) #0.125
mean(na.omit(info.jump.trees.muscle.final1.fish$dup.jump.prob)) #0.169
mean(na.omit(info.jump.trees.muscle.final1.fish$spe.jump.prob)) #0.1869

pval.jump.muscle.fish <- paired.wilcox(info.jump.trees.muscle.final1.fish$spe.jump.prob,info.jump.trees.muscle.final1.fish$dup.jump.prob) #p-value = 2.1e-7


## For median data
dtfr.muscle.fish<-NULL
dtfr.muscle.fish<-data.frame(prop.jump=info.jump.trees.muscle.final1.fish$spe.jump.prob,events="speciation")
dtfr.muscle.fish<-rbind(dtfr.muscle.fish,data.frame(prop.jump=info.jump.trees.muscle.final1.fish$dup.jump.prob,events="duplication"))
dtfr.muscle.fish$Exp.abs <-abs(dtfr.muscle.fish$prop.jump)
dtfr.muscle.fish$Event<-factor(dtfr.muscle.fish$events, levels=c("speciation","duplication"))


library(tidyverse)
median.data.jump.muscle.muscle.fish<-NULL
median.data.jump.muscle.muscle.fish<- dtfr.muscle.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.muscle.muscle.fish$pic.round<-paste0("Median = ",round(median.data.jump.muscle.muscle.fish$Exp.abs,4), sep="")
median.data.jump.muscle.muscle.fish$count<- paste0(median.data.jump.muscle.muscle.fish$Event,"\n\n","(n = ",median.data.jump.muscle.muscle.fish$freq,")")

plot4B<-boxplot.new2(dtfr.muscle.fish,pval.jump.muscle.fish,"Exp",median.data.jump.muscle.muscle.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 5 teleosts"))))+
  #ylab(expression(bold(paste("Proportions of jump events for teleosts ")))) +
  ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.muscle.fish, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.muscle.muscle.fish$count) ## check

########### Testing OC for all fish specific nodes of 3495 trees supporting jump##################
nodes.contrast.muscle.jumps.fish<-nodes.contrast.muscle.jumps.15spe[which(nodes.contrast.muscle.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 3495 trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.muscle.jumps.fish.removing.jump.node <- NULL
nodes.contrast.muscle.jumps.fish$index<-nodes.contrast.muscle.jumps.fish$index.tree
test.output.muscle.new.fish$node<-test.output.muscle.new.fish$From
nodes.contrast.muscle.jumps.fish.removing.jump.node <-anti_join(nodes.contrast.muscle.jumps.fish,test.output.muscle.new.fish,by=c("index","node"))
nodes.contrast.muscle.jumps.fish.removing.jump.node <- nodes.contrast.muscle.jumps.fish.removing.jump.node[which(!is.na(nodes.contrast.muscle.jumps.fish.removing.jump.node$pic)),]
nodes.contrast.muscle.jumps.fish.removing.jump.node$Event <- factor(nodes.contrast.muscle.jumps.fish.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.muscle.jumps.fish.removing.jump.node$Exp.abs<-nodes.contrast.muscle.jumps.fish.removing.jump.node$pic

rm(test.output.muscle.new.fish)

########### Testing OC for all fish specific nodes of 1427 trees that do not support trait jump for teleosts ##################
nodes.contrast.muscle.nojumps.fish<-nodes.contrast.muscle.nojumps.15spe[which(nodes.contrast.muscle.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")

################## Chunk 4: Analysis on contrasts standardized sure 3RWGD trees for traits ##########################

### Step 1a: Considering (Sure.WGD.trees) 1159 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for muscle tissue for sure WGD trees 
count<-0
Calibrated.muscle.3Rexp1<-lapply(Sure.WGD.trees, diagnostic.plot.test.exp.all.tissue, "TPM.muscle")
Calibrated.muscle.3Rexp2<-Calibrated.muscle.3Rexp1[!is.na(Calibrated.muscle.3Rexp1)] 
Calibrated.3R.standardized.muscle<-Calibrated.muscle.3Rexp2[ ! sapply(Calibrated.muscle.3Rexp2, is.null) ]## finally we obtained 725/1159 trees passing the diagnostic tests
rm(Calibrated.muscle.3Rexp1)
rm(Calibrated.muscle.3Rexp2)

### Step 1b: Considering (Calibrated.3R.standardized.muscle) 725 trees passing the diagnostic test for muscle expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.3R.standardized.muscle,tree.index.collect))
rm(count)

test.output.3Rexp.muscle<-NULL
test.output.3Rexp.muscle<-merge(index.info,test.output.muscle1,by="Gene") ## 2125 observations
index.trees.jump<-unique(test.output.3Rexp.muscle$tree.num) ##586 trees
Calibrated.3RWGD.muscle.jump<-Calibrated.3R.standardized.muscle[c(index.trees.jump)] ##586/725=80.83% Brownian trees supported trait jump model for muscle expressions
Calibrated.3RWGD.muscle.nojump<-Calibrated.3R.standardized.muscle[-c(index.trees.jump)] ## 139/725=31.11% trees with no support for jump in trait
rm(index.info)
rm(index.trees.jump)

### Step2: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.3Rexp.muscle))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.3Rexp.muscle)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 2125 obs
}

test.output.3RWGD.muscle.processed<-bind_cols(test.output.3Rexp.muscle,processed.fout) ## 2125 obs for all the 15 species
rm(i)
rm(test.output.3Rexp.muscle)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

##Since our calibrated tree shows both the FishWGD and duplication nodes, 
##we want to modify the node events of the jump output to separate out the effect of FishWGD from the SSDs in those trees.

jump.3RWGD.muscle.trees.summary<-summary.tree(Calibrated.3RWGD.muscle.jump)
jump.3RWGD.muscle.trees.summary$Gene<-jump.3RWGD.muscle.trees.summary$index.tree
jump.3RWGD.muscle.trees.summary$From<-jump.3RWGD.muscle.trees.summary$node
test.output.3RWGD.muscle.processed.mod<-merge(test.output.3RWGD.muscle.processed,jump.3RWGD.muscle.trees.summary,by=c("Gene","From"))
test.output.3RWGD.muscle.processed.mod<-test.output.3RWGD.muscle.processed.mod[c(1,3,2,4:11,16)]
test.output.3RWGD.muscle.processed.mod$node.event<-test.output.3RWGD.muscle.processed.mod$events
colnames(test.output.3RWGD.muscle.processed.mod)[10]<-"pic"
test.output.3RWGD.muscle.processed.mod<-test.output.3RWGD.muscle.processed.mod[c(1:11)] ## 2125 obs

## Step 3: Analysis on all 15 vertebrates clades
test.output.3RWGD.muscle.processed.new<-unique(test.output.3RWGD.muscle.processed.mod[c(6,3,7,8,10,2,1)])  #1850 obs

########### Proportion of jump events in 15 vertebrates species in the 586 3RWGD trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.3RWGD.muscle<-NULL
info.jump.3RWGD.muscle<-bind_rows(lapply(Calibrated.3RWGD.muscle.jump, tree.data.collection.3RWGD))
info.jump.3RWGD.muscle$spe.num<-info.jump.3RWGD.muscle$internal.events-(info.jump.3RWGD.muscle$dup.num + info.jump.3RWGD.muscle$WGD.num)
info.jump.3RWGD.muscle<-info.jump.3RWGD.muscle[c(11,14,3,12,5,4,13)]


## Counting frquency of jump for "speciation" and "duplication" events
test.dup.3RWGD.muscle<-
  test.output.3RWGD.muscle.processed.mod %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.3RWGD.muscle)<-c("index.old","jump.dup")

test.WGD.3RWGD.muscle<-
  test.output.3RWGD.muscle.processed.mod %>%
  filter(node.event=="FishWGD" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.WGD.3RWGD.muscle)<-c("index.old","jump.WGD")

test.spe.3RWGD.muscle<-
  test.output.3RWGD.muscle.processed.mod %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.3RWGD.muscle)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.3RWGD.muscle1<-full_join(test.spe.3RWGD.muscle,test.dup.3RWGD.muscle)
merged.3RWGD.muscle<-full_join(merged.3RWGD.muscle1,test.WGD.3RWGD.muscle)
info.jump.3RWGD.muscle.final1<-merge(info.jump.3RWGD.muscle,merged.3RWGD.muscle, by = c("index.old"))
info.jump.3RWGD.muscle.final1[is.na(info.jump.3RWGD.muscle.final1)] <- 0
rm(merged.3RWGD.muscle1)

## Have to consider trees where atleast one WGD and one speciation exists
info.jump.3RWGD.muscle.final<-info.jump.3RWGD.muscle.final1[which(info.jump.3RWGD.muscle.final1$spe.num>0 & info.jump.3RWGD.muscle.final1$WGD.num>0),]
rm(info.jump.3RWGD.muscle.final1)

## Proportion of jumps
info.jump.3RWGD.muscle.final$dup.jump.prob<-(info.jump.3RWGD.muscle.final$jump.dup/(info.jump.3RWGD.muscle.final$dup.num*2))
info.jump.3RWGD.muscle.final$spe.jump.prob<-(info.jump.3RWGD.muscle.final$jump.spe/(info.jump.3RWGD.muscle.final$spe.num*2))
info.jump.3RWGD.muscle.final$WGD.jump.prob<-(info.jump.3RWGD.muscle.final$jump.WGD/(info.jump.3RWGD.muscle.final$WGD.num*2))

## Since,many WGD trees didnot have a duplication node, it produces NA for dup.jump.prob
## To get rid of this
info.jump.3RWGD.muscle.final2<-info.jump.3RWGD.muscle.final[complete.cases(info.jump.3RWGD.muscle.final),] ##227 obs

## Stats
median(info.jump.3RWGD.muscle.final2$dup.jump.prob) #0
median(info.jump.3RWGD.muscle.final2$spe.jump.prob) #0.038
median(info.jump.3RWGD.muscle.final2$WGD.jump.prob) #0
mean(info.jump.3RWGD.muscle.final2$dup.jump.prob) #0.1495
mean(info.jump.3RWGD.muscle.final2$spe.jump.prob) #0.086
mean(info.jump.3RWGD.muscle.final2$WGD.jump.prob) #0.046

pval.jump.3RWGD.muscle.spe2dup <- paired.wilcox(info.jump.3RWGD.muscle.final2$spe.jump.prob,info.jump.3RWGD.muscle.final2$dup.jump.prob) #p-value = 1.23e-1
pval.jump.3RWGD.muscle.spe2WGD <- paired.wilcox(info.jump.3RWGD.muscle.final2$spe.jump.prob,info.jump.3RWGD.muscle.final2$WGD.jump.prob) #p-value = 2.2e-16

## For median data
dtfr.3RWGD.muscle.15spe<-NULL
dtfr.3RWGD.muscle.15spe<-data.frame(prop.jump=info.jump.3RWGD.muscle.final2$spe.jump.prob,events="speciation")
dtfr.3RWGD.muscle.15spe<-rbind(dtfr.3RWGD.muscle.15spe,data.frame(prop.jump=info.jump.3RWGD.muscle.final2$dup.jump.prob,events="duplication"))
dtfr.3RWGD.muscle.15spe<-rbind(dtfr.3RWGD.muscle.15spe,data.frame(prop.jump=info.jump.3RWGD.muscle.final2$WGD.jump.prob,events="FishWGD"))
dtfr.3RWGD.muscle.15spe$muscle.abs <-abs(dtfr.3RWGD.muscle.15spe$prop.jump)
dtfr.3RWGD.muscle.15spe$Event<-factor(dtfr.3RWGD.muscle.15spe$events, levels=c("speciation","duplication","FishWGD"))


median.data.jump.3RWGD.muscle.15spe<-NULL
median.data.jump.3RWGD.muscle.15spe<- dtfr.3RWGD.muscle.15spe%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.3RWGD.muscle.15spe$pic.round<-paste0("Median = ",round(median.data.jump.3RWGD.muscle.15spe$Exp.abs,4), sep="")
median.data.jump.3RWGD.muscle.15spe$count<- paste0(median.data.jump.3RWGD.muscle.15spe$Event,"\n\n","(n = ",median.data.jump.3RWGD.muscle.15spe$freq,")")

plot3R1<-boxplot.new.3RWGD2(dtfr.3RWGD.muscle.15spe,pval.jump.3RWGD.muscle.spe2dup,pval.jump.3RWGD.muscle.spe2WGD,"Exp.3RWGD",median.data.jump.3RWGD.muscle.15spe) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.05, xmax = 1.7, ymin = 0.67, ymax =0.67, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.4, y = 0.75, label= pval.jump.3RWGD.muscle.spe2dup, fontface = 4)+
  annotate("rect", xmin = 1.05, xmax = 2.8, ymin = 0.88, ymax =0.88, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.7, y = 0.95, label= pval.jump.3RWGD.muscle.spe2WGD, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.3RWGD.muscle.15spe$count) ## check

# Create Data
data.jump.3RWGD.exp <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(586,139)
)

# Compute the position of labels
data.jump.3RWGD.muscle <- data.jump.3RWGD.muscle %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.3RWGD.muscle$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.3RWGD.muscle$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.3RWGD.muscle<-ggplot(data.jump.3RWGD.muscle, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 586 WGD trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.3RWGD.muscle.jump<-summary.tree(Calibrated.3RWGD.muscle.jump)
summary.3RWGD.muscle.jump$Exp.abs <- abs(summary.3RWGD.muscle.jump$pic)
summary.3RWGD.muscle.jumps<-summary.3RWGD.muscle.jump[,!(names(summary.3RWGD.muscle.jump) %in% c("pic"))]
nodes.contrast.3RWGD.muscle.jumps.15spe <- summary.3RWGD.muscle.jumps[which(!is.na(summary.3RWGD.muscle.jumps$Exp.abs)),]
nodes.contrast.3RWGD.muscle.jumps.15spe$Event <- factor(nodes.contrast.3RWGD.muscle.jumps.15spe$events, levels=c("speciation", "duplication", "FishWGD"))
rm(summary.3RWGD.muscle.jump)
rm(summary.3RWGD.muscle.jumps)

############ Testing OC using all 15 species for 586 WGD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.3RWGD.muscle.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.3RWGD.muscle.jumps.15spe$index<-nodes.contrast.3RWGD.muscle.jumps.15spe$index.tree
test.output.3RWGD.muscle.processed.mod$node<-test.output.3RWGD.muscle.processed.mod$From
nodes.contrast.3RWGD.muscle.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.3RWGD.muscle.jumps.15spe,test.output.3RWGD.muscle.processed.mod,by=c("index","node"))
nodes.contrast.3RWGD.muscle.jumps.15spe.removing.jump.node <- nodes.contrast.3RWGD.muscle.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.3RWGD.muscle.jumps.15spe.removing.jump.node$Exp.abs)),]
nodes.contrast.3RWGD.muscle.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.3RWGD.muscle.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication","FishWGD"))
nodes.contrast.3RWGD.muscle.jumps.15spe.removing.jump.node$pic<-abs(nodes.contrast.3RWGD.muscle.jumps.15spe.removing.jump.node$Exp.abs)

########### Testing OC for all nodes of 139 3RWGD exp trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.3RWGD.muscle.nojump<-summary.tree(Calibrated.3RWGD.muscle.nojump)
summary.3RWGD.muscle.nojump$Exp.abs <- abs(summary.3RWGD.muscle.nojump$pic.muscle)
summary.3RWGD.muscle.nojumps<-summary.3RWGD.muscle.nojump[,!(names(summary.3RWGD.muscle.nojump) %in% c("pic.muscle"))]
nodes.contrast.3RWGD.muscle.nojumps.15spe <- summary.3RWGD.muscle.nojumps[which(!is.na(summary.3RWGD.muscle.nojumps$Exp.abs)),]
nodes.contrast.3RWGD.muscle.nojumps.15spe$Event <- factor(nodes.contrast.3RWGD.muscle.nojumps.15spe$events, levels=c("speciation", "duplication","FishWGD"))
rm(summary.3RWGD.muscle.nojump)
rm(summary.3RWGD.muscle.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.3RWGD.muscle.new.fish<-test.output.3RWGD.muscle.processed.mod[which(test.output.3RWGD.muscle.processed.mod$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##1015 obs

## Step 6: Analysis on fish specific clades
test.output.3RWGD.muscle.processed.new.fish<-unique(test.output.3RWGD.muscle.new.fish[c(6,3,7,8,10,2,1)])  #858 obs

########### Proportion of jump events in teleosts in the 586 WGD trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.3RWGD.fish.muscle<-NULL
info.jump.3RWGD.fish.muscle<-bind_rows(lapply(Calibrated.3RWGD.muscle.jump, tree.data.collection.3RWGD.fish))
info.jump.3RWGD.fish.muscle$spe.num<-info.jump.3RWGD.fish.muscle$internal.events.fish-(info.jump.3RWGD.fish.muscle$dup.num + info.jump.3RWGD.fish.muscle$WGD.num)
info.jump.3RWGD.fish.muscle<-info.jump.3RWGD.fish.muscle[c(11,14,3,12,5,4,13)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.3RWGD.fish.muscle<-
  test.output.3RWGD.muscle.new.fish %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.3RWGD.fish.muscle)<-c("index.old","jump.dup")

test.WGD.3RWGD.fish.muscle<-
  test.output.3RWGD.muscle.new.fish %>%
  filter(node.event=="FishWGD" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.WGD.3RWGD.fish.muscle)<-c("index.old","jump.WGD")

test.spe.3RWGD.fish.muscle<-
  test.output.3RWGD.muscle.new.fish %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.3RWGD.fish.muscle)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.3RWGD.muscle.fish1<-full_join(test.spe.3RWGD.fish.muscle,test.dup.3RWGD.fish.muscle)
merged.3RWGD.muscle.fish<-full_join(merged.3RWGD.muscle.fish1,test.WGD.3RWGD.fish.muscle)
info.jump.muscle.final.3RWGD.fish<-merge(info.jump.3RWGD.fish.muscle,merged.3RWGD.muscle.fish, by = c("index.old"))
info.jump.muscle.final.3RWGD.fish[is.na(info.jump.muscle.final.3RWGD.fish)] <- 0
rm(merged.3RWGD.muscle.fish1)

## Have to consider trees where atleast one WGD and one speciation exists
info.jump.muscle.final1.3RWGD.fish<-info.jump.muscle.final.3RWGD.fish[which(info.jump.muscle.final.3RWGD.fish$spe.num>0 & info.jump.muscle.final.3RWGD.fish$WGD.num>0),]
rm(info.jump.muscle.final.3RWGD.fish)

## Proportion of jumps
info.jump.muscle.final1.3RWGD.fish$dup.jump.prob<-(info.jump.muscle.final1.3RWGD.fish$jump.dup/(info.jump.muscle.final1.3RWGD.fish$dup.num*2))
info.jump.muscle.final1.3RWGD.fish$spe.jump.prob<-(info.jump.muscle.final1.3RWGD.fish$jump.spe/(info.jump.muscle.final1.3RWGD.fish$spe.num*2))
info.jump.muscle.final1.3RWGD.fish$WGD.jump.prob<-(info.jump.muscle.final1.3RWGD.fish$jump.WGD/(info.jump.muscle.final1.3RWGD.fish$WGD.num*2))

## Since,many WGD trees didnot have a duplication node, it produces NA for dup.jump.prob
## To get rid of this
info.jump.muscle.final2.3RWGD.fish<-info.jump.muscle.final1.3RWGD.fish[complete.cases(info.jump.muscle.final1.3RWGD.fish),] #97 obs

## Stats
median(na.omit(info.jump.muscle.final2.3RWGD.fish$dup.jump.prob)) #0
median(na.omit(info.jump.muscle.final2.3RWGD.fish$spe.jump.prob)) #0.1
median(na.omit(info.jump.muscle.final2.3RWGD.fish$WGD.jump.prob)) #0
mean(na.omit(info.jump.muscle.final2.3RWGD.fish$dup.jump.prob)) #0.1669
mean(na.omit(info.jump.muscle.final2.3RWGD.fish$spe.jump.prob)) #0.1856
mean(na.omit(info.jump.muscle.final2.3RWGD.fish$WGD.jump.prob)) #0.0989

pval.jump.3RWGD.muscle.fish.spe2dup <- paired.wilcox(info.jump.muscle.final2.3RWGD.fish$spe.jump.prob,info.jump.muscle.final2.3RWGD.fish$dup.jump.prob) #p-value =2.68e-1
pval.jump.3RWGD.muscle.fish.spe2WGD <- paired.wilcox(info.jump.muscle.final2.3RWGD.fish$spe.jump.prob,info.jump.muscle.final2.3RWGD.fish$WGD.jump.prob) #p-value = 5.02e-7

## For median data
dtfr.3RWGD.muscle.fish<-NULL
dtfr.3RWGD.muscle.fish<-data.frame(prop.jump=info.jump.muscle.final2.3RWGD.fish$spe.jump.prob,events="speciation")
dtfr.3RWGD.muscle.fish<-rbind(dtfr.3RWGD.muscle.fish,data.frame(prop.jump=info.jump.muscle.final2.3RWGD.fish$dup.jump.prob,events="duplication"))
dtfr.3RWGD.muscle.fish<-rbind(dtfr.3RWGD.muscle.fish,data.frame(prop.jump=info.jump.muscle.final2.3RWGD.fish$WGD.jump.prob,events="FishWGD"))
dtfr.3RWGD.muscle.fish$Exp.abs <-abs(dtfr.3RWGD.muscle.fish$prop.jump)
dtfr.3RWGD.muscle.fish$Event<-factor(dtfr.3RWGD.muscle.fish$events, levels=c("speciation","duplication","FishWGD"))


median.data.jump.3RWGD.muscle.fish<-NULL
median.data.jump.3RWGD.muscle.fish<- dtfr.3RWGD.muscle.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.3RWGD.muscle.fish$pic.round<-paste0("Median = ",round(median.data.jump.3RWGD.muscle.fish$muscle.abs,4), sep="")
median.data.jump.3RWGD.muscle.fish$count<- paste0(median.data.jump.3RWGD.muscle.fish$Event,"\n","(n = ",median.data.jump.3RWGD.muscle.fish$freq,")")

plot3R2<-boxplot.new.3RWGD2(dtfr.3RWGD.muscle.fish,pval.jump.3RWGD.muscle.fish.spe2dup,pval.jump.3RWGD.muscle.fish.spe2WGD,"Exp.3RWGD",median.data.jump.3RWGD.muscle.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("Jump in exp levels for 5 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.05, xmax = 1.7, ymin = 0.67, ymax =0.67, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.4, y = 0.75, label= pval.jump.3RWGD.muscle.fish.spe2dup, fontface = 4)+
  annotate("rect", xmin = 1.05, xmax = 2.8, ymin = 0.88, ymax =0.88, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.7, y = 0.95, label= pval.jump.3RWGD.muscle.fish.spe2WGD, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.3RWGD.muscle.fish$count) ## check

########### Testing OC for all fish specific nodes of 586 WGD trees supporting jump##################
nodes.contrast.3RWGD.muscle.jumps.fish<-nodes.contrast.3RWGD.muscle.jumps.15spe[which(nodes.contrast.3RWGD.muscle.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 586 WGD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.3RWGD.muscle.jumps.fish.removing.jump.node<-nodes.contrast.3RWGD.muscle.jumps.15spe.removing.jump.node[which(nodes.contrast.3RWGD.muscle.jumps.15spe.removing.jump.node$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

########### Testing OC for all fish specific nodes of 139 3RWGD trees that do not support trait jump for teleosts ##################
nodes.contrast.3RWGD.muscle.nojumps.fish<-nodes.contrast.3RWGD.muscle.nojumps.15spe[which(nodes.contrast.3RWGD.muscle.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")

################## Chunk 5: Analysis on contrasts standardized sure SSD trees for traits ##########################

### Step 1a: Considering (sure.SSD.trees) 4139 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for muscle tissue for sure SSD trees 
count<-0
Calibrated.muscle.SSD1<-lapply(sure.SSD.trees, diagnostic.plot.test.exp.all.tissue, "TPM.muscle")
Calibrated.muscle.SSD2<-Calibrated.muscle.SSD1[!is.na(Calibrated.muscle.SSD1)] 
Calibrated.standardized.SSD.muscle<-Calibrated.muscle.SSD2[ ! sapply(Calibrated.muscle.SSD2, is.null) ]## finally we obtained 2792/4139 trees passing the diagnostic tests
rm(Calibrated.muscle.SSD1)
rm(Calibrated.muscle.SSD2)

### Step 1b: Considering (Calibrated.standardized.SSD.muscle) 2792 trees passing the diagnostic test for muscle expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.standardized.SSD.muscle,tree.index.collect))
rm(count)

test.output.SSD.muscle<-NULL
test.output.SSD.muscle<-merge(index.info,test.output.muscle1,by="Gene") ## 9899 observations
index.trees.jump<-unique(test.output.SSD.muscle$tree.num) ##2319 trees
Calibrated.SSD.muscle.jump<-Calibrated.standardized.SSD.muscle[c(index.trees.jump)] ##2319/2792=83.06% Brownian trees supported trait jump model for muscle expressions
Calibrated.SSD.muscle.nojump<-Calibrated.standardized.SSD.muscle[-c(index.trees.jump)] ## 473/2792=16.94%  trees with no support for jump in tau
rm(index.info)
rm(index.trees.jump)

### Step2: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.SSD.muscle))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.SSD.muscle)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 13748 obs
}

test.output.SSD.muscle.processed<-bind_cols(test.output.SSD.muscle,processed.fout) ## 13748 obs for all the 15 species
test.output.SSD.muscle.processed<-test.output.SSD.muscle.processed[(which(!(test.output.SSD.muscle.processed$clade.name %in% c("Clupeocephala","Osteoglossocephalai")))),] ##9466 obs
rm(i)
rm(test.output.SSD.muscle)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

## Step 3: Analysis on all 15 vertebrates clades
test.output.SSD.muscle.processed.mod<-unique(test.output.SSD.muscle.processed[c(6,3,7,8,10,2,1)])  #7601 obs

########### Proportion of jump events in 15 vertebrates species in the 2319 SSD trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.SSD.muscle<-NULL
count<-0
info.jump.SSD.muscle<-bind_rows(lapply(Calibrated.SSD.muscle.jump, tree.data.collection))
info.jump.SSD.muscle$spe.num<-info.jump.SSD.muscle$internal.events-info.jump.SSD.muscle$dup.num
info.jump.SSD.muscle<-info.jump.SSD.muscle[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.SSD.muscle<-
  test.output.SSD.muscle.processed %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.SSD.muscle)<-c("index.old","jump.dup")

test.spe.SSD.muscle<-
  test.output.SSD.muscle.processed %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.SSD.muscle)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.SSD.muscle<-full_join(test.spe.SSD.muscle,test.dup.SSD.muscle)
info.jump.muscle.final.SSD.all<-merge(info.jump.SSD.muscle,merged.SSD.muscle, by = c("index.old"))
info.jump.muscle.final.SSD.all[is.na(info.jump.muscle.final.SSD.all)] <- 0

## Have to consider trees where atleast one duplication and speciation exists
info.jump.muscle.final1.SSD.all<-info.jump.muscle.final.SSD.all[which(info.jump.muscle.final.SSD.all$dup.num>0 & info.jump.muscle.final.SSD.all$spe.num>0),]
rm(info.jump.muscle.final.SSD.all)

## Proportion of jumps
info.jump.muscle.final1.SSD.all$dup.jump.prob<-(info.jump.muscle.final1.SSD.all$jump.dup/(info.jump.muscle.final1.SSD.all$dup.num*2))
info.jump.muscle.final1.SSD.all$spe.jump.prob<-(info.jump.muscle.final1.SSD.all$jump.spe/(info.jump.muscle.final1.SSD.all$spe.num*2))

## Stats
median(na.omit(info.jump.muscle.final1.SSD.all$dup.jump.prob)) #0
median(na.omit(info.jump.muscle.final1.SSD.all$spe.jump.prob)) #0.0454
mean(na.omit(info.jump.muscle.final1.SSD.all$dup.jump.prob)) #0.2005
mean(na.omit(info.jump.muscle.final1.SSD.all$spe.jump.prob)) #0.1033

pval.jump.SSD.muscle <- paired.wilcox(info.jump.muscle.final1.SSD.all$spe.jump.prob,info.jump.muscle.final1.SSD.all$dup.jump.prob) #p-value = 2.2e-16


########### Testing OC for all nodes of 2319 SSD trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.SSD.muscle.jump<-summary.tree(Calibrated.SSD.muscle.jump)
summary.SSD.muscle.jump$Exp.abs <- abs(summary.SSD.muscle.jump$pic.muscle)
summary.SSD.muscle.jumps<-summary.SSD.muscle.jump[,!(names(summary.SSD.muscle.jump) %in% c("pic.muscle"))]
nodes.contrast.SSD.muscle.jumps.15spe <- summary.SSD.muscle.jumps[which(!is.na(summary.SSD.muscle.jumps$Exp.abs)),]
nodes.contrast.SSD.muscle.jumps.15spe$Event <- factor(nodes.contrast.SSD.muscle.jumps.15spe$events, levels=c("speciation", "duplication"))
nodes.contrast.SSD.muscle.jumps.15spe<-nodes.contrast.SSD.muscle.jumps.15spe[(which(!(nodes.contrast.SSD.muscle.jumps.15spe$label %in% c("Clupeocephala","Osteoglossocephalai")))),]
rm(summary.SSD.muscle.jump)
rm(summary.SSD.muscle.jumps)

############ Testing OC using all 15 species for 2319 SSD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.SSD.muscle.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.SSD.muscle.jumps.15spe$index<-nodes.contrast.SSD.muscle.jumps.15spe$index.tree
test.output.SSD.muscle.processed$node<-test.output.SSD.muscle.processed$From
nodes.contrast.SSD.muscle.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.SSD.muscle.jumps.15spe,test.output.SSD.muscle.processed,by=c("index","node"))
nodes.contrast.SSD.muscle.jumps.15spe.removing.jump.node <- nodes.contrast.SSD.muscle.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.SSD.muscle.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.SSD.muscle.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.SSD.muscle.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.SSD.muscle.jumps.15spe.removing.jump.node$Exp.abs<-nodes.contrast.SSD.muscle.jumps.15spe.removing.jump.node$pic
nodes.contrast.SSD.muscle.jumps.15spe.removing.jump.node<-nodes.contrast.SSD.muscle.jumps.15spe.removing.jump.node[(which(!(nodes.contrast.SSD.muscle.jumps.15spe.removing.jump.node$label %in% c("Clupeocephala","Osteoglossocephalai")))),]

########### Testing OC for all nodes of 473 SSD exp trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.SSD.muscle.nojump<-summary.tree(Calibrated.SSD.muscle.nojump)
summary.SSD.muscle.nojump$Exp.abs <- abs(summary.SSD.muscle.nojump$pic)
summary.SSD.muscle.nojumps<-summary.SSD.muscle.nojump[,!(names(summary.SSD.muscle.nojump) %in% c("pic"))]
nodes.contrast.SSD.muscle.nojumps.15spe <- summary.SSD.muscle.nojumps[which(!is.na(summary.SSD.muscle.nojumps$Exp.abs)),]
nodes.contrast.SSD.muscle.nojumps.15spe$Event <- factor(nodes.contrast.SSD.muscle.nojumps.15spe$events, levels=c("speciation", "duplication"))
nodes.contrast.SSD.muscle.nojumps.15spe<-nodes.contrast.SSD.muscle.nojumps.15spe[(which(!(nodes.contrast.SSD.muscle.nojumps.15spe$label %in% c("Clupeocephala","Osteoglossocephalai")))),] 
rm(summary.SSD.muscle.nojump)
rm(summary.SSD.muscle.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.SSD.muscle.new.fish<-test.output.SSD.muscle.processed[which(test.output.SSD.muscle.processed$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##3013 obs

## Step 6: Analysis on fish specific clades
test.output.SSD.muscle.new.fish.mod<-unique(test.output.SSD.muscle.new.fish[c(6,3,7,8,10,2,1)])  #2344 obs


########### Proportion of jump events in teleosts in the 2319 SSD trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.SSD.muscle.trees.fish<-NULL
count<-0
info.jump.SSD.muscle.trees.fish<-bind_rows(lapply(Calibrated.SSD.muscle.jump, tree.data.collection.fish))
info.jump.SSD.muscle.trees.fish$spe.num<-info.jump.SSD.muscle.trees.fish$internal.events.fish-info.jump.SSD.muscle.trees.fish$dup.num
info.jump.SSD.muscle.trees.fish<-info.jump.SSD.muscle.trees.fish[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.SSD.muscle.fish<-
  test.output.SSD.muscle.new.fish  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.SSD.muscle.fish)<-c("index.old","jump.dup")


test.spe.SSD.muscle.fish<-
  test.output.SSD.muscle.new.fish %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.SSD.muscle.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.SSD.muscle.fish<-full_join(test.spe.SSD.muscle.fish,test.dup.SSD.muscle.fish)
info.jump.trees.SSD.muscle.final.fish<-merge(info.jump.SSD.muscle.trees.fish,merged.SSD.muscle.fish, by = c("index.old"))
info.jump.trees.SSD.muscle.final.fish[is.na(info.jump.trees.SSD.muscle.final.fish)] <- 0

## Have to consider trees where atleast one fish specific duplication and speciation exists
info.jump.trees.SSD.muscle.final1.fish<-info.jump.trees.SSD.muscle.final.fish[which(info.jump.trees.SSD.muscle.final.fish$dup.num>0 & info.jump.trees.SSD.muscle.final.fish$spe.num>0),]

## Proportion of jumps
info.jump.trees.SSD.muscle.final1.fish$dup.jump.prob<-(info.jump.trees.SSD.muscle.final1.fish$jump.dup/(info.jump.trees.SSD.muscle.final1.fish$dup.num*2))
info.jump.trees.SSD.muscle.final1.fish$spe.jump.prob<-(info.jump.trees.SSD.muscle.final1.fish$jump.spe/(info.jump.trees.SSD.muscle.final1.fish$spe.num*2))


median(na.omit(info.jump.trees.SSD.muscle.final1.fish$dup.jump.prob)) #0.167
median(na.omit(info.jump.trees.SSD.muscle.final1.fish$spe.jump.prob)) #0.125
mean(na.omit(info.jump.trees.SSD.muscle.final1.fish$dup.jump.prob)) #0.2547
mean(na.omit(info.jump.trees.SSD.muscle.final1.fish$spe.jump.prob)) #0.1718

########### Testing OC for all fish specific nodes of 2319 SSD trees supporting jump##################
nodes.contrast.SSD.muscle.jumps.fish<-nodes.contrast.SSD.muscle.jumps.15spe[which(nodes.contrast.SSD.muscle.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] #9206 obs

############ Testing OC for teleosts in 2319 SSD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 

nodes.contrast.SSD.muscle.jumps.fish.removing.jump.node <- NULL
nodes.contrast.SSD.muscle.jumps.fish$index<-nodes.contrast.SSD.muscle.jumps.fish$index.tree
test.output.SSD.muscle.new.fish$node<-test.output.SSD.muscle.new.fish$From
nodes.contrast.SSD.muscle.jumps.fish.removing.jump.node <-anti_join(nodes.contrast.SSD.muscle.jumps.fish,test.output.SSD.muscle.new.fish,by=c("index","node"))
nodes.contrast.SSD.muscle.jumps.fish.removing.jump.node <- nodes.contrast.SSD.muscle.jumps.fish.removing.jump.node[which(!is.na(nodes.contrast.SSD.muscle.jumps.fish.removing.jump.node$pic)),]
nodes.contrast.SSD.muscle.jumps.fish.removing.jump.node$Event <- factor(nodes.contrast.SSD.muscle.jumps.fish.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.SSD.muscle.jumps.fish.removing.jump.node$Exp.abs<-nodes.contrast.SSD.muscle.jumps.fish.removing.jump.node$pic
nodes.contrast.SSD.muscle.jumps.fish.removing.jump.node<-nodes.contrast.SSD.muscle.jumps.fish.removing.jump.node[(which(!(nodes.contrast.SSD.muscle.jumps.fish.removing.jump.node$label %in% c("Clupeocephala","Osteoglossocephalai")))),]

########### Testing OC for all fish specific nodes of 473 SSD exp trees that do not support trait jump for teleosts ##################
nodes.contrast.SSD.muscle.nojumps.fish<-nodes.contrast.SSD.muscle.nojumps.15spe[which(nodes.contrast.SSD.muscle.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] #1735 obs

#save.image("norm_exp_15spe_4tips_new.RData")

######## Chunk6: Testing the ortholog conjecture for all the contrast standardized SSD trees for muscle ###########

########### Testing OC for all nodes of 2792 trees supporting jump for 15 species ##################

## Calculation of nodes contrasts for muscle
Count<-0
Calibrated.standardized.SSD.muscle<-lapply(Calibrated.standardized.SSD.muscle,contrast.calc,trait="TPM.muscle")

## Summarizing contrasts for normalized exp levels 
summary.muscle.SSD<-summary.tree(Calibrated.standardized.SSD.muscle)
summary.muscle.SSD$Exp.abs <- abs(summary.muscle.SSD$pic_TPM.muscle)
#summary.muscle.all<-summary.muscle.all[,!(names(summary.muscle.all) %in% c("pic_TPM.muscle"))]
nodes.contrast.muscle.SSD.15spe <- summary.muscle.SSD[which(!is.na(summary.muscle.SSD$Exp.abs)),]
nodes.contrast.muscle.SSD.15spe$Event <- factor(nodes.contrast.muscle.SSD.15spe$events, levels=c("speciation", "duplication"))
rm(summary.muscle.SSD)

########### Testing OC for all nodes of 4922 trees supporting jump for 5 teleosts ##################
nodes.contrast.muscle.SSD.fish<-nodes.contrast.muscle.SSD.15spe[which(nodes.contrast.muscle.SSD.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")
