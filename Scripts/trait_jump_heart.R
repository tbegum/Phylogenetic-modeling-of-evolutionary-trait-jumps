

### This code is written to work with "trait jump model" for heart

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

test.output.heart<-read.csv(paste0(folder1,"To_pablo/heart/results_heart.csv", sep=""),sep=",",header=T) ##183764 obs

## Only considering trees with posterior probability of jump more than 70%
test.output.heart1<-test.output.heart[which(test.output.heart$JumpProbability>=0.7),] ##34343 obs

########################## Chunk 3: Analysis on all contrasts standardized trees for traits ###############################
### Step 1a: Considering (trees.all.interest) 6923 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for heart tissue for all trees 
count<-0
Calibrated.heart.exp1<-lapply(trees.all.interest, diagnostic.plot.test.exp.all.tissue, "TPM.heart")
Calibrated.heart.exp2<-Calibrated.heart.exp1[!is.na(Calibrated.heart.exp1)] 
Calibrated.all.standardized.heart<-Calibrated.heart.exp2[ ! sapply(Calibrated.heart.exp2, is.null) ]## finally we obtained 4922/6923 tree data passing the diagnostic tests
rm(Calibrated.heart.exp1)
rm(Calibrated.heart.exp2)

### Step 1b: Considering (Calibrated.all.standardized.heart) 4922 tree data passing the diagnostic test for heart expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.all.standardized.heart,tree.index.collect))
rm(count)

test.output.heart.all<-NULL
test.output.heart.all<-merge(index.info,test.output.heart1,by="Gene") ## 21836 observations
index.trees.jump<-unique(test.output.heart.all$tree.num) ##3495 trees
Calibrated.standardized.heart.jump<-Calibrated.all.standardized.heart[c(index.trees.jump)] ##3495/4922=71% Brownian trees supported trait jump model for heart expressions 
Calibrated.standardized.heart.nojump<-Calibrated.all.standardized.heart[-c(index.trees.jump)] ## 1427/4922=29% trees with no support for jump in heart expression
rm(index.info)
rm(index.trees.jump)

### Step2: Calculating PICs of required trait for trees to match the index number of the trees on which levolution was run
Count<-0
trees.all.interest<-lapply(trees.all.interest,contrast.calc,trait="TPM.heart")
# Count<-0
# Calibrated.standardized.heart.jump<-lapply(Calibrated.standardized.heart.jump,contrast.calc,trait="TPM.heart") 
# rm(Count)

### Step3: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.heart.all))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.heart.all)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 21836 obs
}
test.output.heart.processed<-bind_cols(test.output.heart.all,processed.fout) ## 21836 obs for all the 15 species
rm(i)
rm(test.output.heart.all)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

## Step 4: Analysis on all 15 vertebrates clades
test.output.heart.new.mod<-unique(test.output.heart.processed[c(6,3,7,8,10,2,1)])  #16995 obs

########### Proportion of jump events in 15 vertebrates species in the 3495 trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.trees.heart<-NULL
count<-0
info.jump.trees.heart<-bind_rows(lapply(Calibrated.standardized.heart.jump, tree.data.collection))
info.jump.trees.heart$spe.num<-info.jump.trees.heart$internal.events-info.jump.trees.heart$dup.num
info.jump.trees.heart<-info.jump.trees.heart[c(11,12,3:5)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.heart<-
  test.output.heart.processed  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.heart)<-c("index.old","jump.dup")

test.spe.heart<-
  test.output.heart.processed %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.heart)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.heart<-full_join(test.spe.heart,test.dup.heart)
info.jump.trees.heart.final<-merge(info.jump.trees.heart,merged.heart, by = c("index.old"))
info.jump.trees.heart.final[is.na(info.jump.trees.heart.final)] <- 0

## Have to consider trees where atleast one duplication and speciation exists
info.jump.trees.heart.final1<-info.jump.trees.heart.final[which(info.jump.trees.heart.final$dup.num>0 & info.jump.trees.heart.final$spe.num>0),]
rm(info.jump.trees.heart.final)

## Proportion of jumps
info.jump.trees.heart.final1$dup.jump.prob<-(info.jump.trees.heart.final1$jump.dup/(info.jump.trees.heart.final1$dup.num*2))
info.jump.trees.heart.final1$spe.jump.prob<-(info.jump.trees.heart.final1$jump.spe/(info.jump.trees.heart.final1$spe.num*2))

## Stats
median(na.omit(info.jump.trees.heart.final1$dup.jump.prob)) #0
median(na.omit(info.jump.trees.heart.final1$spe.jump.prob)) #0.083
mean(na.omit(info.jump.trees.heart.final1$dup.jump.prob)) #0.2449
mean(na.omit(info.jump.trees.heart.final1$spe.jump.prob)) #0.1487

pval.jump.heart <- paired.wilcox(info.jump.trees.heart.final1$spe.jump.prob,info.jump.trees.heart.final1$dup.jump.prob) #p-value = 2.2e-16

## For median data
dtfr.heart<-NULL
dtfr.heart<-data.frame(prop.jump=info.jump.trees.heart.final1$spe.jump.prob,events="speciation")
dtfr.heart<-rbind(dtfr.heart,data.frame(prop.jump=info.jump.trees.heart.final1$dup.jump.prob,events="duplication"))
dtfr.heart$Exp.abs <-abs(dtfr.heart$prop.jump)
dtfr.heart$Event<-factor(dtfr.heart$events, levels=c("speciation","duplication"))

library(tidyverse)
median.data.jump.heart<-NULL
median.data.jump.heart<- dtfr.heart%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.heart$pic.round<-paste0("Median = ",round(median.data.jump.heart$Exp.abs,4), sep="")
median.data.jump.heart$count<- paste0(median.data.jump.heart$Event,"\n\n","(n = ",median.data.jump.heart$freq,")")

plot4A<-boxplot.new2(dtfr.heart,pval.jump.heart,"Exp",median.data.jump.heart) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.heart, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.heart$count) ## check

# Create Data
data.jump.heart <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(3495,1427)
)

# Compute the position of labels
data.jump.heart <- data.jump.heart %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.heart$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.heart$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.heart<-ggplot(data.jump.heart, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 3495 trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.heart.jump<-summary.tree(Calibrated.standardized.heart.jump)
summary.heart.jump$Exp.abs <- abs(summary.heart.jump$pic.heart)
summary.heart.jumps<-summary.heart.jump[,!(names(summary.heart.jump) %in% c("pic.heart"))]
nodes.contrast.heart.jumps.15spe <- summary.heart.jumps[which(!is.na(summary.heart.jumps$Exp.abs)),]
nodes.contrast.heart.jumps.15spe$Event <- factor(nodes.contrast.heart.jumps.15spe$events, levels=c("speciation", "duplication"))
rm(summary.heart.jump)
rm(summary.heart.jumps)

############ Testing OC using all 15 species for 3495 trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.heart.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.heart.jumps.15spe$index<-nodes.contrast.heart.jumps.15spe$index.tree
test.output.heart.processed$node<-test.output.heart.processed$From
nodes.contrast.heart.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.heart.jumps.15spe,test.output.heart.processed,by=c("index","node"))
nodes.contrast.heart.jumps.15spe.removing.jump.node <- nodes.contrast.heart.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.heart.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.heart.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.heart.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.heart.jumps.15spe.removing.jump.node$Exp.abs<-nodes.contrast.heart.jumps.15spe.removing.jump.node$pic


########### Testing OC for all nodes of 1427 trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.heart.nojump<-summary.tree(Calibrated.standardized.heart.nojump)
summary.heart.nojump$Exp.abs <- abs(summary.heart.nojump$pic.heart)
summary.heart.nojumps<-summary.heart.nojump[,!(names(summary.heart.nojump) %in% c("pic.heart"))]
nodes.contrast.heart.nojumps.15spe <- summary.heart.nojump[which(!is.na(summary.heart.nojump$Exp.abs)),]
nodes.contrast.heart.nojumps.15spe$Event <- factor(nodes.contrast.heart.nojumps.15spe$events, levels=c("speciation", "duplication"))
rm(summary.heart.nojump)
rm(summary.heart.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.heart.new.fish<-test.output.heart.processed[which(test.output.heart.processed$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##8102 obs

## Step 6: Analysis on fish specific clades
test.output.heart.new.fish.mod<-unique(test.output.heart.new.fish[c(6,3,7,8,10,2,1)])  #6193 obs

########### Proportion of jump events in teleosts in the 3495 trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.trees.heart.fish<-NULL
count<-0
info.jump.trees.heart.fish<-bind_rows(lapply(Calibrated.standardized.heart.jump, tree.data.collection.fish))
info.jump.trees.heart.fish$spe.num<-info.jump.trees.heart.fish$internal.events.fish-info.jump.trees.heart.fish$dup.num
info.jump.trees.heart.fish<-info.jump.trees.heart.fish[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.heart.fish<-
  test.output.heart.new.fish  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.heart.fish)<-c("index.old","jump.dup")

test.spe.heart.fish<-
  test.output.heart.new.fish %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.heart.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.heart.fish<-full_join(test.spe.heart.fish,test.dup.heart.fish)
info.jump.trees.heart.final.fish<-merge(info.jump.trees.heart.fish,merged.heart.fish, by = c("index.old"))
info.jump.trees.heart.final.fish[is.na(info.jump.trees.heart.final.fish)] <- 0

## Have to consider trees where atleast one fish specific duplication and speciation exists
info.jump.trees.heart.final1.fish<-info.jump.trees.heart.final.fish[which(info.jump.trees.heart.final.fish$dup.num>0 & info.jump.trees.heart.final.fish$spe.num>0),]

## Proportion of jumps
info.jump.trees.heart.final1.fish$dup.jump.prob<-(info.jump.trees.heart.final1.fish$jump.dup/(info.jump.trees.heart.final1.fish$dup.num*2))
info.jump.trees.heart.final1.fish$spe.jump.prob<-(info.jump.trees.heart.final1.fish$jump.spe/(info.jump.trees.heart.final1.fish$spe.num*2))


median(na.omit(info.jump.trees.heart.final1.fish$dup.jump.prob)) #0.0
median(na.omit(info.jump.trees.heart.final1.fish$spe.jump.prob)) #0.1667
mean(na.omit(info.jump.trees.heart.final1.fish$dup.jump.prob)) #0.2263
mean(na.omit(info.jump.trees.heart.final1.fish$spe.jump.prob)) #0.2467

pval.jump.heart.fish <- paired.wilcox(info.jump.trees.heart.final1.fish$spe.jump.prob,info.jump.trees.heart.final1.fish$dup.jump.prob) #p-value = 2.02e-7


## For median data
dtfr.heart.fish<-NULL
dtfr.heart.fish<-data.frame(prop.jump=info.jump.trees.heart.final1.fish$spe.jump.prob,events="speciation")
dtfr.heart.fish<-rbind(dtfr.heart.fish,data.frame(prop.jump=info.jump.trees.heart.final1.fish$dup.jump.prob,events="duplication"))
dtfr.heart.fish$Exp.abs <-abs(dtfr.heart.fish$prop.jump)
dtfr.heart.fish$Event<-factor(dtfr.heart.fish$events, levels=c("speciation","duplication"))


library(tidyverse)
median.data.jump.heart.heart.fish<-NULL
median.data.jump.heart.heart.fish<- dtfr.heart.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.heart.heart.fish$pic.round<-paste0("Median = ",round(median.data.jump.heart.heart.fish$Exp.abs,4), sep="")
median.data.jump.heart.heart.fish$count<- paste0(median.data.jump.heart.heart.fish$Event,"\n\n","(n = ",median.data.jump.heart.heart.fish$freq,")")

plot4B<-boxplot.new2(dtfr.heart.fish,pval.jump.heart.fish,"Exp",median.data.jump.heart.heart.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 5 teleosts"))))+
  #ylab(expression(bold(paste("Proportions of jump events for teleosts ")))) +
  ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.heart.fish, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.heart.heart.fish$count) ## check

########### Testing OC for all fish specific nodes of 3495 trees supporting jump##################
nodes.contrast.heart.jumps.fish<-nodes.contrast.heart.jumps.15spe[which(nodes.contrast.heart.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 3495 trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.heart.jumps.fish.removing.jump.node <- NULL
nodes.contrast.heart.jumps.fish$index<-nodes.contrast.heart.jumps.fish$index.tree
test.output.heart.new.fish$node<-test.output.heart.new.fish$From
nodes.contrast.heart.jumps.fish.removing.jump.node <-anti_join(nodes.contrast.heart.jumps.fish,test.output.heart.new.fish,by=c("index","node"))
nodes.contrast.heart.jumps.fish.removing.jump.node <- nodes.contrast.heart.jumps.fish.removing.jump.node[which(!is.na(nodes.contrast.heart.jumps.fish.removing.jump.node$pic)),]
nodes.contrast.heart.jumps.fish.removing.jump.node$Event <- factor(nodes.contrast.heart.jumps.fish.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.heart.jumps.fish.removing.jump.node$Exp.abs<-nodes.contrast.heart.jumps.fish.removing.jump.node$pic

rm(test.output.heart.new.fish)

########### Testing OC for all fish specific nodes of 1427 trees that do not support trait jump for teleosts ##################
nodes.contrast.heart.nojumps.fish<-nodes.contrast.heart.nojumps.15spe[which(nodes.contrast.heart.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")

################## Chunk 4: Analysis on contrasts standardized sure 3RWGD trees for traits ##########################

### Step 1a: Considering (Sure.WGD.trees) 1159 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for heart tissue for sure WGD trees 
count<-0
Calibrated.heart.3Rexp1<-lapply(Sure.WGD.trees, diagnostic.plot.test.exp.all.tissue, "TPM.heart")
Calibrated.heart.3Rexp2<-Calibrated.heart.3Rexp1[!is.na(Calibrated.heart.3Rexp1)] 
Calibrated.3R.standardized.heart<-Calibrated.heart.3Rexp2[ ! sapply(Calibrated.heart.3Rexp2, is.null) ]## finally we obtained 807/1159 trees passing the diagnostic tests
rm(Calibrated.heart.3Rexp1)
rm(Calibrated.heart.3Rexp2)

### Step 1b: Considering (Calibrated.3R.standardized.heart) 807 trees passing the diagnostic test for heart expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.3R.standardized.heart,tree.index.collect))
rm(count)

test.output.3Rexp.heart<-NULL
test.output.3Rexp.heart<-merge(index.info,test.output.heart1,by="Gene") ## 3525 observations
index.trees.jump<-unique(test.output.3Rexp.heart$tree.num) ##556 trees
Calibrated.3RWGD.heart.jump<-Calibrated.3R.standardized.heart[c(index.trees.jump)] ##556/807=68.89% Brownian trees supported trait jump model for heart expressions
Calibrated.3RWGD.heart.nojump<-Calibrated.3R.standardized.heart[-c(index.trees.jump)] ## 251/807=31.11% trees with no support for jump in trait
rm(index.info)
rm(index.trees.jump)

### Step2: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.3Rexp.heart))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.3Rexp.heart)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 3525 obs
}

test.output.3RWGD.heart.processed<-bind_cols(test.output.3Rexp.heart,processed.fout) ## 3525 obs for all the 15 species
rm(i)
rm(test.output.3Rexp.heart)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

##Since our calibrated tree shows both the FishWGD and duplication nodes, 
##we want to modify the node events of the jump output to separate out the effect of FishWGD from the SSDs in those trees.

jump.3RWGD.heart.trees.summary<-summary.tree(Calibrated.3RWGD.heart.jump)
jump.3RWGD.heart.trees.summary$Gene<-jump.3RWGD.heart.trees.summary$index.tree
jump.3RWGD.heart.trees.summary$From<-jump.3RWGD.heart.trees.summary$node
test.output.3RWGD.heart.processed.mod<-merge(test.output.3RWGD.heart.processed,jump.3RWGD.heart.trees.summary,by=c("Gene","From"))
test.output.3RWGD.heart.processed.mod<-test.output.3RWGD.heart.processed.mod[c(1,3,2,4:11,16)]
test.output.3RWGD.heart.processed.mod$node.event<-test.output.3RWGD.heart.processed.mod$events
colnames(test.output.3RWGD.heart.processed.mod)[10]<-"pic"
test.output.3RWGD.heart.processed.mod<-test.output.3RWGD.heart.processed.mod[c(1:11)] ## 3525 obs

## Step 3: Analysis on all 15 vertebrates clades
test.output.3RWGD.heart.processed.new<-unique(test.output.3RWGD.heart.processed.mod[c(6,3,7,8,10,2,1)])  #2813 obs

########### Proportion of jump events in 15 vertebrates species in the 556 3RWGD trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.3RWGD.heart<-NULL
info.jump.3RWGD.heart<-bind_rows(lapply(Calibrated.3RWGD.heart.jump, tree.data.collection.3RWGD))
info.jump.3RWGD.heart$spe.num<-info.jump.3RWGD.heart$internal.events-(info.jump.3RWGD.heart$dup.num + info.jump.3RWGD.heart$WGD.num)
info.jump.3RWGD.heart<-info.jump.3RWGD.heart[c(11,14,3,12,5,4,13)]


## Counting frquency of jump for "speciation" and "duplication" events
test.dup.3RWGD.heart<-
  test.output.3RWGD.heart.processed.mod %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.3RWGD.heart)<-c("index.old","jump.dup")

test.WGD.3RWGD.heart<-
  test.output.3RWGD.heart.processed.mod %>%
  filter(node.event=="FishWGD" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.WGD.3RWGD.heart)<-c("index.old","jump.WGD")

test.spe.3RWGD.heart<-
  test.output.3RWGD.heart.processed.mod %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.3RWGD.heart)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.3RWGD.heart1<-full_join(test.spe.3RWGD.heart,test.dup.3RWGD.heart)
merged.3RWGD.heart<-full_join(merged.3RWGD.heart1,test.WGD.3RWGD.heart)
info.jump.3RWGD.heart.final1<-merge(info.jump.3RWGD.heart,merged.3RWGD.heart, by = c("index.old"))
info.jump.3RWGD.heart.final1[is.na(info.jump.3RWGD.heart.final1)] <- 0
rm(merged.3RWGD.heart1)

## Have to consider trees where atleast one WGD and one speciation exists
info.jump.3RWGD.heart.final<-info.jump.3RWGD.heart.final1[which(info.jump.3RWGD.heart.final1$spe.num>0 & info.jump.3RWGD.heart.final1$WGD.num>0),]
rm(info.jump.3RWGD.heart.final1)

## Proportion of jumps
info.jump.3RWGD.heart.final$dup.jump.prob<-(info.jump.3RWGD.heart.final$jump.dup/(info.jump.3RWGD.heart.final$dup.num*2))
info.jump.3RWGD.heart.final$spe.jump.prob<-(info.jump.3RWGD.heart.final$jump.spe/(info.jump.3RWGD.heart.final$spe.num*2))
info.jump.3RWGD.heart.final$WGD.jump.prob<-(info.jump.3RWGD.heart.final$jump.WGD/(info.jump.3RWGD.heart.final$WGD.num*2))

## Since,many WGD trees didnot have a duplication node, it produces NA for dup.jump.prob
## To get rid of this
info.jump.3RWGD.heart.final2<-info.jump.3RWGD.heart.final[complete.cases(info.jump.3RWGD.heart.final),] ##227 obs

## Stats
median(info.jump.3RWGD.heart.final2$dup.jump.prob) #0.1
median(info.jump.3RWGD.heart.final2$spe.jump.prob) #0.0909
median(info.jump.3RWGD.heart.final2$WGD.jump.prob) #0
mean(info.jump.3RWGD.heart.final2$dup.jump.prob) #0.2965
mean(info.jump.3RWGD.heart.final2$spe.jump.prob) #0.167
mean(info.jump.3RWGD.heart.final2$WGD.jump.prob) #0.1596

pval.jump.3RWGD.heart.spe2dup <- paired.wilcox(info.jump.3RWGD.heart.final2$spe.jump.prob,info.jump.3RWGD.heart.final2$dup.jump.prob) #p-value = 2.536e-8
pval.jump.3RWGD.heart.spe2WGD <- paired.wilcox(info.jump.3RWGD.heart.final2$spe.jump.prob,info.jump.3RWGD.heart.final2$WGD.jump.prob) #p-value = 3.57e-7

## For median data
dtfr.3RWGD.heart.15spe<-NULL
dtfr.3RWGD.heart.15spe<-data.frame(prop.jump=info.jump.3RWGD.heart.final2$spe.jump.prob,events="speciation")
dtfr.3RWGD.heart.15spe<-rbind(dtfr.3RWGD.heart.15spe,data.frame(prop.jump=info.jump.3RWGD.heart.final2$dup.jump.prob,events="duplication"))
dtfr.3RWGD.heart.15spe<-rbind(dtfr.3RWGD.heart.15spe,data.frame(prop.jump=info.jump.3RWGD.heart.final2$WGD.jump.prob,events="FishWGD"))
dtfr.3RWGD.heart.15spe$heart.abs <-abs(dtfr.3RWGD.heart.15spe$prop.jump)
dtfr.3RWGD.heart.15spe$Event<-factor(dtfr.3RWGD.heart.15spe$events, levels=c("speciation","duplication","FishWGD"))


median.data.jump.3RWGD.heart.15spe<-NULL
median.data.jump.3RWGD.heart.15spe<- dtfr.3RWGD.heart.15spe%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.3RWGD.heart.15spe$pic.round<-paste0("Median = ",round(median.data.jump.3RWGD.heart.15spe$Exp.abs,4), sep="")
median.data.jump.3RWGD.heart.15spe$count<- paste0(median.data.jump.3RWGD.heart.15spe$Event,"\n\n","(n = ",median.data.jump.3RWGD.heart.15spe$freq,")")

plot3R1<-boxplot.new.3RWGD2(dtfr.3RWGD.heart.15spe,pval.jump.3RWGD.heart.spe2dup,pval.jump.3RWGD.heart.spe2WGD,"Exp.3RWGD",median.data.jump.3RWGD.heart.15spe) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.05, xmax = 1.7, ymin = 0.67, ymax =0.67, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.4, y = 0.75, label= pval.jump.3RWGD.heart.spe2dup, fontface = 4)+
  annotate("rect", xmin = 1.05, xmax = 2.8, ymin = 0.88, ymax =0.88, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.7, y = 0.95, label= pval.jump.3RWGD.heart.spe2WGD, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.3RWGD.heart.15spe$count) ## check

# Create Data
data.jump.3RWGD.exp <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(556,251)
)

# Compute the position of labels
data.jump.3RWGD.heart <- data.jump.3RWGD.heart %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.3RWGD.heart$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.3RWGD.heart$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.3RWGD.heart<-ggplot(data.jump.3RWGD.heart, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 556 WGD trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.3RWGD.heart.jump<-summary.tree(Calibrated.3RWGD.heart.jump)
summary.3RWGD.heart.jump$Exp.abs <- abs(summary.3RWGD.heart.jump$pic)
summary.3RWGD.heart.jumps<-summary.3RWGD.heart.jump[,!(names(summary.3RWGD.heart.jump) %in% c("pic"))]
nodes.contrast.3RWGD.heart.jumps.15spe <- summary.3RWGD.heart.jumps[which(!is.na(summary.3RWGD.heart.jumps$Exp.abs)),]
nodes.contrast.3RWGD.heart.jumps.15spe$Event <- factor(nodes.contrast.3RWGD.heart.jumps.15spe$events, levels=c("speciation", "duplication", "FishWGD"))
rm(summary.3RWGD.heart.jump)
rm(summary.3RWGD.heart.jumps)

############ Testing OC using all 15 species for 556 WGD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.3RWGD.heart.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.3RWGD.heart.jumps.15spe$index<-nodes.contrast.3RWGD.heart.jumps.15spe$index.tree
test.output.3RWGD.heart.processed.mod$node<-test.output.3RWGD.heart.processed.mod$From
nodes.contrast.3RWGD.heart.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.3RWGD.heart.jumps.15spe,test.output.3RWGD.heart.processed.mod,by=c("index","node"))
nodes.contrast.3RWGD.heart.jumps.15spe.removing.jump.node <- nodes.contrast.3RWGD.heart.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.3RWGD.heart.jumps.15spe.removing.jump.node$Exp.abs)),]
nodes.contrast.3RWGD.heart.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.3RWGD.heart.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication","FishWGD"))
nodes.contrast.3RWGD.heart.jumps.15spe.removing.jump.node$pic<-abs(nodes.contrast.3RWGD.heart.jumps.15spe.removing.jump.node$Exp.abs)

########### Testing OC for all nodes of 251 3RWGD exp trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.3RWGD.heart.nojump<-summary.tree(Calibrated.3RWGD.heart.nojump)
summary.3RWGD.heart.nojump$Exp.abs <- abs(summary.3RWGD.heart.nojump$pic.heart)
summary.3RWGD.heart.nojumps<-summary.3RWGD.heart.nojump[,!(names(summary.3RWGD.heart.nojump) %in% c("pic.heart"))]
nodes.contrast.3RWGD.heart.nojumps.15spe <- summary.3RWGD.heart.nojumps[which(!is.na(summary.3RWGD.heart.nojumps$Exp.abs)),]
nodes.contrast.3RWGD.heart.nojumps.15spe$Event <- factor(nodes.contrast.3RWGD.heart.nojumps.15spe$events, levels=c("speciation", "duplication","FishWGD"))
rm(summary.3RWGD.heart.nojump)
rm(summary.3RWGD.heart.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.3RWGD.heart.new.fish<-test.output.3RWGD.heart.processed.mod[which(test.output.3RWGD.heart.processed.mod$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##1636 obs

## Step 6: Analysis on fish specific clades
test.output.3RWGD.heart.processed.new.fish<-unique(test.output.3RWGD.heart.new.fish[c(6,3,7,8,10,2,1)])  #1273obs

########### Proportion of jump events in teleosts in the 556 WGD trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.3RWGD.fish.heart<-NULL
info.jump.3RWGD.fish.heart<-bind_rows(lapply(Calibrated.3RWGD.heart.jump, tree.data.collection.3RWGD.fish))
info.jump.3RWGD.fish.heart$spe.num<-info.jump.3RWGD.fish.heart$internal.events.fish-(info.jump.3RWGD.fish.heart$dup.num + info.jump.3RWGD.fish.heart$WGD.num)
info.jump.3RWGD.fish.heart<-info.jump.3RWGD.fish.heart[c(11,14,3,12,5,4,13)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.3RWGD.fish.heart<-
  test.output.3RWGD.heart.new.fish %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.3RWGD.fish.heart)<-c("index.old","jump.dup")

test.WGD.3RWGD.fish.heart<-
  test.output.3RWGD.heart.new.fish %>%
  filter(node.event=="FishWGD" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.WGD.3RWGD.fish.heart)<-c("index.old","jump.WGD")

test.spe.3RWGD.fish.heart<-
  test.output.3RWGD.heart.new.fish %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.3RWGD.fish.heart)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.3RWGD.heart.fish1<-full_join(test.spe.3RWGD.fish.heart,test.dup.3RWGD.fish.heart)
merged.3RWGD.heart.fish<-full_join(merged.3RWGD.heart.fish1,test.WGD.3RWGD.fish.heart)
info.jump.heart.final.3RWGD.fish<-merge(info.jump.3RWGD.fish.heart,merged.3RWGD.heart.fish, by = c("index.old"))
info.jump.heart.final.3RWGD.fish[is.na(info.jump.heart.final.3RWGD.fish)] <- 0
rm(merged.3RWGD.heart.fish1)

## Have to consider trees where atleast one WGD and one speciation exists
info.jump.heart.final1.3RWGD.fish<-info.jump.heart.final.3RWGD.fish[which(info.jump.heart.final.3RWGD.fish$spe.num>0 & info.jump.heart.final.3RWGD.fish$WGD.num>0),]
rm(info.jump.heart.final.3RWGD.fish)

## Proportion of jumps
info.jump.heart.final1.3RWGD.fish$dup.jump.prob<-(info.jump.heart.final1.3RWGD.fish$jump.dup/(info.jump.heart.final1.3RWGD.fish$dup.num*2))
info.jump.heart.final1.3RWGD.fish$spe.jump.prob<-(info.jump.heart.final1.3RWGD.fish$jump.spe/(info.jump.heart.final1.3RWGD.fish$spe.num*2))
info.jump.heart.final1.3RWGD.fish$WGD.jump.prob<-(info.jump.heart.final1.3RWGD.fish$jump.WGD/(info.jump.heart.final1.3RWGD.fish$WGD.num*2))

## Since,many WGD trees didnot have a duplication node, it produces NA for dup.jump.prob
## To get rid of this
info.jump.heart.final2.3RWGD.fish<-info.jump.heart.final1.3RWGD.fish[complete.cases(info.jump.heart.final1.3RWGD.fish),] #97 obs

## Stats
median(na.omit(info.jump.heart.final2.3RWGD.fish$dup.jump.prob)) #0
median(na.omit(info.jump.heart.final2.3RWGD.fish$spe.jump.prob)) #0.2
median(na.omit(info.jump.heart.final2.3RWGD.fish$WGD.jump.prob)) #0
mean(na.omit(info.jump.heart.final2.3RWGD.fish$dup.jump.prob)) #0.253
mean(na.omit(info.jump.heart.final2.3RWGD.fish$spe.jump.prob)) #0.286
mean(na.omit(info.jump.heart.final2.3RWGD.fish$WGD.jump.prob)) #0.1582

pval.jump.3RWGD.heart.fish.spe2dup <- paired.wilcox(info.jump.heart.final2.3RWGD.fish$spe.jump.prob,info.jump.heart.final2.3RWGD.fish$dup.jump.prob) #p-value =2.79e-1
pval.jump.3RWGD.heart.fish.spe2WGD <- paired.wilcox(info.jump.heart.final2.3RWGD.fish$spe.jump.prob,info.jump.heart.final2.3RWGD.fish$WGD.jump.prob) #p-value = 5.39e-9

## For median data
dtfr.3RWGD.heart.fish<-NULL
dtfr.3RWGD.heart.fish<-data.frame(prop.jump=info.jump.heart.final2.3RWGD.fish$spe.jump.prob,events="speciation")
dtfr.3RWGD.heart.fish<-rbind(dtfr.3RWGD.heart.fish,data.frame(prop.jump=info.jump.heart.final2.3RWGD.fish$dup.jump.prob,events="duplication"))
dtfr.3RWGD.heart.fish<-rbind(dtfr.3RWGD.heart.fish,data.frame(prop.jump=info.jump.heart.final2.3RWGD.fish$WGD.jump.prob,events="FishWGD"))
dtfr.3RWGD.heart.fish$Exp.abs <-abs(dtfr.3RWGD.heart.fish$prop.jump)
dtfr.3RWGD.heart.fish$Event<-factor(dtfr.3RWGD.heart.fish$events, levels=c("speciation","duplication","FishWGD"))


median.data.jump.3RWGD.heart.fish<-NULL
median.data.jump.3RWGD.heart.fish<- dtfr.3RWGD.heart.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.3RWGD.heart.fish$pic.round<-paste0("Median = ",round(median.data.jump.3RWGD.heart.fish$heart.abs,4), sep="")
median.data.jump.3RWGD.heart.fish$count<- paste0(median.data.jump.3RWGD.heart.fish$Event,"\n","(n = ",median.data.jump.3RWGD.heart.fish$freq,")")

plot3R2<-boxplot.new.3RWGD2(dtfr.3RWGD.heart.fish,pval.jump.3RWGD.heart.fish.spe2dup,pval.jump.3RWGD.heart.fish.spe2WGD,"Exp.3RWGD",median.data.jump.3RWGD.heart.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("Jump in exp levels for 5 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.05, xmax = 1.7, ymin = 0.67, ymax =0.67, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.4, y = 0.75, label= pval.jump.3RWGD.heart.fish.spe2dup, fontface = 4)+
  annotate("rect", xmin = 1.05, xmax = 2.8, ymin = 0.88, ymax =0.88, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.7, y = 0.95, label= pval.jump.3RWGD.heart.fish.spe2WGD, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.3RWGD.heart.fish$count) ## check

########### Testing OC for all fish specific nodes of 556 WGD trees supporting jump##################
nodes.contrast.3RWGD.heart.jumps.fish<-nodes.contrast.3RWGD.heart.jumps.15spe[which(nodes.contrast.3RWGD.heart.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 556 WGD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.3RWGD.heart.jumps.fish.removing.jump.node<-nodes.contrast.3RWGD.heart.jumps.15spe.removing.jump.node[which(nodes.contrast.3RWGD.heart.jumps.15spe.removing.jump.node$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

########### Testing OC for all fish specific nodes of 251 3RWGD trees that do not support trait jump for teleosts ##################
nodes.contrast.3RWGD.heart.nojumps.fish<-nodes.contrast.3RWGD.heart.nojumps.15spe[which(nodes.contrast.3RWGD.heart.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")

################## Chunk 5: Analysis on contrasts standardized sure SSD trees for traits ##########################

### Step 1a: Considering (sure.SSD.trees) 4139 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for heart tissue for sure SSD trees 
count<-0
Calibrated.heart.SSD1<-lapply(sure.SSD.trees, diagnostic.plot.test.exp.all.tissue, "TPM.heart")
Calibrated.heart.SSD2<-Calibrated.heart.SSD1[!is.na(Calibrated.heart.SSD1)] 
Calibrated.standardized.SSD.heart<-Calibrated.heart.SSD2[ ! sapply(Calibrated.heart.SSD2, is.null) ]## finally we obtained 2976/4139 trees passing the diagnostic tests
rm(Calibrated.heart.SSD1)
rm(Calibrated.heart.SSD2)

### Step 1b: Considering (Calibrated.standardized.SSD.heart) 2976 trees passing the diagnostic test for heart expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.standardized.SSD.heart,tree.index.collect))
rm(count)

test.output.SSD.heart<-NULL
test.output.SSD.heart<-merge(index.info,test.output.heart1,by="Gene") ## 13748 observations
index.trees.jump<-unique(test.output.SSD.heart$tree.num) ##2135 trees
Calibrated.SSD.heart.jump<-Calibrated.standardized.SSD.heart[c(index.trees.jump)] ##2135/2976=71.74% Brownian trees supported trait jump model for heart expressions
Calibrated.SSD.heart.nojump<-Calibrated.standardized.SSD.heart[-c(index.trees.jump)] ## 841/2976=28.26%  trees with no support for jump in tau
rm(index.info)
rm(index.trees.jump)

### Step2: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.SSD.heart))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.SSD.heart)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 13748 obs
}

test.output.SSD.heart.processed<-bind_cols(test.output.SSD.heart,processed.fout) ## 13748 obs for all the 15 species
test.output.SSD.heart.processed<-test.output.SSD.heart.processed[(which(!(test.output.SSD.heart.processed$clade.name %in% c("Clupeocephala","Osteoglossocephalai")))),] ##13178 obs
rm(i)
rm(test.output.SSD.heart)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

## Step 3: Analysis on all 15 vertebrates clades
test.output.SSD.heart.processed.mod<-unique(test.output.SSD.heart.processed[c(6,3,7,8,10,2,1)])  #10006 obs

########### Proportion of jump events in 15 vertebrates species in the 2135 SSD trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.SSD.heart<-NULL
count<-0
info.jump.SSD.heart<-bind_rows(lapply(Calibrated.SSD.heart.jump, tree.data.collection))
info.jump.SSD.heart$spe.num<-info.jump.SSD.heart$internal.events-info.jump.SSD.heart$dup.num
info.jump.SSD.heart<-info.jump.SSD.heart[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.SSD.heart<-
  test.output.SSD.heart.processed %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.SSD.heart)<-c("index.old","jump.dup")

test.spe.SSD.heart<-
  test.output.SSD.heart.processed %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.SSD.heart)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.SSD.heart<-full_join(test.spe.SSD.heart,test.dup.SSD.heart)
info.jump.heart.final.SSD.all<-merge(info.jump.SSD.heart,merged.SSD.heart, by = c("index.old"))
info.jump.heart.final.SSD.all[is.na(info.jump.heart.final.SSD.all)] <- 0

## Have to consider trees where atleast one duplication and speciation exists
info.jump.heart.final1.SSD.all<-info.jump.heart.final.SSD.all[which(info.jump.heart.final.SSD.all$dup.num>0 & info.jump.heart.final.SSD.all$spe.num>0),]
rm(info.jump.heart.final.SSD.all)

## Proportion of jumps
info.jump.heart.final1.SSD.all$dup.jump.prob<-(info.jump.heart.final1.SSD.all$jump.dup/(info.jump.heart.final1.SSD.all$dup.num*2))
info.jump.heart.final1.SSD.all$spe.jump.prob<-(info.jump.heart.final1.SSD.all$jump.spe/(info.jump.heart.final1.SSD.all$spe.num*2))

## Stats
median(na.omit(info.jump.heart.final1.SSD.all$dup.jump.prob)) #0.1667
median(na.omit(info.jump.heart.final1.SSD.all$spe.jump.prob)) #0.077
mean(na.omit(info.jump.heart.final1.SSD.all$dup.jump.prob)) #0.2955
mean(na.omit(info.jump.heart.final1.SSD.all$spe.jump.prob)) #0.1469


########### Testing OC for all nodes of 2135 SSD trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.SSD.heart.jump<-summary.tree(Calibrated.SSD.heart.jump)
summary.SSD.heart.jump$Exp.abs <- abs(summary.SSD.heart.jump$pic.heart)
summary.SSD.heart.jumps<-summary.SSD.heart.jump[,!(names(summary.SSD.heart.jump) %in% c("pic.heart"))]
nodes.contrast.SSD.heart.jumps.15spe <- summary.SSD.heart.jumps[which(!is.na(summary.SSD.heart.jumps$Exp.abs)),]
nodes.contrast.SSD.heart.jumps.15spe$Event <- factor(nodes.contrast.SSD.heart.jumps.15spe$events, levels=c("speciation", "duplication"))
nodes.contrast.SSD.heart.jumps.15spe<-nodes.contrast.SSD.heart.jumps.15spe[(which(!(nodes.contrast.SSD.heart.jumps.15spe$label %in% c("Clupeocephala","Osteoglossocephalai")))),]
rm(summary.SSD.heart.jump)
rm(summary.SSD.heart.jumps)

############ Testing OC using all 15 species for 2135 SSD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.SSD.heart.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.SSD.heart.jumps.15spe$index<-nodes.contrast.SSD.heart.jumps.15spe$index.tree
test.output.SSD.heart.processed$node<-test.output.SSD.heart.processed$From
nodes.contrast.SSD.heart.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.SSD.heart.jumps.15spe,test.output.SSD.heart.processed,by=c("index","node"))
nodes.contrast.SSD.heart.jumps.15spe.removing.jump.node <- nodes.contrast.SSD.heart.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.SSD.heart.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.SSD.heart.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.SSD.heart.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.SSD.heart.jumps.15spe.removing.jump.node$Exp.abs<-nodes.contrast.SSD.heart.jumps.15spe.removing.jump.node$pic
nodes.contrast.SSD.heart.jumps.15spe.removing.jump.node<-nodes.contrast.SSD.heart.jumps.15spe.removing.jump.node[(which(!(nodes.contrast.SSD.heart.jumps.15spe.removing.jump.node$label %in% c("Clupeocephala","Osteoglossocephalai")))),]

########### Testing OC for all nodes of 841 SSD exp trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.SSD.heart.nojump<-summary.tree(Calibrated.SSD.heart.nojump)
summary.SSD.heart.nojump$Exp.abs <- abs(summary.SSD.heart.nojump$pic)
summary.SSD.heart.nojumps<-summary.SSD.heart.nojump[,!(names(summary.SSD.heart.nojump) %in% c("pic"))]
nodes.contrast.SSD.heart.nojumps.15spe <- summary.SSD.heart.nojumps[which(!is.na(summary.SSD.heart.nojumps$Exp.abs)),]
nodes.contrast.SSD.heart.nojumps.15spe$Event <- factor(nodes.contrast.SSD.heart.nojumps.15spe$events, levels=c("speciation", "duplication"))
nodes.contrast.SSD.heart.nojumps.15spe<-nodes.contrast.SSD.heart.nojumps.15spe[(which(!(nodes.contrast.SSD.heart.nojumps.15spe$label %in% c("Clupeocephala","Osteoglossocephalai")))),] 
rm(summary.SSD.heart.nojump)
rm(summary.SSD.heart.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.SSD.heart.new.fish<-test.output.SSD.heart.processed[which(test.output.SSD.heart.processed$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##3897 obs

## Step 6: Analysis on fish specific clades
test.output.SSD.heart.new.fish.mod<-unique(test.output.SSD.heart.new.fish[c(6,3,7,8,10,2,1)])  #2826 obs


########### Proportion of jump events in teleosts in the 2135 SSD trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.SSD.heart.trees.fish<-NULL
count<-0
info.jump.SSD.heart.trees.fish<-bind_rows(lapply(Calibrated.SSD.heart.jump, tree.data.collection.fish))
info.jump.SSD.heart.trees.fish$spe.num<-info.jump.SSD.heart.trees.fish$internal.events.fish-info.jump.SSD.heart.trees.fish$dup.num
info.jump.SSD.heart.trees.fish<-info.jump.SSD.heart.trees.fish[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.SSD.heart.fish<-
  test.output.SSD.heart.new.fish  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.SSD.heart.fish)<-c("index.old","jump.dup")


test.spe.SSD.heart.fish<-
  test.output.SSD.heart.new.fish %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.SSD.heart.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.SSD.heart.fish<-full_join(test.spe.SSD.heart.fish,test.dup.SSD.heart.fish)
info.jump.trees.SSD.heart.final.fish<-merge(info.jump.SSD.heart.trees.fish,merged.SSD.heart.fish, by = c("index.old"))
info.jump.trees.SSD.heart.final.fish[is.na(info.jump.trees.SSD.heart.final.fish)] <- 0

## Have to consider trees where atleast one fish specific duplication and speciation exists
info.jump.trees.SSD.heart.final1.fish<-info.jump.trees.SSD.heart.final.fish[which(info.jump.trees.SSD.heart.final.fish$dup.num>0 & info.jump.trees.SSD.heart.final.fish$spe.num>0),]

## Proportion of jumps
info.jump.trees.SSD.heart.final1.fish$dup.jump.prob<-(info.jump.trees.SSD.heart.final1.fish$jump.dup/(info.jump.trees.SSD.heart.final1.fish$dup.num*2))
info.jump.trees.SSD.heart.final1.fish$spe.jump.prob<-(info.jump.trees.SSD.heart.final1.fish$jump.spe/(info.jump.trees.SSD.heart.final1.fish$spe.num*2))


median(na.omit(info.jump.trees.SSD.heart.final1.fish$dup.jump.prob)) #0.25
median(na.omit(info.jump.trees.SSD.heart.final1.fish$spe.jump.prob)) #0.167
mean(na.omit(info.jump.trees.SSD.heart.final1.fish$dup.jump.prob)) #0.2989
mean(na.omit(info.jump.trees.SSD.heart.final1.fish$spe.jump.prob)) #0.2161

pval.jump.SSD.heart.fish <- paired.wilcox(info.jump.trees.SSD.heart.final1.fish$spe.jump.prob,info.jump.trees.SSD.heart.final1.fish$dup.jump.prob) #p-value = 6.05e-12


########### Testing OC for all fish specific nodes of 2135 SSD trees supporting jump##################
nodes.contrast.SSD.heart.jumps.fish<-nodes.contrast.SSD.heart.jumps.15spe[which(nodes.contrast.SSD.heart.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] #8522 obs

############ Testing OC for teleosts in 2135 SSD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 

nodes.contrast.SSD.heart.jumps.fish.removing.jump.node <- NULL
nodes.contrast.SSD.heart.jumps.fish$index<-nodes.contrast.SSD.heart.jumps.fish$index.tree
test.output.SSD.heart.new.fish$node<-test.output.SSD.heart.new.fish$From
nodes.contrast.SSD.heart.jumps.fish.removing.jump.node <-anti_join(nodes.contrast.SSD.heart.jumps.fish,test.output.SSD.heart.new.fish,by=c("index","node"))
nodes.contrast.SSD.heart.jumps.fish.removing.jump.node <- nodes.contrast.SSD.heart.jumps.fish.removing.jump.node[which(!is.na(nodes.contrast.SSD.heart.jumps.fish.removing.jump.node$pic)),]
nodes.contrast.SSD.heart.jumps.fish.removing.jump.node$Event <- factor(nodes.contrast.SSD.heart.jumps.fish.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.SSD.heart.jumps.fish.removing.jump.node$Exp.abs<-nodes.contrast.SSD.heart.jumps.fish.removing.jump.node$pic
nodes.contrast.SSD.heart.jumps.fish.removing.jump.node<-nodes.contrast.SSD.heart.jumps.fish.removing.jump.node[(which(!(nodes.contrast.SSD.heart.jumps.fish.removing.jump.node$label %in% c("Clupeocephala","Osteoglossocephalai")))),]

########### Testing OC for all fish specific nodes of 841 SSD exp trees that do not support trait jump for teleosts ##################
nodes.contrast.SSD.heart.nojumps.fish<-nodes.contrast.SSD.heart.nojumps.15spe[which(nodes.contrast.SSD.heart.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] #3120 obs

#save.image("norm_exp_15spe_4tips_new.RData")


######## Chunk6: Testing the ortholog conjecture for all the contrast standardized SSD trees for heart ###########

########### Testing OC for all nodes of 2976 trees supporting jump for 15 species ##################

## Calculation of nodes contrasts for heart
Count<-0
Calibrated.standardized.SSD.heart<-lapply(Calibrated.standardized.SSD.heart,contrast.calc,trait="TPM.heart")

## Summarizing contrasts for normalized exp levels 
summary.heart.SSD<-summary.tree(Calibrated.standardized.SSD.heart)
summary.heart.SSD$Exp.abs <- abs(summary.heart.SSD$pic_TPM.heart)
#summary.heart.all<-summary.heart.all[,!(names(summary.heart.all) %in% c("pic_TPM.heart"))]
nodes.contrast.heart.SSD.15spe <- summary.heart.SSD[which(!is.na(summary.heart.SSD$Exp.abs)),]
nodes.contrast.heart.SSD.15spe$Event <- factor(nodes.contrast.heart.SSD.15spe$events, levels=c("speciation", "duplication"))
rm(summary.heart.SSD)

########### Testing OC for all nodes of 4922 trees supporting jump for 5 teleosts ##################
nodes.contrast.heart.SSD.fish<-nodes.contrast.heart.SSD.15spe[which(nodes.contrast.heart.SSD.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")


