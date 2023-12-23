

### This code is written to work with "trait jump model" for testis

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

test.output.testis<-read.csv(paste0(folder1,"To_pablo/testis/results_testis.csv", sep=""),sep=",",header=T) ##173266 obs

## Only considering trees with posterior probability of jump more than 70%
test.output.testis1<-test.output.testis[which(test.output.testis$JumpProbability>=0.7),] ##35651 obs


########################## Chunk 3: Analysis on all contrasts standardized trees for traits ###############################

### Step 1a: Considering (trees.all.interest) 6923 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for testis tissue for all trees 
count<-0
Calibrated.testis.exp1<-lapply(trees.all.interest, diagnostic.plot.test.exp.all.tissue, "TPM.testis")
Calibrated.testis.exp2<-Calibrated.testis.exp1[!is.na(Calibrated.testis.exp1)] 
Calibrated.all.standardized.testis<-Calibrated.testis.exp2[ ! sapply(Calibrated.testis.exp2, is.null) ]## finally we obtained 3805/6923 tree data passing the diagnostic tests
rm(Calibrated.testis.exp1)
rm(Calibrated.testis.exp2)

### Step 1b: Considering (Calibrated.all.standardized.testis) 3805 tree data passing the diagnostic test for testis expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.all.standardized.testis,tree.index.collect))
rm(count)

test.output.testis.all<-NULL
test.output.testis.all<-merge(index.info,test.output.testis1,by="Gene") ## 15192 observations
index.trees.jump<-unique(test.output.testis.all$tree.num) ##2377 trees
Calibrated.standardized.testis.jump<-Calibrated.all.standardized.testis[c(index.trees.jump)] ##2377/3805=62.47% Brownian trees supported trait jump model for testis expressions 
Calibrated.standardized.testis.nojump<-Calibrated.all.standardized.testis[-c(index.trees.jump)] ## 1428/3805=37.53% trees with no support for jump in testis expression
rm(index.info)
rm(index.trees.jump)

### Step2: Calculating PICs of required trait for trees to match the index number of the trees on which levolution was run
Count<-0
trees.all.interest<-lapply(trees.all.interest,contrast.calc,trait="TPM.testis")
# Count<-0
# Calibrated.standardized.testis.jump<-lapply(Calibrated.standardized.testis.jump,contrast.calc,trait="TPM.testis") 
# rm(Count)

### Step3: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.testis.all))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.testis.all)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 15192 obs
}

test.output.testis.processed<-bind_cols(test.output.testis.all,processed.fout) ## 15192 obs for all the 15 species
rm(i)
rm(test.output.testis.all)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

## Step 4: Analysis on all 15 vertebrates clades
test.output.testis.new.mod<-unique(test.output.testis.processed[c(6,3,7,8,10,2,1)])  #11781 obs

########### Proportion of jump events in 15 vertebrates species in the 2377 trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.trees.testis<-NULL
count<-0
info.jump.trees.testis<-bind_rows(lapply(Calibrated.standardized.testis.jump, tree.data.collection))
info.jump.trees.testis$spe.num<-info.jump.trees.testis$internal.events-info.jump.trees.testis$dup.num
info.jump.trees.testis<-info.jump.trees.testis[c(11,12,3:5)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.testis<-
  test.output.testis.processed  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.testis)<-c("index.old","jump.dup")

test.spe.testis<-
  test.output.testis.processed %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.testis)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.testis<-full_join(test.spe.testis,test.dup.testis)
info.jump.trees.testis.final<-merge(info.jump.trees.testis,merged.testis, by = c("index.old"))
info.jump.trees.testis.final[is.na(info.jump.trees.testis.final)] <- 0

## Have to consider trees where atleast one duplication and speciation exists
info.jump.trees.testis.final1<-info.jump.trees.testis.final[which(info.jump.trees.testis.final$dup.num>0 & info.jump.trees.testis.final$spe.num>0),]
rm(info.jump.trees.testis.final)

## Proportion of jumps
info.jump.trees.testis.final1$dup.jump.prob<-(info.jump.trees.testis.final1$jump.dup/(info.jump.trees.testis.final1$dup.num*2))
info.jump.trees.testis.final1$spe.jump.prob<-(info.jump.trees.testis.final1$jump.spe/(info.jump.trees.testis.final1$spe.num*2))

## Stats
median(na.omit(info.jump.trees.testis.final1$dup.jump.prob)) #0.2428
median(na.omit(info.jump.trees.testis.final1$spe.jump.prob)) #0.1
mean(na.omit(info.jump.trees.testis.final1$dup.jump.prob)) #0.2705
mean(na.omit(info.jump.trees.testis.final1$spe.jump.prob)) #0.1711

pval.jump.testis <- paired.wilcox(info.jump.trees.testis.final1$spe.jump.prob,info.jump.trees.testis.final1$dup.jump.prob) #p-value = 2.2e-16

## For median data
dtfr.testis<-NULL
dtfr.testis<-data.frame(prop.jump=info.jump.trees.testis.final1$spe.jump.prob,events="speciation")
dtfr.testis<-rbind(dtfr.testis,data.frame(prop.jump=info.jump.trees.testis.final1$dup.jump.prob,events="duplication"))
dtfr.testis$Exp.abs <-abs(dtfr.testis$prop.jump)
dtfr.testis$Event<-factor(dtfr.testis$events, levels=c("speciation","duplication"))

library(tidyverse)
median.data.jump.testis<-NULL
median.data.jump.testis<- dtfr.testis%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.testis$pic.round<-paste0("Median = ",round(median.data.jump.testis$Exp.abs,4), sep="")
median.data.jump.testis$count<- paste0(median.data.jump.testis$Event,"\n\n","(n = ",median.data.jump.testis$freq,")")

plot4A<-boxplot.new2(dtfr.testis,pval.jump.testis,"Exp",median.data.jump.testis) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.testis, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.testis$count) ## check

# Create Data
data.jump.testis <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(2377,1428)
)

# Compute the position of labels
data.jump.testis <- data.jump.testis %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.testis$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.testis$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.testis<-ggplot(data.jump.testis, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 2377 trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.testis.jump<-summary.tree(Calibrated.standardized.testis.jump)
summary.testis.jump$Exp.abs <- abs(summary.testis.jump$pic.testis)
summary.testis.jumps<-summary.testis.jump[,!(names(summary.testis.jump) %in% c("pic.testis"))]
nodes.contrast.testis.jumps.15spe <- summary.testis.jumps[which(!is.na(summary.testis.jumps$Exp.abs)),]
nodes.contrast.testis.jumps.15spe$Event <- factor(nodes.contrast.testis.jumps.15spe$events, levels=c("speciation", "duplication"))
rm(summary.testis.jump)
rm(summary.testis.jumps)

############ Testing OC using all 15 species for 2377 trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.testis.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.testis.jumps.15spe$index<-nodes.contrast.testis.jumps.15spe$index.tree
test.output.testis.processed$node<-test.output.testis.processed$From
nodes.contrast.testis.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.testis.jumps.15spe,test.output.testis.processed,by=c("index","node"))
nodes.contrast.testis.jumps.15spe.removing.jump.node <- nodes.contrast.testis.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.testis.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.testis.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.testis.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.testis.jumps.15spe.removing.jump.node$Exp.abs<-nodes.contrast.testis.jumps.15spe.removing.jump.node$pic


########### Testing OC for all nodes of 1428 trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.testis.nojump<-summary.tree(Calibrated.standardized.testis.nojump)
summary.testis.nojump$Exp.abs <- abs(summary.testis.nojump$pic.testis)
summary.testis.nojumps<-summary.testis.nojump[,!(names(summary.testis.nojump) %in% c("pic.testis"))]
nodes.contrast.testis.nojumps.15spe <- summary.testis.nojump[which(!is.na(summary.testis.nojump$Exp.abs)),]
nodes.contrast.testis.nojumps.15spe$Event <- factor(nodes.contrast.testis.nojumps.15spe$events, levels=c("speciation", "duplication"))
rm(summary.testis.nojump)
rm(summary.testis.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.testis.new.fish<-test.output.testis.processed[which(test.output.testis.processed$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##6042 obs

## Step 6: Analysis on fish specific clades
test.output.testis.new.fish.mod<-unique(test.output.testis.new.fish[c(6,3,7,8,10,2,1)])  #4566 obs

########### Proportion of jump events in teleosts in the 2377 trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.trees.testis.fish<-NULL
count<-0
info.jump.trees.testis.fish<-bind_rows(lapply(Calibrated.standardized.testis.jump, tree.data.collection.fish))
info.jump.trees.testis.fish$spe.num<-info.jump.trees.testis.fish$internal.events.fish-info.jump.trees.testis.fish$dup.num
info.jump.trees.testis.fish<-info.jump.trees.testis.fish[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.testis.fish<-
  test.output.testis.new.fish  %>%
  filter(node.event=="duplication") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.testis.fish)<-c("index.old","jump.dup")

test.spe.testis.fish<-
  test.output.testis.new.fish %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.testis.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.testis.fish<-full_join(test.spe.testis.fish,test.dup.testis.fish)
info.jump.trees.testis.final.fish<-merge(info.jump.trees.testis.fish,merged.testis.fish, by = c("index.old"))
info.jump.trees.testis.final.fish[is.na(info.jump.trees.testis.final.fish)] <- 0

## Have to consider trees where atleast one fish specific duplication and speciation exists
info.jump.trees.testis.final1.fish<-info.jump.trees.testis.final.fish[which(info.jump.trees.testis.final.fish$dup.num>0 & info.jump.trees.testis.final.fish$spe.num>0),]

## Proportion of jumps
info.jump.trees.testis.final1.fish$dup.jump.prob<-(info.jump.trees.testis.final1.fish$jump.dup/(info.jump.trees.testis.final1.fish$dup.num*2))
info.jump.trees.testis.final1.fish$spe.jump.prob<-(info.jump.trees.testis.final1.fish$jump.spe/(info.jump.trees.testis.final1.fish$spe.num*2))


median(na.omit(info.jump.trees.testis.final1.fish$dup.jump.prob)) #0.1667
median(na.omit(info.jump.trees.testis.final1.fish$spe.jump.prob)) #0.2042
mean(na.omit(info.jump.trees.testis.final1.fish$dup.jump.prob)) #0.2834
mean(na.omit(info.jump.trees.testis.final1.fish$spe.jump.prob)) #0.2958

pval.jump.testis.fish <- paired.wilcox(info.jump.trees.testis.final1.fish$spe.jump.prob,info.jump.trees.testis.final1.fish$dup.jump.prob) #p-value = 9.87e-02


## For median data
dtfr.testis.fish<-NULL
dtfr.testis.fish<-data.frame(prop.jump=info.jump.trees.testis.final1.fish$spe.jump.prob,events="speciation")
dtfr.testis.fish<-rbind(dtfr.testis.fish,data.frame(prop.jump=info.jump.trees.testis.final1.fish$dup.jump.prob,events="duplication"))
dtfr.testis.fish$Exp.abs <-abs(dtfr.testis.fish$prop.jump)
dtfr.testis.fish$Event<-factor(dtfr.testis.fish$events, levels=c("speciation","duplication"))


library(tidyverse)
median.data.jump.testis.testis.fish<-NULL
median.data.jump.testis.testis.fish<- dtfr.testis.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.testis.testis.fish$pic.round<-paste0("Median = ",round(median.data.jump.testis.testis.fish$Exp.abs,4), sep="")
median.data.jump.testis.testis.fish$count<- paste0(median.data.jump.testis.testis.fish$Event,"\n\n","(n = ",median.data.jump.testis.testis.fish$freq,")")

plot4B<-boxplot.new2(dtfr.testis.fish,pval.jump.testis.fish,"Exp",median.data.jump.testis.testis.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 5 teleosts"))))+
  #ylab(expression(bold(paste("Proportions of jump events for teleosts ")))) +
  ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.testis.fish, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.testis.testis.fish$count) ## check

########### Testing OC for all fish specific nodes of 2377 trees supporting jump##################
nodes.contrast.testis.jumps.fish<-nodes.contrast.testis.jumps.15spe[which(nodes.contrast.testis.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 2377 trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.testis.jumps.fish.removing.jump.node <- NULL
nodes.contrast.testis.jumps.fish$index<-nodes.contrast.testis.jumps.fish$index.tree
test.output.testis.new.fish$node<-test.output.testis.new.fish$From
nodes.contrast.testis.jumps.fish.removing.jump.node <-anti_join(nodes.contrast.testis.jumps.fish,test.output.testis.new.fish,by=c("index","node"))
nodes.contrast.testis.jumps.fish.removing.jump.node <- nodes.contrast.testis.jumps.fish.removing.jump.node[which(!is.na(nodes.contrast.testis.jumps.fish.removing.jump.node$pic)),]
nodes.contrast.testis.jumps.fish.removing.jump.node$Event <- factor(nodes.contrast.testis.jumps.fish.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.testis.jumps.fish.removing.jump.node$Exp.abs<-nodes.contrast.testis.jumps.fish.removing.jump.node$pic

rm(test.output.testis.new.fish)

########### Testing OC for all fish specific nodes of 1428 trees that do not support trait jump for teleosts ##################
nodes.contrast.testis.nojumps.fish<-nodes.contrast.testis.nojumps.15spe[which(nodes.contrast.testis.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")

################## Chunk 4: Analysis on contrasts standardized sure 3RWGD trees for traits ##########################

### Step 1a: Considering (Sure.WGD.trees) 1159 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for testis tissue for sure WGD trees 
count<-0
Calibrated.testis.3Rexp1<-lapply(Sure.WGD.trees, diagnostic.plot.test.exp.all.tissue, "TPM.testis")
Calibrated.testis.3Rexp2<-Calibrated.testis.3Rexp1[!is.na(Calibrated.testis.3Rexp1)] 
Calibrated.3R.standardized.testis<-Calibrated.testis.3Rexp2[ ! sapply(Calibrated.testis.3Rexp2, is.null) ]## finally we obtained 647/1159 trees passing the diagnostic tests
rm(Calibrated.testis.3Rexp1)
rm(Calibrated.testis.3Rexp2)

### Step 1b: Considering (Calibrated.3R.standardized.testis) 647 trees passing the diagnostic test for testis expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.3R.standardized.testis,tree.index.collect))
rm(count)

test.output.3Rexp.testis<-NULL
test.output.3Rexp.testis<-merge(index.info,test.output.testis1,by="Gene") ## 2813 observations
index.trees.jump<-unique(test.output.3Rexp.testis$tree.num) ##408 trees
Calibrated.3RWGD.testis.jump<-Calibrated.3R.standardized.testis[c(index.trees.jump)] ##408/647=63.06% Brownian trees supported trait jump model for testis expressions
Calibrated.3RWGD.testis.nojump<-Calibrated.3R.standardized.testis[-c(index.trees.jump)] ## 239/647=36.94% trees with no support for jump in trait
rm(index.info)
rm(index.trees.jump)

### Step2: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.3Rexp.testis))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.3Rexp.testis)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 2716 obs
}

test.output.3RWGD.testis.processed<-bind_cols(test.output.3Rexp.testis,processed.fout) ## 2813 obs for all the 15 species
rm(i)
rm(test.output.3Rexp.testis)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

##Since our calibrated tree shows both the FishWGD and duplication nodes, 
##we want to modify the node events of the jump output to separate out the effect of FishWGD from the SSDs in those trees.

jump.3RWGD.testis.trees.summary<-summary.tree(Calibrated.3RWGD.testis.jump)
jump.3RWGD.testis.trees.summary$Gene<-jump.3RWGD.testis.trees.summary$index.tree
jump.3RWGD.testis.trees.summary$From<-jump.3RWGD.testis.trees.summary$node
test.output.3RWGD.testis.processed.mod<-merge(test.output.3RWGD.testis.processed,jump.3RWGD.testis.trees.summary,by=c("Gene","From"))
test.output.3RWGD.testis.processed.mod<-test.output.3RWGD.testis.processed.mod[c(1,3,2,4:11,16)]
test.output.3RWGD.testis.processed.mod$node.event<-test.output.3RWGD.testis.processed.mod$events
colnames(test.output.3RWGD.testis.processed.mod)[10]<-"pic"
test.output.3RWGD.testis.processed.mod<-test.output.3RWGD.testis.processed.mod[c(1:11)] ## 2716 obs

## Step 3: Analysis on all 15 vertebrates clades
test.output.3RWGD.testis.processed.new<-unique(test.output.3RWGD.testis.processed.mod[c(6,3,7,8,10,2,1)])  #2139 obs

########### Proportion of jump events in 15 vertebrates species in the 408 3RWGD trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.3RWGD.testis<-NULL
info.jump.3RWGD.testis<-bind_rows(lapply(Calibrated.3RWGD.testis.jump, tree.data.collection.3RWGD))
info.jump.3RWGD.testis$spe.num<-info.jump.3RWGD.testis$internal.events-(info.jump.3RWGD.testis$dup.num + info.jump.3RWGD.testis$WGD.num)
info.jump.3RWGD.testis<-info.jump.3RWGD.testis[c(11,14,3,12,5,4,13)]


## Counting frquency of jump for "speciation" and "duplication" events
test.dup.3RWGD.testis<-
  test.output.3RWGD.testis.processed.mod %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.3RWGD.testis)<-c("index.old","jump.dup")

test.WGD.3RWGD.testis<-
  test.output.3RWGD.testis.processed.mod %>%
  filter(node.event=="FishWGD" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.WGD.3RWGD.testis)<-c("index.old","jump.WGD")

test.spe.3RWGD.testis<-
  test.output.3RWGD.testis.processed.mod %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.3RWGD.testis)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.3RWGD.testis1<-full_join(test.spe.3RWGD.testis,test.dup.3RWGD.testis)
merged.3RWGD.testis<-full_join(merged.3RWGD.testis1,test.WGD.3RWGD.testis)
info.jump.3RWGD.testis.final1<-merge(info.jump.3RWGD.testis,merged.3RWGD.testis, by = c("index.old"))
info.jump.3RWGD.testis.final1[is.na(info.jump.3RWGD.testis.final1)] <- 0
rm(merged.3RWGD.testis1)

## Have to consider trees where atleast one WGD and one speciation exists
info.jump.3RWGD.testis.final<-info.jump.3RWGD.testis.final1[which(info.jump.3RWGD.testis.final1$spe.num>0 & info.jump.3RWGD.testis.final1$WGD.num>0),]
rm(info.jump.3RWGD.testis.final1)

## Proportion of jumps
info.jump.3RWGD.testis.final$dup.jump.prob<-(info.jump.3RWGD.testis.final$jump.dup/(info.jump.3RWGD.testis.final$dup.num*2))
info.jump.3RWGD.testis.final$spe.jump.prob<-(info.jump.3RWGD.testis.final$jump.spe/(info.jump.3RWGD.testis.final$spe.num*2))
info.jump.3RWGD.testis.final$WGD.jump.prob<-(info.jump.3RWGD.testis.final$jump.WGD/(info.jump.3RWGD.testis.final$WGD.num*2))

## Since,many WGD trees didnot have a duplication node, it produces NA for dup.jump.prob
## To get rid of this
info.jump.3RWGD.testis.final2<-info.jump.3RWGD.testis.final[complete.cases(info.jump.3RWGD.testis.final),] ##227 obs

## Stats
median(info.jump.3RWGD.testis.final2$dup.jump.prob) #0.167
median(info.jump.3RWGD.testis.final2$spe.jump.prob) #0.093
median(info.jump.3RWGD.testis.final2$WGD.jump.prob) #0
mean(info.jump.3RWGD.testis.final2$dup.jump.prob) #0.2727
mean(info.jump.3RWGD.testis.final2$spe.jump.prob) #0.1797
mean(info.jump.3RWGD.testis.final2$WGD.jump.prob) #0.1627

pval.jump.3RWGD.testis.spe2dup <- paired.wilcox(info.jump.3RWGD.testis.final2$spe.jump.prob,info.jump.3RWGD.testis.final2$dup.jump.prob) #p-value = 6.65e-4
pval.jump.3RWGD.testis.spe2WGD <- paired.wilcox(info.jump.3RWGD.testis.final2$spe.jump.prob,info.jump.3RWGD.testis.final2$WGD.jump.prob) #p-value = 2.21e-4

## For median data
dtfr.3RWGD.testis.15spe<-NULL
dtfr.3RWGD.testis.15spe<-data.frame(prop.jump=info.jump.3RWGD.testis.final2$spe.jump.prob,events="speciation")
dtfr.3RWGD.testis.15spe<-rbind(dtfr.3RWGD.testis.15spe,data.frame(prop.jump=info.jump.3RWGD.testis.final2$dup.jump.prob,events="duplication"))
dtfr.3RWGD.testis.15spe<-rbind(dtfr.3RWGD.testis.15spe,data.frame(prop.jump=info.jump.3RWGD.testis.final2$WGD.jump.prob,events="FishWGD"))
dtfr.3RWGD.testis.15spe$testis.abs <-abs(dtfr.3RWGD.testis.15spe$prop.jump)
dtfr.3RWGD.testis.15spe$Event<-factor(dtfr.3RWGD.testis.15spe$events, levels=c("speciation","duplication","FishWGD"))


median.data.jump.3RWGD.testis.15spe<-NULL
median.data.jump.3RWGD.testis.15spe<- dtfr.3RWGD.testis.15spe%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.3RWGD.testis.15spe$pic.round<-paste0("Median = ",round(median.data.jump.3RWGD.testis.15spe$Exp.abs,4), sep="")
median.data.jump.3RWGD.testis.15spe$count<- paste0(median.data.jump.3RWGD.testis.15spe$Event,"\n\n","(n = ",median.data.jump.3RWGD.testis.15spe$freq,")")

plot3R1<-boxplot.new.3RWGD2(dtfr.3RWGD.testis.15spe,pval.jump.3RWGD.testis.spe2dup,pval.jump.3RWGD.testis.spe2WGD,"Exp.3RWGD",median.data.jump.3RWGD.testis.15spe) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.05, xmax = 1.7, ymin = 0.67, ymax =0.67, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.4, y = 0.75, label= pval.jump.3RWGD.testis.spe2dup, fontface = 4)+
  annotate("rect", xmin = 1.05, xmax = 2.8, ymin = 0.88, ymax =0.88, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.7, y = 0.95, label= pval.jump.3RWGD.testis.spe2WGD, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.3RWGD.testis.15spe$count) ## check

# Create Data
data.jump.3RWGD.exp <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(445,355)
)

# Compute the position of labels
data.jump.3RWGD.testis <- data.jump.3RWGD.testis %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.3RWGD.testis$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.3RWGD.testis$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.3RWGD.testis<-ggplot(data.jump.3RWGD.testis, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 408 WGD trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.3RWGD.testis.jump<-summary.tree(Calibrated.3RWGD.testis.jump)
summary.3RWGD.testis.jump$Exp.abs <- abs(summary.3RWGD.testis.jump$pic)
summary.3RWGD.testis.jumps<-summary.3RWGD.testis.jump[,!(names(summary.3RWGD.testis.jump) %in% c("pic"))]
nodes.contrast.3RWGD.testis.jumps.15spe <- summary.3RWGD.testis.jumps[which(!is.na(summary.3RWGD.testis.jumps$Exp.abs)),]
nodes.contrast.3RWGD.testis.jumps.15spe$Event <- factor(nodes.contrast.3RWGD.testis.jumps.15spe$events, levels=c("speciation", "duplication", "FishWGD"))
rm(summary.3RWGD.testis.jump)
rm(summary.3RWGD.testis.jumps)

############ Testing OC using all 15 species for 408 WGD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.3RWGD.testis.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.3RWGD.testis.jumps.15spe$index<-nodes.contrast.3RWGD.testis.jumps.15spe$index.tree
test.output.3RWGD.testis.processed.mod$node<-test.output.3RWGD.testis.processed.mod$From
nodes.contrast.3RWGD.testis.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.3RWGD.testis.jumps.15spe,test.output.3RWGD.testis.processed.mod,by=c("index","node"))
nodes.contrast.3RWGD.testis.jumps.15spe.removing.jump.node <- nodes.contrast.3RWGD.testis.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.3RWGD.testis.jumps.15spe.removing.jump.node$Exp.abs)),]
nodes.contrast.3RWGD.testis.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.3RWGD.testis.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication","FishWGD"))
nodes.contrast.3RWGD.testis.jumps.15spe.removing.jump.node$pic<-abs(nodes.contrast.3RWGD.testis.jumps.15spe.removing.jump.node$Exp.abs)

########### Testing OC for all nodes of 239 3RWGD exp trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.3RWGD.testis.nojump<-summary.tree(Calibrated.3RWGD.testis.nojump)
summary.3RWGD.testis.nojump$Exp.abs <- abs(summary.3RWGD.testis.nojump$pic.testis)
summary.3RWGD.testis.nojumps<-summary.3RWGD.testis.nojump[,!(names(summary.3RWGD.testis.nojump) %in% c("pic.testis"))]
nodes.contrast.3RWGD.testis.nojumps.15spe <- summary.3RWGD.testis.nojumps[which(!is.na(summary.3RWGD.testis.nojumps$Exp.abs)),]
nodes.contrast.3RWGD.testis.nojumps.15spe$Event <- factor(nodes.contrast.3RWGD.testis.nojumps.15spe$events, levels=c("speciation", "duplication","FishWGD"))
rm(summary.3RWGD.testis.nojump)
rm(summary.3RWGD.testis.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.3RWGD.testis.new.fish<-test.output.3RWGD.testis.processed.mod[which(test.output.3RWGD.testis.processed.mod$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##1384 obs

## Step 6: Analysis on fish specific clades
test.output.3RWGD.testis.processed.new.fish<-unique(test.output.3RWGD.testis.new.fish[c(6,3,7,8,10,2,1)])  #1022obs

########### Proportion of jump events in teleosts in the 408 WGD trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.3RWGD.fish.testis<-NULL
info.jump.3RWGD.fish.testis<-bind_rows(lapply(Calibrated.3RWGD.testis.jump, tree.data.collection.3RWGD.fish))
info.jump.3RWGD.fish.testis$spe.num<-info.jump.3RWGD.fish.testis$internal.events.fish-(info.jump.3RWGD.fish.testis$dup.num + info.jump.3RWGD.fish.testis$WGD.num)
info.jump.3RWGD.fish.testis<-info.jump.3RWGD.fish.testis[c(11,14,3,12,5,4,13)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.3RWGD.fish.testis<-
  test.output.3RWGD.testis.new.fish %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.3RWGD.fish.testis)<-c("index.old","jump.dup")

test.WGD.3RWGD.fish.testis<-
  test.output.3RWGD.testis.new.fish %>%
  filter(node.event=="FishWGD" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.WGD.3RWGD.fish.testis)<-c("index.old","jump.WGD")

test.spe.3RWGD.fish.testis<-
  test.output.3RWGD.testis.new.fish %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.3RWGD.fish.testis)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.3RWGD.testis.fish1<-full_join(test.spe.3RWGD.fish.testis,test.dup.3RWGD.fish.testis)
merged.3RWGD.testis.fish<-full_join(merged.3RWGD.testis.fish1,test.WGD.3RWGD.fish.testis)
info.jump.testis.final.3RWGD.fish<-merge(info.jump.3RWGD.fish.testis,merged.3RWGD.testis.fish, by = c("index.old"))
info.jump.testis.final.3RWGD.fish[is.na(info.jump.testis.final.3RWGD.fish)] <- 0
rm(merged.3RWGD.testis.fish1)

## Have to consider trees where atleast one WGD and one speciation exists
info.jump.testis.final1.3RWGD.fish<-info.jump.testis.final.3RWGD.fish[which(info.jump.testis.final.3RWGD.fish$spe.num>0 & info.jump.testis.final.3RWGD.fish$WGD.num>0),]
rm(info.jump.testis.final.3RWGD.fish)

## Proportion of jumps
info.jump.testis.final1.3RWGD.fish$dup.jump.prob<-(info.jump.testis.final1.3RWGD.fish$jump.dup/(info.jump.testis.final1.3RWGD.fish$dup.num*2))
info.jump.testis.final1.3RWGD.fish$spe.jump.prob<-(info.jump.testis.final1.3RWGD.fish$jump.spe/(info.jump.testis.final1.3RWGD.fish$spe.num*2))
info.jump.testis.final1.3RWGD.fish$WGD.jump.prob<-(info.jump.testis.final1.3RWGD.fish$jump.WGD/(info.jump.testis.final1.3RWGD.fish$WGD.num*2))

## Since,many WGD trees didnot have a duplication node, it produces NA for dup.jump.prob
## To get rid of this
info.jump.testis.final2.3RWGD.fish<-info.jump.testis.final1.3RWGD.fish[complete.cases(info.jump.testis.final1.3RWGD.fish),] #97 obs

## Stats
median(na.omit(info.jump.testis.final2.3RWGD.fish$dup.jump.prob)) #0.25
median(na.omit(info.jump.testis.final2.3RWGD.fish$spe.jump.prob)) #0.174
median(na.omit(info.jump.testis.final2.3RWGD.fish$WGD.jump.prob)) #0
mean(na.omit(info.jump.testis.final2.3RWGD.fish$dup.jump.prob)) #0.2861
mean(na.omit(info.jump.testis.final2.3RWGD.fish$spe.jump.prob)) #0.289
mean(na.omit(info.jump.testis.final2.3RWGD.fish$WGD.jump.prob)) #0.1742

pval.jump.3RWGD.testis.fish.spe2dup <- paired.wilcox(info.jump.testis.final2.3RWGD.fish$spe.jump.prob,info.jump.testis.final2.3RWGD.fish$dup.jump.prob) #p-value = 7.9e-1
pval.jump.3RWGD.testis.fish.spe2WGD <- paired.wilcox(info.jump.testis.final2.3RWGD.fish$spe.jump.prob,info.jump.testis.final2.3RWGD.fish$WGD.jump.prob) #p-value = 8.61e-4

## For median data
dtfr.3RWGD.testis.fish<-NULL
dtfr.3RWGD.testis.fish<-data.frame(prop.jump=info.jump.testis.final2.3RWGD.fish$spe.jump.prob,events="speciation")
dtfr.3RWGD.testis.fish<-rbind(dtfr.3RWGD.testis.fish,data.frame(prop.jump=info.jump.testis.final2.3RWGD.fish$dup.jump.prob,events="duplication"))
dtfr.3RWGD.testis.fish<-rbind(dtfr.3RWGD.testis.fish,data.frame(prop.jump=info.jump.testis.final2.3RWGD.fish$WGD.jump.prob,events="FishWGD"))
dtfr.3RWGD.testis.fish$Exp.abs <-abs(dtfr.3RWGD.testis.fish$prop.jump)
dtfr.3RWGD.testis.fish$Event<-factor(dtfr.3RWGD.testis.fish$events, levels=c("speciation","duplication","FishWGD"))


median.data.jump.3RWGD.testis.fish<-NULL
median.data.jump.3RWGD.testis.fish<- dtfr.3RWGD.testis.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.3RWGD.testis.fish$pic.round<-paste0("Median = ",round(median.data.jump.3RWGD.testis.fish$testis.abs,4), sep="")
median.data.jump.3RWGD.testis.fish$count<- paste0(median.data.jump.3RWGD.testis.fish$Event,"\n","(n = ",median.data.jump.3RWGD.testis.fish$freq,")")

plot3R2<-boxplot.new.3RWGD2(dtfr.3RWGD.testis.fish,pval.jump.3RWGD.testis.fish.spe2dup,pval.jump.3RWGD.testis.fish.spe2WGD,"Exp.3RWGD",median.data.jump.3RWGD.testis.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("Jump in exp levels for 5 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.05, xmax = 1.7, ymin = 0.67, ymax =0.67, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.4, y = 0.75, label= pval.jump.3RWGD.testis.fish.spe2dup, fontface = 4)+
  annotate("rect", xmin = 1.05, xmax = 2.8, ymin = 0.88, ymax =0.88, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.7, y = 0.95, label= pval.jump.3RWGD.testis.fish.spe2WGD, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.3RWGD.testis.fish$count) ## check

########### Testing OC for all fish specific nodes of 408 WGD trees supporting jump##################
nodes.contrast.3RWGD.testis.jumps.fish<-nodes.contrast.3RWGD.testis.jumps.15spe[which(nodes.contrast.3RWGD.testis.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 408 WGD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.3RWGD.testis.jumps.fish.removing.jump.node<-nodes.contrast.3RWGD.testis.jumps.15spe.removing.jump.node[which(nodes.contrast.3RWGD.testis.jumps.15spe.removing.jump.node$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

########### Testing OC for all fish specific nodes of 239 3RWGD trees that do not support trait jump for teleosts ##################
nodes.contrast.3RWGD.testis.nojumps.fish<-nodes.contrast.3RWGD.testis.nojumps.15spe[which(nodes.contrast.3RWGD.testis.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")

################## Chunk 5: Analysis on contrasts standardized sure SSD trees for traits ##########################

### Step 1a: Considering (sure.SSD.trees) 4139 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for testis tissue for sure SSD trees 
count<-0
Calibrated.testis.SSD1<-lapply(sure.SSD.trees, diagnostic.plot.test.exp.all.tissue, "TPM.testis")
Calibrated.testis.SSD2<-Calibrated.testis.SSD1[!is.na(Calibrated.testis.SSD1)] 
Calibrated.standardized.SSD.testis<-Calibrated.testis.SSD2[ ! sapply(Calibrated.testis.SSD2, is.null) ]## finally we obtained 2636/4139 trees passing the diagnostic tests
rm(Calibrated.testis.SSD1)
rm(Calibrated.testis.SSD2)

### Step 1b: Considering (Calibrated.standardized.SSD.testis) 2636 trees passing the diagnostic test for testis expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.standardized.SSD.testis,tree.index.collect))
rm(count)

test.output.SSD.testis<-NULL
test.output.SSD.testis<-merge(index.info,test.output.testis1,by="Gene") ## 8901 observations
index.trees.jump<-unique(test.output.SSD.testis$tree.num) ##1441 trees
Calibrated.SSD.testis.jump<-Calibrated.standardized.SSD.testis[c(index.trees.jump)] ##1441/2636=54.67% Brownian trees supported trait jump model for testis expressions
Calibrated.SSD.testis.nojump<-Calibrated.standardized.SSD.testis[-c(index.trees.jump)] ## 831/2636=31.52%  trees with no support for jump in tau
rm(index.info)
rm(index.trees.jump)

### Step2: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.SSD.testis))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.SSD.testis)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 8901 obs
}

test.output.SSD.testis.processed<-bind_cols(test.output.SSD.testis,processed.fout) ## 8901 obs for all the 15 species
test.output.SSD.testis.processed<-test.output.SSD.testis.processed[(which(!(test.output.SSD.testis.processed$clade.name %in% c("Clupeocephala","Osteoglossocephalai")))),] ##8484 obs
rm(i)
rm(test.output.SSD.testis)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

## Step 3: Analysis on all 15 vertebrates clades
test.output.SSD.testis.processed.mod<-unique(test.output.SSD.testis.processed[c(6,3,7,8,10,2,1)])  #6504 obs

########### Proportion of jump events in 15 vertebrates species in the 1441 SSD trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.SSD.testis<-NULL
count<-0
info.jump.SSD.testis<-bind_rows(lapply(Calibrated.SSD.testis.jump, tree.data.collection))
info.jump.SSD.testis$spe.num<-info.jump.SSD.testis$internal.events-info.jump.SSD.testis$dup.num
info.jump.SSD.testis<-info.jump.SSD.testis[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.SSD.testis<-
  test.output.SSD.testis.processed %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.SSD.testis)<-c("index.old","jump.dup")

test.spe.SSD.testis<-
  test.output.SSD.testis.processed%>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.SSD.testis)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.SSD.testis<-full_join(test.spe.SSD.testis,test.dup.SSD.testis)
info.jump.testis.final.SSD.all<-merge(info.jump.SSD.testis,merged.SSD.testis, by = c("index.old"))
info.jump.testis.final.SSD.all[is.na(info.jump.testis.final.SSD.all)] <- 0

## Have to consider trees where atleast one duplication and speciation exists
info.jump.testis.final1.SSD.all<-info.jump.testis.final.SSD.all[which(info.jump.testis.final.SSD.all$dup.num>0 & info.jump.testis.final.SSD.all$spe.num>0),]
rm(info.jump.testis.final.SSD.all)

## Proportion of jumps
info.jump.testis.final1.SSD.all$dup.jump.prob<-(info.jump.testis.final1.SSD.all$jump.dup/(info.jump.testis.final1.SSD.all$dup.num*2))
info.jump.testis.final1.SSD.all$spe.jump.prob<-(info.jump.testis.final1.SSD.all$jump.spe/(info.jump.testis.final1.SSD.all$spe.num*2))

## Stats
median(na.omit(info.jump.testis.final1.SSD.all$dup.jump.prob)) #0.25
median(na.omit(info.jump.testis.final1.SSD.all$spe.jump.prob)) #0.083
mean(na.omit(info.jump.testis.final1.SSD.all$dup.jump.prob)) #0.317
mean(na.omit(info.jump.testis.final1.SSD.all$spe.jump.prob)) #0.156

pval.jump.SSD.testis <- paired.wilcox(info.jump.testis.final1.SSD.all$spe.jump.prob,info.jump.testis.final1.SSD.all$dup.jump.prob) #p-value = 2.2e-16

########### Testing OC for all nodes of 1441 SSD trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.SSD.testis.jump<-summary.tree(Calibrated.SSD.testis.jump)
summary.SSD.testis.jump$Exp.abs <- abs(summary.SSD.testis.jump$pic.testis)
summary.SSD.testis.jumps<-summary.SSD.testis.jump[,!(names(summary.SSD.testis.jump) %in% c("pic.testis"))]
nodes.contrast.SSD.testis.jumps.15spe <- summary.SSD.testis.jumps[which(!is.na(summary.SSD.testis.jumps$Exp.abs)),]
nodes.contrast.SSD.testis.jumps.15spe$Event <- factor(nodes.contrast.SSD.testis.jumps.15spe$events, levels=c("speciation", "duplication"))
nodes.contrast.SSD.testis.jumps.15spe<-nodes.contrast.SSD.testis.jumps.15spe[(which(!(nodes.contrast.SSD.testis.jumps.15spe$label %in% c("Clupeocephala","Osteoglossocephalai")))),]
rm(summary.SSD.testis.jump)
rm(summary.SSD.testis.jumps)

############ Testing OC using all 15 species for 1441 SSD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.SSD.testis.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.SSD.testis.jumps.15spe$index<-nodes.contrast.SSD.testis.jumps.15spe$index.tree
test.output.SSD.testis.processed$node<-test.output.SSD.testis.processed$From
nodes.contrast.SSD.testis.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.SSD.testis.jumps.15spe,test.output.SSD.testis.processed,by=c("index","node"))
nodes.contrast.SSD.testis.jumps.15spe.removing.jump.node <- nodes.contrast.SSD.testis.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.SSD.testis.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.SSD.testis.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.SSD.testis.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.SSD.testis.jumps.15spe.removing.jump.node$Exp.abs<-nodes.contrast.SSD.testis.jumps.15spe.removing.jump.node$pic
nodes.contrast.SSD.testis.jumps.15spe.removing.jump.node<-nodes.contrast.SSD.testis.jumps.15spe.removing.jump.node[(which(!(nodes.contrast.SSD.testis.jumps.15spe.removing.jump.node$label %in% c("Clupeocephala","Osteoglossocephalai")))),]

########### Testing OC for all nodes of 831 SSD exp trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.SSD.testis.nojump<-summary.tree(Calibrated.SSD.testis.nojump)
summary.SSD.testis.nojump$Exp.abs <- abs(summary.SSD.testis.nojump$pic)
summary.SSD.testis.nojumps<-summary.SSD.testis.nojump[,!(names(summary.SSD.testis.nojump) %in% c("pic"))]
nodes.contrast.SSD.testis.nojumps.15spe <- summary.SSD.testis.nojumps[which(!is.na(summary.SSD.testis.nojumps$Exp.abs)),]
nodes.contrast.SSD.testis.nojumps.15spe$Event <- factor(nodes.contrast.SSD.testis.nojumps.15spe$events, levels=c("speciation", "duplication"))
nodes.contrast.SSD.testis.nojumps.15spe<-nodes.contrast.SSD.testis.nojumps.15spe[(which(!(nodes.contrast.SSD.testis.nojumps.15spe$label %in% c("Clupeocephala","Osteoglossocephalai")))),] 
rm(summary.SSD.testis.nojump)
rm(summary.SSD.testis.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.SSD.testis.new.fish<-test.output.SSD.testis.processed[which(test.output.SSD.testis.processed$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##7101 obs

## Step 6: Analysis on fish specific clades
test.output.SSD.testis.new.fish.mod<-unique(test.output.SSD.testis.new.fish[c(6,3,7,8,10,2,1)])  #2254 obs


########### Proportion of jump events in teleosts in the 1441 SSD trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.SSD.testis.trees.fish<-NULL
count<-0
info.jump.SSD.testis.trees.fish<-bind_rows(lapply(Calibrated.SSD.testis.jump, tree.data.collection.fish))
info.jump.SSD.testis.trees.fish$spe.num<-info.jump.SSD.testis.trees.fish$internal.events.fish-info.jump.SSD.testis.trees.fish$dup.num
info.jump.SSD.testis.trees.fish<-info.jump.SSD.testis.trees.fish[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.SSD.testis.fish<-
  test.output.SSD.testis.new.fish  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.SSD.testis.fish)<-c("index.old","jump.dup")


test.spe.SSD.testis.fish<-
  test.output.SSD.testis.new.fish %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.SSD.testis.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.SSD.testis.fish<-full_join(test.spe.SSD.testis.fish,test.dup.SSD.testis.fish)
info.jump.trees.SSD.testis.final.fish<-merge(info.jump.SSD.testis.trees.fish,merged.SSD.testis.fish, by = c("index.old"))
info.jump.trees.SSD.testis.final.fish[is.na(info.jump.trees.SSD.testis.final.fish)] <- 0

## Have to consider trees where atleast one fish specific duplication and speciation exists
info.jump.trees.SSD.testis.final1.fish<-info.jump.trees.SSD.testis.final.fish[which(info.jump.trees.SSD.testis.final.fish$dup.num>0 & info.jump.trees.SSD.testis.final.fish$spe.num>0),]

## Proportion of jumps
info.jump.trees.SSD.testis.final1.fish$dup.jump.prob<-(info.jump.trees.SSD.testis.final1.fish$jump.dup/(info.jump.trees.SSD.testis.final1.fish$dup.num*2))
info.jump.trees.SSD.testis.final1.fish$spe.jump.prob<-(info.jump.trees.SSD.testis.final1.fish$jump.spe/(info.jump.trees.SSD.testis.final1.fish$spe.num*2))


median(na.omit(info.jump.trees.SSD.testis.final1.fish$dup.jump.prob)) #0.33
median(na.omit(info.jump.trees.SSD.testis.final1.fish$spe.jump.prob)) #0.1667
mean(na.omit(info.jump.trees.SSD.testis.final1.fish$dup.jump.prob)) #0.357
mean(na.omit(info.jump.trees.SSD.testis.final1.fish$spe.jump.prob)) #0.2239

pval.jump.SSD.testis.fish <- paired.wilcox(info.jump.trees.SSD.testis.final1.fish$spe.jump.prob,info.jump.trees.SSD.testis.final1.fish$dup.jump.prob) #p-value < 2.2e-16

########### Testing OC for all fish specific nodes of 1441 SSD trees supporting jump##################
nodes.contrast.SSD.testis.jumps.fish<-nodes.contrast.SSD.testis.jumps.15spe[which(nodes.contrast.SSD.testis.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 1441 SSD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 

nodes.contrast.SSD.testis.jumps.fish.removing.jump.node <- NULL
nodes.contrast.SSD.testis.jumps.fish$index<-nodes.contrast.SSD.testis.jumps.fish$index.tree
test.output.SSD.testis.new.fish$node<-test.output.SSD.testis.new.fish$From
nodes.contrast.SSD.testis.jumps.fish.removing.jump.node <-anti_join(nodes.contrast.SSD.testis.jumps.fish,test.output.SSD.testis.new.fish,by=c("index","node"))
nodes.contrast.SSD.testis.jumps.fish.removing.jump.node <- nodes.contrast.SSD.testis.jumps.fish.removing.jump.node[which(!is.na(nodes.contrast.SSD.testis.jumps.fish.removing.jump.node$pic)),]
nodes.contrast.SSD.testis.jumps.fish.removing.jump.node$Event <- factor(nodes.contrast.SSD.testis.jumps.fish.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.SSD.testis.jumps.fish.removing.jump.node$Exp.abs<-nodes.contrast.SSD.testis.jumps.fish.removing.jump.node$pic
nodes.contrast.SSD.testis.jumps.fish.removing.jump.node<-nodes.contrast.SSD.testis.jumps.fish.removing.jump.node[(which(!(nodes.contrast.SSD.testis.jumps.fish.removing.jump.node$label %in% c("Clupeocephala","Osteoglossocephalai")))),]

########### Testing OC for all fish specific nodes of 831 SSD exp trees that do not support trait jump for teleosts ##################
nodes.contrast.SSD.testis.nojumps.fish<-nodes.contrast.SSD.testis.nojumps.15spe[which(nodes.contrast.SSD.testis.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")


