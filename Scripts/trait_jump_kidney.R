

### This code is written to work with "trait jump model" for kidney

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

test.output.kidney<-read.csv(paste0(folder1,"To_pablo/kidney/results_kidney.csv", sep=""),sep=",",header=T) ##174178 obs

## Only considering trees with posterior probability of jump more than 70%
test.output.kidney1<-test.output.kidney[which(test.output.kidney$JumpProbability>=0.7),] ##34054 obs

########################## Chunk 3: Analysis on all contrasts standardized trees for traits ###############################
### Step 1a: Considering (trees.all.interest) 6923 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for kidney tissue for all trees 
count<-0
Calibrated.kidney.exp1<-lapply(trees.all.interest, diagnostic.plot.test.exp.all.tissue, "TPM.kidney")
Calibrated.kidney.exp2<-Calibrated.kidney.exp1[!is.na(Calibrated.kidney.exp1)] 
Calibrated.all.standardized.kidney<-Calibrated.kidney.exp2[ ! sapply(Calibrated.kidney.exp2, is.null) ]## finally we obtained 4622/6923 tree data passing the diagnostic tests
rm(Calibrated.kidney.exp1)
rm(Calibrated.kidney.exp2)

### Step 1b: Considering (Calibrated.all.standardized.kidney) 4622 tree data passing the diagnostic test for kidney expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.all.standardized.kidney,tree.index.collect))
rm(count)

test.output.kidney.all<-NULL
test.output.kidney.all<-merge(index.info,test.output.kidney1,by="Gene") ## 19530 observations
index.trees.jump<-unique(test.output.kidney.all$tree.num) ##2987 trees
Calibrated.standardized.kidney.jump<-Calibrated.all.standardized.kidney[c(index.trees.jump)] ##2987/4622=64.62% Brownian trees supported trait jump model for kidney expressions 
Calibrated.standardized.kidney.nojump<-Calibrated.all.standardized.kidney[-c(index.trees.jump)] ## 1635/4622=35.37% trees with no support for jump in kidney expression
rm(index.info)
rm(index.trees.jump)

### Step2: Calculating PICs of required trait for trees to match the index number of the trees on which levolution was run
Count<-0
trees.all.interest<-lapply(trees.all.interest,contrast.calc,trait="TPM.kidney")
# Count<-0
# Calibrated.standardized.kidney.jump<-lapply(Calibrated.standardized.kidney.jump,contrast.calc,trait="TPM.kidney") 
# rm(Count)

### Step3: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.kidney.all))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.kidney.all)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 19530 obs
}
test.output.kidney.processed<-bind_cols(test.output.kidney.all,processed.fout) ## 19530 obs for all the 15 species
rm(i)
rm(test.output.kidney.all)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

## Step 4: Analysis on all 15 vertebrates clades
test.output.kidney.new.mod<-unique(test.output.kidney.processed[c(6,3,7,8,10,2,1)])  #14815 obs

########### Proportion of jump events in 15 vertebrates species in the 2987 trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.trees.kidney<-NULL
count<-0
info.jump.trees.kidney<-bind_rows(lapply(Calibrated.standardized.kidney.jump, tree.data.collection))
info.jump.trees.kidney$spe.num<-info.jump.trees.kidney$internal.events-info.jump.trees.kidney$dup.num
info.jump.trees.kidney<-info.jump.trees.kidney[c(11,12,3:5)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.kidney<-
  test.output.kidney.processed  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.kidney)<-c("index.old","jump.dup")

test.spe.kidney<-
  test.output.kidney.processed %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.kidney)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.kidney<-full_join(test.spe.kidney,test.dup.kidney)
info.jump.trees.kidney.final<-merge(info.jump.trees.kidney,merged.kidney, by = c("index.old"))
info.jump.trees.kidney.final[is.na(info.jump.trees.kidney.final)] <- 0

## Have to consider trees where atleast one duplication and speciation exists
info.jump.trees.kidney.final1<-info.jump.trees.kidney.final[which(info.jump.trees.kidney.final$dup.num>0 & info.jump.trees.kidney.final$spe.num>0),]
rm(info.jump.trees.kidney.final)

## Proportion of jumps
info.jump.trees.kidney.final1$dup.jump.prob<-(info.jump.trees.kidney.final1$jump.dup/(info.jump.trees.kidney.final1$dup.num*2))
info.jump.trees.kidney.final1$spe.jump.prob<-(info.jump.trees.kidney.final1$jump.spe/(info.jump.trees.kidney.final1$spe.num*2))

## Stats
median(na.omit(info.jump.trees.kidney.final1$dup.jump.prob)) #1667
median(na.omit(info.jump.trees.kidney.final1$spe.jump.prob)) #0.0769
mean(na.omit(info.jump.trees.kidney.final1$dup.jump.prob)) #0.2897
mean(na.omit(info.jump.trees.kidney.final1$spe.jump.prob)) #0.1532

pval.jump.kidney <- paired.wilcox(info.jump.trees.kidney.final1$spe.jump.prob,info.jump.trees.kidney.final1$dup.jump.prob) #p-value = 2.2e-16

## For median data
dtfr.kidney<-NULL
dtfr.kidney<-data.frame(prop.jump=info.jump.trees.kidney.final1$spe.jump.prob,events="speciation")
dtfr.kidney<-rbind(dtfr.kidney,data.frame(prop.jump=info.jump.trees.kidney.final1$dup.jump.prob,events="duplication"))
dtfr.kidney$Exp.abs <-abs(dtfr.kidney$prop.jump)
dtfr.kidney$Event<-factor(dtfr.kidney$events, levels=c("speciation","duplication"))

library(tidyverse)
median.data.jump.kidney<-NULL
median.data.jump.kidney<- dtfr.kidney%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.kidney$pic.round<-paste0("Median = ",round(median.data.jump.kidney$Exp.abs,4), sep="")
median.data.jump.kidney$count<- paste0(median.data.jump.kidney$Event,"\n\n","(n = ",median.data.jump.kidney$freq,")")

plot4A<-boxplot.new2(dtfr.kidney,pval.jump.kidney,"Exp",median.data.jump.kidney) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.kidney, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.kidney$count) ## check

# Create Data
data.jump.kidney <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(2987,1635)
)

# Compute the position of labels
data.jump.kidney <- data.jump.kidney %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.kidney$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.kidney$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.kidney<-ggplot(data.jump.kidney, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 2987 trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.kidney.jump<-summary.tree(Calibrated.standardized.kidney.jump)
summary.kidney.jump$Exp.abs <- abs(summary.kidney.jump$pic.kidney)
summary.kidney.jumps<-summary.kidney.jump[,!(names(summary.kidney.jump) %in% c("pic.kidney"))]
nodes.contrast.kidney.jumps.15spe <- summary.kidney.jumps[which(!is.na(summary.kidney.jumps$Exp.abs)),]
nodes.contrast.kidney.jumps.15spe$Event <- factor(nodes.contrast.kidney.jumps.15spe$events, levels=c("speciation", "duplication"))
rm(summary.kidney.jump)
rm(summary.kidney.jumps)

############ Testing OC using all 15 species for 2987 trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.kidney.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.kidney.jumps.15spe$index<-nodes.contrast.kidney.jumps.15spe$index.tree
test.output.kidney.processed$node<-test.output.kidney.processed$From
nodes.contrast.kidney.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.kidney.jumps.15spe,test.output.kidney.processed,by=c("index","node"))
nodes.contrast.kidney.jumps.15spe.removing.jump.node <- nodes.contrast.kidney.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.kidney.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.kidney.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.kidney.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.kidney.jumps.15spe.removing.jump.node$Exp.abs<-nodes.contrast.kidney.jumps.15spe.removing.jump.node$pic


########### Testing OC for all nodes of 1635 trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.kidney.nojump<-summary.tree(Calibrated.standardized.kidney.nojump)
summary.kidney.nojump$Exp.abs <- abs(summary.kidney.nojump$pic.kidney)
summary.kidney.nojumps<-summary.kidney.nojump[,!(names(summary.kidney.nojump) %in% c("pic.kidney"))]
nodes.contrast.kidney.nojumps.15spe <- summary.kidney.nojump[which(!is.na(summary.kidney.nojump$Exp.abs)),]
nodes.contrast.kidney.nojumps.15spe$Event <- factor(nodes.contrast.kidney.nojumps.15spe$events, levels=c("speciation", "duplication"))
rm(summary.kidney.nojump)
rm(summary.kidney.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.kidney.new.fish<-test.output.kidney.processed[which(test.output.kidney.processed$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##8677 obs

## Step 6: Analysis on fish specific clades
test.output.kidney.new.fish.mod<-unique(test.output.kidney.new.fish[c(6,3,7,8,10,2,1)])  #7581 obs

########### Proportion of jump events in teleosts in the 2987 trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.trees.kidney.fish<-NULL
count<-0
info.jump.trees.kidney.fish<-bind_rows(lapply(Calibrated.standardized.kidney.jump, tree.data.collection.fish))
info.jump.trees.kidney.fish$spe.num<-info.jump.trees.kidney.fish$internal.events.fish-info.jump.trees.kidney.fish$dup.num
info.jump.trees.kidney.fish<-info.jump.trees.kidney.fish[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.kidney.fish<-
  test.output.kidney.new.fish  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.kidney.fish)<-c("index.old","jump.dup")

test.spe.kidney.fish<-
  test.output.kidney.new.fish%>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.kidney.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.kidney.fish<-full_join(test.spe.kidney.fish,test.dup.kidney.fish)
info.jump.trees.kidney.final.fish<-merge(info.jump.trees.kidney.fish,merged.kidney.fish, by = c("index.old"))
info.jump.trees.kidney.final.fish[is.na(info.jump.trees.kidney.final.fish)] <- 0

## Have to consider trees where atleast one fish specific duplication and speciation exists
info.jump.trees.kidney.final1.fish<-info.jump.trees.kidney.final.fish[which(info.jump.trees.kidney.final.fish$dup.num>0 & info.jump.trees.kidney.final.fish$spe.num>0),]

## Proportion of jumps
info.jump.trees.kidney.final1.fish$dup.jump.prob<-(info.jump.trees.kidney.final1.fish$jump.dup/(info.jump.trees.kidney.final1.fish$dup.num*2))
info.jump.trees.kidney.final1.fish$spe.jump.prob<-(info.jump.trees.kidney.final1.fish$jump.spe/(info.jump.trees.kidney.final1.fish$spe.num*2))


median(na.omit(info.jump.trees.kidney.final1.fish$dup.jump.prob)) #0.125
median(na.omit(info.jump.trees.kidney.final1.fish$spe.jump.prob)) #0.1667
mean(na.omit(info.jump.trees.kidney.final1.fish$dup.jump.prob)) #0.2534
mean(na.omit(info.jump.trees.kidney.final1.fish$spe.jump.prob)) #0.2411

pval.jump.kidney.fish <- paired.wilcox(info.jump.trees.kidney.final1.fish$spe.jump.prob,info.jump.trees.kidney.final1.fish$dup.jump.prob) #p-value = 2.32e-1


## For median data
dtfr.kidney.fish<-NULL
dtfr.kidney.fish<-data.frame(prop.jump=info.jump.trees.kidney.final1.fish$spe.jump.prob,events="speciation")
dtfr.kidney.fish<-rbind(dtfr.kidney.fish,data.frame(prop.jump=info.jump.trees.kidney.final1.fish$dup.jump.prob,events="duplication"))
dtfr.kidney.fish$Exp.abs <-abs(dtfr.kidney.fish$prop.jump)
dtfr.kidney.fish$Event<-factor(dtfr.kidney.fish$events, levels=c("speciation","duplication"))


library(tidyverse)
median.data.jump.kidney.kidney.fish<-NULL
median.data.jump.kidney.kidney.fish<- dtfr.kidney.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.kidney.kidney.fish$pic.round<-paste0("Median = ",round(median.data.jump.kidney.kidney.fish$Exp.abs,4), sep="")
median.data.jump.kidney.kidney.fish$count<- paste0(median.data.jump.kidney.kidney.fish$Event,"\n\n","(n = ",median.data.jump.kidney.kidney.fish$freq,")")

plot4B<-boxplot.new2(dtfr.kidney.fish,pval.jump.kidney.fish,"Exp",median.data.jump.kidney.kidney.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 5 teleosts"))))+
  #ylab(expression(bold(paste("Proportions of jump events for teleosts ")))) +
  ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.kidney.fish, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.kidney.kidney.fish$count) ## check

########### Testing OC for all fish specific nodes of 2987 trees supporting jump##################
nodes.contrast.kidney.jumps.fish<-nodes.contrast.kidney.jumps.15spe[which(nodes.contrast.kidney.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 2987 trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.kidney.jumps.fish.removing.jump.node <- NULL
nodes.contrast.kidney.jumps.fish$index<-nodes.contrast.kidney.jumps.fish$index.tree
test.output.kidney.new.fish$node<-test.output.kidney.new.fish$From
nodes.contrast.kidney.jumps.fish.removing.jump.node <-anti_join(nodes.contrast.kidney.jumps.fish,test.output.kidney.new.fish,by=c("index","node"))
nodes.contrast.kidney.jumps.fish.removing.jump.node <- nodes.contrast.kidney.jumps.fish.removing.jump.node[which(!is.na(nodes.contrast.kidney.jumps.fish.removing.jump.node$pic)),]
nodes.contrast.kidney.jumps.fish.removing.jump.node$Event <- factor(nodes.contrast.kidney.jumps.fish.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.kidney.jumps.fish.removing.jump.node$Exp.abs<-nodes.contrast.kidney.jumps.fish.removing.jump.node$pic

rm(test.output.kidney.new.fish)

########### Testing OC for all fish specific nodes of 1635 trees that do not support trait jump for teleosts ##################
nodes.contrast.kidney.nojumps.fish<-nodes.contrast.kidney.nojumps.15spe[which(nodes.contrast.kidney.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")

################## Chunk 4: Analysis on contrasts standardized sure 3RWGD trees for traits ##########################

### Step 1a: Considering (Sure.WGD.trees) 1159 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for kidney tissue for sure WGD trees 
count<-0
Calibrated.kidney.3Rexp1<-lapply(Sure.WGD.trees, diagnostic.plot.test.exp.all.tissue, "TPM.kidney")
Calibrated.kidney.3Rexp2<-Calibrated.kidney.3Rexp1[!is.na(Calibrated.kidney.3Rexp1)] 
Calibrated.3R.standardized.kidney<-Calibrated.kidney.3Rexp2[ ! sapply(Calibrated.kidney.3Rexp2, is.null) ]## finally we obtained 818/1159 trees passing the diagnostic tests
rm(Calibrated.kidney.3Rexp1)
rm(Calibrated.kidney.3Rexp2)

### Step 1b: Considering (Calibrated.3R.standardized.kidney) 892 trees passing the diagnostic test for kidney expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.3R.standardized.kidney,tree.index.collect))
rm(count)

test.output.3Rexp.kidney<-NULL
test.output.3Rexp.kidney<-merge(index.info,test.output.kidney1,by="Gene") ## 3185 observations
index.trees.jump<-unique(test.output.3Rexp.kidney$tree.num) ##495 trees
Calibrated.3RWGD.kidney.jump<-Calibrated.3R.standardized.kidney[c(index.trees.jump)] ##495/818=60.51% Brownian trees supported trait jump model for kidney expressions
Calibrated.3RWGD.kidney.nojump<-Calibrated.3R.standardized.kidney[-c(index.trees.jump)] ## 323/818=39.48% trees with no support for jump in trait
rm(index.info)
rm(index.trees.jump)

### Step2: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.3Rexp.kidney))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.3Rexp.kidney)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 3185 obs
}

test.output.3RWGD.kidney.processed<-bind_cols(test.output.3Rexp.kidney,processed.fout) ## 3185 obs for all the 15 species
rm(i)
rm(test.output.3Rexp.kidney)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

##Since our calibrated tree shows both the FishWGD and duplication nodes, 
##we want to modify the node events of the jump output to separate out the effect of FishWGD from the SSDs in those trees.

jump.3RWGD.kidney.trees.summary<-summary.tree(Calibrated.3RWGD.kidney.jump)
jump.3RWGD.kidney.trees.summary$Gene<-jump.3RWGD.kidney.trees.summary$index.tree
jump.3RWGD.kidney.trees.summary$From<-jump.3RWGD.kidney.trees.summary$node
test.output.3RWGD.kidney.processed.mod<-merge(test.output.3RWGD.kidney.processed,jump.3RWGD.kidney.trees.summary,by=c("Gene","From"))
test.output.3RWGD.kidney.processed.mod<-test.output.3RWGD.kidney.processed.mod[c(1,3,2,4:11,16)]
test.output.3RWGD.kidney.processed.mod$node.event<-test.output.3RWGD.kidney.processed.mod$events
colnames(test.output.3RWGD.kidney.processed.mod)[10]<-"pic"
test.output.3RWGD.kidney.processed.mod<-test.output.3RWGD.kidney.processed.mod[c(1:11)] ## 3185 obs

## Step 3: Analysis on all 15 vertebrates clades
test.output.3RWGD.kidney.processed.new<-unique(test.output.3RWGD.kidney.processed.mod[c(6,3,7,8,10,2,1)])  #2513 obs

########### Proportion of jump events in 15 vertebrates species in the 495 3RWGD trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.3RWGD.kidney<-NULL
info.jump.3RWGD.kidney<-bind_rows(lapply(Calibrated.3RWGD.kidney.jump, tree.data.collection.3RWGD))
info.jump.3RWGD.kidney$spe.num<-info.jump.3RWGD.kidney$internal.events-(info.jump.3RWGD.kidney$dup.num + info.jump.3RWGD.kidney$WGD.num)
info.jump.3RWGD.kidney<-info.jump.3RWGD.kidney[c(11,14,3,12,5,4,13)]


## Counting frquency of jump for "speciation" and "duplication" events
test.dup.3RWGD.kidney<-
  test.output.3RWGD.kidney.processed.mod %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.3RWGD.kidney)<-c("index.old","jump.dup")

test.WGD.3RWGD.kidney<-
  test.output.3RWGD.kidney.processed.mod %>%
  filter(node.event=="FishWGD" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.WGD.3RWGD.kidney)<-c("index.old","jump.WGD")

test.spe.3RWGD.kidney<-
  test.output.3RWGD.kidney.processed.mod %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.3RWGD.kidney)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.3RWGD.kidney1<-full_join(test.spe.3RWGD.kidney,test.dup.3RWGD.kidney)
merged.3RWGD.kidney<-full_join(merged.3RWGD.kidney1,test.WGD.3RWGD.kidney)
info.jump.3RWGD.kidney.final1<-merge(info.jump.3RWGD.kidney,merged.3RWGD.kidney, by = c("index.old"))
info.jump.3RWGD.kidney.final1[is.na(info.jump.3RWGD.kidney.final1)] <- 0
rm(merged.3RWGD.kidney1)

## Have to consider trees where atleast one WGD and one speciation exists
info.jump.3RWGD.kidney.final<-info.jump.3RWGD.kidney.final1[which(info.jump.3RWGD.kidney.final1$spe.num>0 & info.jump.3RWGD.kidney.final1$WGD.num>0),]
rm(info.jump.3RWGD.kidney.final1)

## Proportion of jumps
info.jump.3RWGD.kidney.final$dup.jump.prob<-(info.jump.3RWGD.kidney.final$jump.dup/(info.jump.3RWGD.kidney.final$dup.num*2))
info.jump.3RWGD.kidney.final$spe.jump.prob<-(info.jump.3RWGD.kidney.final$jump.spe/(info.jump.3RWGD.kidney.final$spe.num*2))
info.jump.3RWGD.kidney.final$WGD.jump.prob<-(info.jump.3RWGD.kidney.final$jump.WGD/(info.jump.3RWGD.kidney.final$WGD.num*2))

## Since,many WGD trees didnot have a duplication node, it produces NA for dup.jump.prob
## To get rid of this
info.jump.3RWGD.kidney.final2<-info.jump.3RWGD.kidney.final[complete.cases(info.jump.3RWGD.kidney.final),] ##227 obs

## Stats
median(info.jump.3RWGD.kidney.final2$dup.jump.prob) #0.25
median(info.jump.3RWGD.kidney.final2$spe.jump.prob) #0.071
median(info.jump.3RWGD.kidney.final2$WGD.jump.prob) #0
mean(info.jump.3RWGD.kidney.final2$dup.jump.prob) #03453
mean(info.jump.3RWGD.kidney.final2$spe.jump.prob) #0.1591
mean(info.jump.3RWGD.kidney.final2$WGD.jump.prob) #0.1479

pval.jump.3RWGD.kidney.spe2dup <- paired.wilcox(info.jump.3RWGD.kidney.final2$spe.jump.prob,info.jump.3RWGD.kidney.final2$dup.jump.prob) #p-value = 2.2e-16
pval.jump.3RWGD.kidney.spe2WGD <- paired.wilcox(info.jump.3RWGD.kidney.final2$spe.jump.prob,info.jump.3RWGD.kidney.final2$WGD.jump.prob) #p-value = 1.33e-6

## For median data
dtfr.3RWGD.kidney.15spe<-NULL
dtfr.3RWGD.kidney.15spe<-data.frame(prop.jump=info.jump.3RWGD.kidney.final2$spe.jump.prob,events="speciation")
dtfr.3RWGD.kidney.15spe<-rbind(dtfr.3RWGD.kidney.15spe,data.frame(prop.jump=info.jump.3RWGD.kidney.final2$dup.jump.prob,events="duplication"))
dtfr.3RWGD.kidney.15spe<-rbind(dtfr.3RWGD.kidney.15spe,data.frame(prop.jump=info.jump.3RWGD.kidney.final2$WGD.jump.prob,events="FishWGD"))
dtfr.3RWGD.kidney.15spe$kidney.abs <-abs(dtfr.3RWGD.kidney.15spe$prop.jump)
dtfr.3RWGD.kidney.15spe$Event<-factor(dtfr.3RWGD.kidney.15spe$events, levels=c("speciation","duplication","FishWGD"))


median.data.jump.3RWGD.kidney.15spe<-NULL
median.data.jump.3RWGD.kidney.15spe<- dtfr.3RWGD.kidney.15spe%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.3RWGD.kidney.15spe$pic.round<-paste0("Median = ",round(median.data.jump.3RWGD.kidney.15spe$Exp.abs,4), sep="")
median.data.jump.3RWGD.kidney.15spe$count<- paste0(median.data.jump.3RWGD.kidney.15spe$Event,"\n\n","(n = ",median.data.jump.3RWGD.kidney.15spe$freq,")")

plot3R1<-boxplot.new.3RWGD2(dtfr.3RWGD.kidney.15spe,pval.jump.3RWGD.kidney.spe2dup,pval.jump.3RWGD.kidney.spe2WGD,"Exp.3RWGD",median.data.jump.3RWGD.kidney.15spe) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.05, xmax = 1.7, ymin = 0.67, ymax =0.67, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.4, y = 0.75, label= pval.jump.3RWGD.kidney.spe2dup, fontface = 4)+
  annotate("rect", xmin = 1.05, xmax = 2.8, ymin = 0.88, ymax =0.88, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.7, y = 0.95, label= pval.jump.3RWGD.kidney.spe2WGD, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.3RWGD.kidney.15spe$count) ## check

# Create Data
data.jump.3RWGD.exp <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(495,323)
)

# Compute the position of labels
data.jump.3RWGD.kidney <- data.jump.3RWGD.kidney %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.3RWGD.kidney$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.3RWGD.kidney$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.3RWGD.kidney<-ggplot(data.jump.3RWGD.kidney, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 495 WGD trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.3RWGD.kidney.jump<-summary.tree(Calibrated.3RWGD.kidney.jump)
summary.3RWGD.kidney.jump$Exp.abs <- abs(summary.3RWGD.kidney.jump$pic)
summary.3RWGD.kidney.jumps<-summary.3RWGD.kidney.jump[,!(names(summary.3RWGD.kidney.jump) %in% c("pic"))]
nodes.contrast.3RWGD.kidney.jumps.15spe <- summary.3RWGD.kidney.jumps[which(!is.na(summary.3RWGD.kidney.jumps$Exp.abs)),]
nodes.contrast.3RWGD.kidney.jumps.15spe$Event <- factor(nodes.contrast.3RWGD.kidney.jumps.15spe$events, levels=c("speciation", "duplication", "FishWGD"))
rm(summary.3RWGD.kidney.jump)
rm(summary.3RWGD.kidney.jumps)

############ Testing OC using all 15 species for 495 WGD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.3RWGD.kidney.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.3RWGD.kidney.jumps.15spe$index<-nodes.contrast.3RWGD.kidney.jumps.15spe$index.tree
test.output.3RWGD.kidney.processed.mod$node<-test.output.3RWGD.kidney.processed.mod$From
nodes.contrast.3RWGD.kidney.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.3RWGD.kidney.jumps.15spe,test.output.3RWGD.kidney.processed.mod,by=c("index","node"))
nodes.contrast.3RWGD.kidney.jumps.15spe.removing.jump.node <- nodes.contrast.3RWGD.kidney.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.3RWGD.kidney.jumps.15spe.removing.jump.node$Exp.abs)),]
nodes.contrast.3RWGD.kidney.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.3RWGD.kidney.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication","FishWGD"))
nodes.contrast.3RWGD.kidney.jumps.15spe.removing.jump.node$pic<-abs(nodes.contrast.3RWGD.kidney.jumps.15spe.removing.jump.node$Exp.abs)

########### Testing OC for all nodes of 495 3RWGD exp trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.3RWGD.kidney.nojump<-summary.tree(Calibrated.3RWGD.kidney.nojump)
summary.3RWGD.kidney.nojump$Exp.abs <- abs(summary.3RWGD.kidney.nojump$pic.kidney)
summary.3RWGD.kidney.nojumps<-summary.3RWGD.kidney.nojump[,!(names(summary.3RWGD.kidney.nojump) %in% c("pic.kidney"))]
nodes.contrast.3RWGD.kidney.nojumps.15spe <- summary.3RWGD.kidney.nojumps[which(!is.na(summary.3RWGD.kidney.nojumps$Exp.abs)),]
nodes.contrast.3RWGD.kidney.nojumps.15spe$Event <- factor(nodes.contrast.3RWGD.kidney.nojumps.15spe$events, levels=c("speciation", "duplication","FishWGD"))
rm(summary.3RWGD.kidney.nojump)
rm(summary.3RWGD.kidney.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.3RWGD.kidney.new.fish<-test.output.3RWGD.kidney.processed.mod[which(test.output.3RWGD.kidney.processed.mod$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##1520 obs

## Step 6: Analysis on fish specific clades
test.output.3RWGD.kidney.processed.new.fish<-unique(test.output.3RWGD.kidney.new.fish[c(6,3,7,8,10,2,1)])  #1180obs

########### Proportion of jump events in teleosts in the 495 WGD trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.3RWGD.fish.kidney<-NULL
info.jump.3RWGD.fish.kidney<-bind_rows(lapply(Calibrated.3RWGD.kidney.jump, tree.data.collection.3RWGD.fish))
info.jump.3RWGD.fish.kidney$spe.num<-info.jump.3RWGD.fish.kidney$internal.events.fish-(info.jump.3RWGD.fish.kidney$dup.num + info.jump.3RWGD.fish.kidney$WGD.num)
info.jump.3RWGD.fish.kidney<-info.jump.3RWGD.fish.kidney[c(11,14,3,12,5,4,13)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.3RWGD.fish.kidney<-
  test.output.3RWGD.kidney.new.fish %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.3RWGD.fish.kidney)<-c("index.old","jump.dup")

test.WGD.3RWGD.fish.kidney<-
  test.output.3RWGD.kidney.new.fish %>%
  filter(node.event=="FishWGD" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.WGD.3RWGD.fish.kidney)<-c("index.old","jump.WGD")

test.spe.3RWGD.fish.kidney<-
  test.output.3RWGD.kidney.new.fish%>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.3RWGD.fish.kidney)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.3RWGD.kidney.fish1<-full_join(test.spe.3RWGD.fish.kidney,test.dup.3RWGD.fish.kidney)
merged.3RWGD.kidney.fish<-full_join(merged.3RWGD.kidney.fish1,test.WGD.3RWGD.fish.kidney)
info.jump.kidney.final.3RWGD.fish<-merge(info.jump.3RWGD.fish.kidney,merged.3RWGD.kidney.fish, by = c("index.old"))
info.jump.kidney.final.3RWGD.fish[is.na(info.jump.kidney.final.3RWGD.fish)] <- 0
rm(merged.3RWGD.kidney.fish1)

## Have to consider trees where atleast one WGD and one speciation exists
info.jump.kidney.final1.3RWGD.fish<-info.jump.kidney.final.3RWGD.fish[which(info.jump.kidney.final.3RWGD.fish$spe.num>0 & info.jump.kidney.final.3RWGD.fish$WGD.num>0),]
rm(info.jump.kidney.final.3RWGD.fish)

## Proportion of jumps
info.jump.kidney.final1.3RWGD.fish$dup.jump.prob<-(info.jump.kidney.final1.3RWGD.fish$jump.dup/(info.jump.kidney.final1.3RWGD.fish$dup.num*2))
info.jump.kidney.final1.3RWGD.fish$spe.jump.prob<-(info.jump.kidney.final1.3RWGD.fish$jump.spe/(info.jump.kidney.final1.3RWGD.fish$spe.num*2))
info.jump.kidney.final1.3RWGD.fish$WGD.jump.prob<-(info.jump.kidney.final1.3RWGD.fish$jump.WGD/(info.jump.kidney.final1.3RWGD.fish$WGD.num*2))

## Since,many WGD trees didnot have a duplication node, it produces NA for dup.jump.prob
## To get rid of this
info.jump.kidney.final2.3RWGD.fish<-info.jump.kidney.final1.3RWGD.fish[complete.cases(info.jump.kidney.final1.3RWGD.fish),] #97 obs

## Stats
median(na.omit(info.jump.kidney.final2.3RWGD.fish$dup.jump.prob)) #0
median(na.omit(info.jump.kidney.final2.3RWGD.fish$spe.jump.prob)) #0.125
median(na.omit(info.jump.kidney.final2.3RWGD.fish$WGD.jump.prob)) #0
mean(na.omit(info.jump.kidney.final2.3RWGD.fish$dup.jump.prob)) #0.2423
mean(na.omit(info.jump.kidney.final2.3RWGD.fish$spe.jump.prob)) #0.2327
mean(na.omit(info.jump.kidney.final2.3RWGD.fish$WGD.jump.prob)) #0.1512

pval.jump.3RWGD.kidney.fish.spe2dup <- paired.wilcox(info.jump.kidney.final2.3RWGD.fish$spe.jump.prob,info.jump.kidney.final2.3RWGD.fish$dup.jump.prob) #p-value =9.04e-1
pval.jump.3RWGD.kidney.fish.spe2WGD <- paired.wilcox(info.jump.kidney.final2.3RWGD.fish$spe.jump.prob,info.jump.kidney.final2.3RWGD.fish$WGD.jump.prob) #p-value = 9.17e-7

## For median data
dtfr.3RWGD.kidney.fish<-NULL
dtfr.3RWGD.kidney.fish<-data.frame(prop.jump=info.jump.kidney.final2.3RWGD.fish$spe.jump.prob,events="speciation")
dtfr.3RWGD.kidney.fish<-rbind(dtfr.3RWGD.kidney.fish,data.frame(prop.jump=info.jump.kidney.final2.3RWGD.fish$dup.jump.prob,events="duplication"))
dtfr.3RWGD.kidney.fish<-rbind(dtfr.3RWGD.kidney.fish,data.frame(prop.jump=info.jump.kidney.final2.3RWGD.fish$WGD.jump.prob,events="FishWGD"))
dtfr.3RWGD.kidney.fish$Exp.abs <-abs(dtfr.3RWGD.kidney.fish$prop.jump)
dtfr.3RWGD.kidney.fish$Event<-factor(dtfr.3RWGD.kidney.fish$events, levels=c("speciation","duplication","FishWGD"))


median.data.jump.3RWGD.kidney.fish<-NULL
median.data.jump.3RWGD.kidney.fish<- dtfr.3RWGD.kidney.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.3RWGD.kidney.fish$pic.round<-paste0("Median = ",round(median.data.jump.3RWGD.kidney.fish$kidney.abs,4), sep="")
median.data.jump.3RWGD.kidney.fish$count<- paste0(median.data.jump.3RWGD.kidney.fish$Event,"\n","(n = ",median.data.jump.3RWGD.kidney.fish$freq,")")

plot3R2<-boxplot.new.3RWGD2(dtfr.3RWGD.kidney.fish,pval.jump.3RWGD.kidney.fish.spe2dup,pval.jump.3RWGD.kidney.fish.spe2WGD,"Exp.3RWGD",median.data.jump.3RWGD.kidney.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("Jump in exp levels for 5 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.05, xmax = 1.7, ymin = 0.67, ymax =0.67, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.4, y = 0.75, label= pval.jump.3RWGD.kidney.fish.spe2dup, fontface = 4)+
  annotate("rect", xmin = 1.05, xmax = 2.8, ymin = 0.88, ymax =0.88, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.7, y = 0.95, label= pval.jump.3RWGD.kidney.fish.spe2WGD, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.3RWGD.kidney.fish$count) ## check

########### Testing OC for all fish specific nodes of 495 WGD trees supporting jump##################
nodes.contrast.3RWGD.kidney.jumps.fish<-nodes.contrast.3RWGD.kidney.jumps.15spe[which(nodes.contrast.3RWGD.kidney.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 495 WGD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.3RWGD.kidney.jumps.fish.removing.jump.node<-nodes.contrast.3RWGD.kidney.jumps.15spe.removing.jump.node[which(nodes.contrast.3RWGD.kidney.jumps.15spe.removing.jump.node$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

########### Testing OC for all fish specific nodes of 323 3RWGD trees that do not support trait jump for teleosts ##################
nodes.contrast.3RWGD.kidney.nojumps.fish<-nodes.contrast.3RWGD.kidney.nojumps.15spe[which(nodes.contrast.3RWGD.kidney.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")

################## Chunk 5: Analysis on contrasts standardized sure SSD trees for traits ##########################

### Step 1a: Considering (sure.SSD.trees) 4139 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for kidney tissue for sure SSD trees 
count<-0
Calibrated.kidney.SSD1<-lapply(sure.SSD.trees, diagnostic.plot.test.exp.all.tissue, "TPM.kidney")
Calibrated.kidney.SSD2<-Calibrated.kidney.SSD1[!is.na(Calibrated.kidney.SSD1)] 
Calibrated.standardized.SSD.kidney<-Calibrated.kidney.SSD2[ ! sapply(Calibrated.kidney.SSD2, is.null) ]## finally we obtained 2701/4139 trees passing the diagnostic tests
rm(Calibrated.kidney.SSD1)
rm(Calibrated.kidney.SSD2)

### Step 1b: Considering (Calibrated.standardized.SSD.kidney) 2701 trees passing the diagnostic test for kidney expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.standardized.SSD.kidney,tree.index.collect))
rm(count)

test.output.SSD.kidney<-NULL
test.output.SSD.kidney<-merge(index.info,test.output.kidney1,by="Gene") ## 11553 observations
index.trees.jump<-unique(test.output.SSD.kidney$tree.num) ##1788 trees
Calibrated.SSD.kidney.jump<-Calibrated.standardized.SSD.kidney[c(index.trees.jump)] ##1788/2701=66.2% Brownian trees supported trait jump model for kidney expressions
Calibrated.SSD.kidney.nojump<-Calibrated.standardized.SSD.kidney[-c(index.trees.jump)] ## 913/2701=33.8%  trees with no support for jump in tau
rm(index.info)
rm(index.trees.jump)

### Step2: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.SSD.kidney))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.SSD.kidney)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 11553 obs
}

test.output.SSD.kidney.processed<-bind_cols(test.output.SSD.kidney,processed.fout) ## 11553 obs for all the 15 species
test.output.SSD.kidney.processed<-test.output.SSD.kidney.processed[(which(!(test.output.SSD.kidney.processed$clade.name %in% c("Clupeocephala","Osteoglossocephalai")))),] ##11036 obs
rm(i)
rm(test.output.SSD.kidney)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

## Step 3: Analysis on all 15 vertebrates clades
test.output.SSD.kidney.processed.mod<-unique(test.output.SSD.kidney.processed[c(6,3,7,8,10,2,1)])  #8242 obs

########### Proportion of jump events in 15 vertebrates species in the 1788 SSD trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.SSD.kidney<-NULL
count<-0
info.jump.SSD.kidney<-bind_rows(lapply(Calibrated.SSD.kidney.jump, tree.data.collection))
info.jump.SSD.kidney$spe.num<-info.jump.SSD.kidney$internal.events-info.jump.SSD.kidney$dup.num
info.jump.SSD.kidney<-info.jump.SSD.kidney[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.SSD.kidney<-
  test.output.SSD.kidney.processed %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.SSD.kidney)<-c("index.old","jump.dup")

test.spe.SSD.kidney<-
  test.output.SSD.kidney.processed%>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.SSD.kidney)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.SSD.kidney<-full_join(test.spe.SSD.kidney,test.dup.SSD.kidney)
info.jump.kidney.final.SSD.all<-merge(info.jump.SSD.kidney,merged.SSD.kidney, by = c("index.old"))
info.jump.kidney.final.SSD.all[is.na(info.jump.kidney.final.SSD.all)] <- 0

## Have to consider trees where atleast one duplication and speciation exists
info.jump.kidney.final1.SSD.all<-info.jump.kidney.final.SSD.all[which(info.jump.kidney.final.SSD.all$dup.num>0 & info.jump.kidney.final.SSD.all$spe.num>0),]
rm(info.jump.kidney.final.SSD.all)

## Proportion of jumps
info.jump.kidney.final1.SSD.all$dup.jump.prob<-(info.jump.kidney.final1.SSD.all$jump.dup/(info.jump.kidney.final1.SSD.all$dup.num*2))
info.jump.kidney.final1.SSD.all$spe.jump.prob<-(info.jump.kidney.final1.SSD.all$jump.spe/(info.jump.kidney.final1.SSD.all$spe.num*2))

## Stats
median(na.omit(info.jump.kidney.final1.SSD.all$dup.jump.prob)) #0.25
median(na.omit(info.jump.kidney.final1.SSD.all$spe.jump.prob)) #0.0769
mean(na.omit(info.jump.kidney.final1.SSD.all$dup.jump.prob)) #0.3442
mean(na.omit(info.jump.kidney.final1.SSD.all$spe.jump.prob)) #0.1445

pval.jump.SSD.kidney <- paired.wilcox(info.jump.kidney.final1.SSD.all$spe.jump.prob,info.jump.kidney.final1.SSD.all$dup.jump.prob) #p-value = 2.2e-16

########### Testing OC for all nodes of 1788 SSD trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.SSD.kidney.jump<-summary.tree(Calibrated.SSD.kidney.jump)
summary.SSD.kidney.jump$Exp.abs <- abs(summary.SSD.kidney.jump$pic.kidney)
summary.SSD.kidney.jumps<-summary.SSD.kidney.jump[,!(names(summary.SSD.kidney.jump) %in% c("pic.kidney"))]
nodes.contrast.SSD.kidney.jumps.15spe <- summary.SSD.kidney.jumps[which(!is.na(summary.SSD.kidney.jumps$Exp.abs)),]
nodes.contrast.SSD.kidney.jumps.15spe$Event <- factor(nodes.contrast.SSD.kidney.jumps.15spe$events, levels=c("speciation", "duplication"))
nodes.contrast.SSD.kidney.jumps.15spe<-nodes.contrast.SSD.kidney.jumps.15spe[(which(!(nodes.contrast.SSD.kidney.jumps.15spe$label %in% c("Clupeocephala","Osteoglossocephalai")))),]
rm(summary.SSD.kidney.jump)
rm(summary.SSD.kidney.jumps)

############ Testing OC using all 15 species for 1788 SSD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.SSD.kidney.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.SSD.kidney.jumps.15spe$index<-nodes.contrast.SSD.kidney.jumps.15spe$index.tree
test.output.SSD.kidney.processed$node<-test.output.SSD.kidney.processed$From
nodes.contrast.SSD.kidney.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.SSD.kidney.jumps.15spe,test.output.SSD.kidney.processed,by=c("index","node"))
nodes.contrast.SSD.kidney.jumps.15spe.removing.jump.node <- nodes.contrast.SSD.kidney.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.SSD.kidney.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.SSD.kidney.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.SSD.kidney.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.SSD.kidney.jumps.15spe.removing.jump.node$Exp.abs<-nodes.contrast.SSD.kidney.jumps.15spe.removing.jump.node$pic
nodes.contrast.SSD.kidney.jumps.15spe.removing.jump.node<-nodes.contrast.SSD.kidney.jumps.15spe.removing.jump.node[(which(!(nodes.contrast.SSD.kidney.jumps.15spe.removing.jump.node$label %in% c("Clupeocephala","Osteoglossocephalai")))),]

########### Testing OC for all nodes of 913 SSD exp trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.SSD.kidney.nojump<-summary.tree(Calibrated.SSD.kidney.nojump)
summary.SSD.kidney.nojump$Exp.abs <- abs(summary.SSD.kidney.nojump$pic)
summary.SSD.kidney.nojumps<-summary.SSD.kidney.nojump[,!(names(summary.SSD.kidney.nojump) %in% c("pic"))]
nodes.contrast.SSD.kidney.nojumps.15spe <- summary.SSD.kidney.nojumps[which(!is.na(summary.SSD.kidney.nojumps$Exp.abs)),]
nodes.contrast.SSD.kidney.nojumps.15spe$Event <- factor(nodes.contrast.SSD.kidney.nojumps.15spe$events, levels=c("speciation", "duplication"))
nodes.contrast.SSD.kidney.nojumps.15spe<-nodes.contrast.SSD.kidney.nojumps.15spe[(which(!(nodes.contrast.SSD.kidney.nojumps.15spe$label %in% c("Clupeocephala","Osteoglossocephalai")))),] 
rm(summary.SSD.kidney.nojump)
rm(summary.SSD.kidney.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.SSD.kidney.new.fish<-test.output.SSD.kidney.processed[which(test.output.SSD.kidney.processed$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##3568 obs

## Step 6: Analysis on fish specific clades
test.output.SSD.kidney.new.fish.mod<-unique(test.output.SSD.kidney.new.fish[c(6,3,7,8,10,2,1)])  #2618 obs


########### Proportion of jump events in teleosts in the 1788 SSD trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.SSD.kidney.trees.fish<-NULL
count<-0
info.jump.SSD.kidney.trees.fish<-bind_rows(lapply(Calibrated.SSD.kidney.jump, tree.data.collection.fish))
info.jump.SSD.kidney.trees.fish$spe.num<-info.jump.SSD.kidney.trees.fish$internal.events.fish-info.jump.SSD.kidney.trees.fish$dup.num
info.jump.SSD.kidney.trees.fish<-info.jump.SSD.kidney.trees.fish[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.SSD.kidney.fish<-
  test.output.SSD.kidney.new.fish  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.SSD.kidney.fish)<-c("index.old","jump.dup")


test.spe.SSD.kidney.fish<-
  test.output.SSD.kidney.new.fish %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.SSD.kidney.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.SSD.kidney.fish<-full_join(test.spe.SSD.kidney.fish,test.dup.SSD.kidney.fish)
info.jump.trees.SSD.kidney.final.fish<-merge(info.jump.SSD.kidney.trees.fish,merged.SSD.kidney.fish, by = c("index.old"))
info.jump.trees.SSD.kidney.final.fish[is.na(info.jump.trees.SSD.kidney.final.fish)] <- 0

## Have to consider trees where atleast one fish specific duplication and speciation exists
info.jump.trees.SSD.kidney.final1.fish<-info.jump.trees.SSD.kidney.final.fish[which(info.jump.trees.SSD.kidney.final.fish$dup.num>0 & info.jump.trees.SSD.kidney.final.fish$spe.num>0),]

## Proportion of jumps
info.jump.trees.SSD.kidney.final1.fish$dup.jump.prob<-(info.jump.trees.SSD.kidney.final1.fish$jump.dup/(info.jump.trees.SSD.kidney.final1.fish$dup.num*2))
info.jump.trees.SSD.kidney.final1.fish$spe.jump.prob<-(info.jump.trees.SSD.kidney.final1.fish$jump.spe/(info.jump.trees.SSD.kidney.final1.fish$spe.num*2))


median(na.omit(info.jump.trees.SSD.kidney.final1.fish$dup.jump.prob)) #0.25
median(na.omit(info.jump.trees.SSD.kidney.final1.fish$spe.jump.prob)) #0.125
mean(na.omit(info.jump.trees.SSD.kidney.final1.fish$dup.jump.prob)) #0.3303
mean(na.omit(info.jump.trees.SSD.kidney.final1.fish$spe.jump.prob)) #0.1968

pval.jump.SSD.kidney.fish <- paired.wilcox(info.jump.trees.SSD.kidney.final1.fish$spe.jump.prob,info.jump.trees.SSD.kidney.final1.fish$dup.jump.prob) #p-value = 2.2e-16


########### Testing OC for all fish specific nodes of 1788 SSD trees supporting jump##################
nodes.contrast.SSD.kidney.jumps.fish<-nodes.contrast.SSD.kidney.jumps.15spe[which(nodes.contrast.SSD.kidney.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] #10684 obs

############ Testing OC for teleosts in 1788 SSD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 

nodes.contrast.SSD.kidney.jumps.fish.removing.jump.node <- NULL
nodes.contrast.SSD.kidney.jumps.fish$index<-nodes.contrast.SSD.kidney.jumps.fish$index.tree
test.output.SSD.kidney.new.fish$node<-test.output.SSD.kidney.new.fish$From
nodes.contrast.SSD.kidney.jumps.fish.removing.jump.node <-anti_join(nodes.contrast.SSD.kidney.jumps.fish,test.output.SSD.kidney.new.fish,by=c("index","node"))
nodes.contrast.SSD.kidney.jumps.fish.removing.jump.node <- nodes.contrast.SSD.kidney.jumps.fish.removing.jump.node[which(!is.na(nodes.contrast.SSD.kidney.jumps.fish.removing.jump.node$pic)),]
nodes.contrast.SSD.kidney.jumps.fish.removing.jump.node$Event <- factor(nodes.contrast.SSD.kidney.jumps.fish.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.SSD.kidney.jumps.fish.removing.jump.node$Exp.abs<-nodes.contrast.SSD.kidney.jumps.fish.removing.jump.node$pic
nodes.contrast.SSD.kidney.jumps.fish.removing.jump.node<-nodes.contrast.SSD.kidney.jumps.fish.removing.jump.node[(which(!(nodes.contrast.SSD.kidney.jumps.fish.removing.jump.node$label %in% c("Clupeocephala","Osteoglossocephalai")))),]

########### Testing OC for all fish specific nodes of 913 SSD exp trees that do not support trait jump for teleosts ##################
nodes.contrast.SSD.kidney.nojumps.fish<-nodes.contrast.SSD.kidney.nojumps.15spe[which(nodes.contrast.SSD.kidney.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] #1886 obs

#save.image("norm_exp_15spe_4tips_new.RData")


