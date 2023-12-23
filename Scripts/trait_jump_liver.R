

### This code is written to work with "trait jump model" for liver

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

test.output.liver<-read.csv(paste0(folder1,"To_pablo/liver/results_liver.csv", sep=""),sep=",",header=T) ##209958 obs

## Only considering trees with posterior probability of jump more than 70%
test.output.liver1<-test.output.liver[which(test.output.liver$JumpProbability>=0.7),] ##25459 obs

########################## Chunk 3: Analysis on all contrasts standardized trees for traits ###############################
### Step 1a: Considering (trees.all.interest) 6923 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for liver tissue for all trees 
count<-0
Calibrated.liver.exp1<-lapply(trees.all.interest, diagnostic.plot.test.exp.all.tissue, "TPM.liver")
Calibrated.liver.exp2<-Calibrated.liver.exp1[!is.na(Calibrated.liver.exp1)] 
Calibrated.all.standardized.liver<-Calibrated.liver.exp2[ ! sapply(Calibrated.liver.exp2, is.null) ]## finally we obtained 5328/6923 tree data passing the diagnostic tests
rm(Calibrated.liver.exp1)
rm(Calibrated.liver.exp2)

### Step 1b: Considering (Calibrated.all.standardized.liver) 5328 tree data passing the diagnostic test for liver expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.all.standardized.liver,tree.index.collect))
rm(count)

test.output.liver.all<-NULL
test.output.liver.all<-merge(index.info,test.output.liver1,by="Gene") ## 17514 observations
index.trees.jump<-unique(test.output.liver.all$tree.num) ##4480 trees
Calibrated.standardized.liver.jump<-Calibrated.all.standardized.liver[c(index.trees.jump)] ##4480/5328=84.08% Brownian trees supported trait jump model for liver expressions 
Calibrated.standardized.liver.nojump<-Calibrated.all.standardized.liver[-c(index.trees.jump)] ## 848/5328=15.92% trees with no support for jump in liver expression
rm(index.info)
rm(index.trees.jump)

### Step2: Calculating PICs of required trait for trees to match the index number of the trees on which levolution was run
Count<-0
trees.all.interest<-lapply(trees.all.interest,contrast.calc,trait="TPM.liver")
# Count<-0
# Calibrated.standardized.liver.jump<-lapply(Calibrated.standardized.liver.jump,contrast.calc,trait="TPM.liver") 
# rm(Count)

### Step3: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.liver.all))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.liver.all)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 17514 obs
}
test.output.liver.processed<-bind_cols(test.output.liver.all,processed.fout) ## 17514 obs for all the 15 species
rm(i)
rm(test.output.liver.all)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

## Step 4: Analysis on all 15 vertebrates clades
test.output.liver.new.mod<-unique(test.output.liver.processed[c(6,3,7,8,10,2,1)])  #14716 obs

########### Proportion of jump events in 15 vertebrates species in the 4480 trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.trees.liver<-NULL
count<-0
info.jump.trees.liver<-bind_rows(lapply(Calibrated.standardized.liver.jump, tree.data.collection))
info.jump.trees.liver$spe.num<-info.jump.trees.liver$internal.events-info.jump.trees.liver$dup.num
info.jump.trees.liver<-info.jump.trees.liver[c(11,12,3:5)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.liver<-
  test.output.liver.processed  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.liver)<-c("index.old","jump.dup")

test.spe.liver<-
  test.output.liver.processed %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.liver)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.liver<-full_join(test.spe.liver,test.dup.liver)
info.jump.trees.liver.final<-merge(info.jump.trees.liver,merged.liver, by = c("index.old"))
info.jump.trees.liver.final[is.na(info.jump.trees.liver.final)] <- 0

## Have to consider trees where atleast one duplication and speciation exists
info.jump.trees.liver.final1<-info.jump.trees.liver.final[which(info.jump.trees.liver.final$dup.num>0 & info.jump.trees.liver.final$spe.num>0),]
rm(info.jump.trees.liver.final)

## Proportion of jumps
info.jump.trees.liver.final1$dup.jump.prob<-(info.jump.trees.liver.final1$jump.dup/(info.jump.trees.liver.final1$dup.num*2))
info.jump.trees.liver.final1$spe.jump.prob<-(info.jump.trees.liver.final1$jump.spe/(info.jump.trees.liver.final1$spe.num*2))

## Stats
median(na.omit(info.jump.trees.liver.final1$dup.jump.prob)) #0
median(na.omit(info.jump.trees.liver.final1$spe.jump.prob)) #0.059
mean(na.omit(info.jump.trees.liver.final1$dup.jump.prob)) #0.1628
mean(na.omit(info.jump.trees.liver.final1$spe.jump.prob)) #0.0978

## surprisingly jump in dup > jump in spe based on the values
pval.jump.liver <- paired.wilcox(info.jump.trees.liver.final1$spe.jump.prob,info.jump.trees.liver.final1$dup.jump.prob) #p-value = 4.47e-7

## For median data
dtfr.liver<-NULL
dtfr.liver<-data.frame(prop.jump=info.jump.trees.liver.final1$spe.jump.prob,events="speciation")
dtfr.liver<-rbind(dtfr.liver,data.frame(prop.jump=info.jump.trees.liver.final1$dup.jump.prob,events="duplication"))
dtfr.liver$Exp.abs <-abs(dtfr.liver$prop.jump)
dtfr.liver$Event<-factor(dtfr.liver$events, levels=c("speciation","duplication"))

library(tidyverse)
median.data.jump.liver<-NULL
median.data.jump.liver<- dtfr.liver%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.liver$pic.round<-paste0("Median = ",round(median.data.jump.liver$Exp.abs,4), sep="")
median.data.jump.liver$count<- paste0(median.data.jump.liver$Event,"\n\n","(n = ",median.data.jump.liver$freq,")")

plot4A<-boxplot.new2(dtfr.liver,pval.jump.liver,"Exp",median.data.jump.liver) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.liver, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.liver$count) ## check

# Create Data
data.jump.liver <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(4480,848)
)

# Compute the position of labels
data.jump.liver <- data.jump.liver %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.liver$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.liver$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.liver<-ggplot(data.jump.liver, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 4480 trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.liver.jump<-summary.tree(Calibrated.standardized.liver.jump)
summary.liver.jump$Exp.abs <- abs(summary.liver.jump$pic.liver)
summary.liver.jumps<-summary.liver.jump[,!(names(summary.liver.jump) %in% c("pic.liver"))]
nodes.contrast.liver.jumps.15spe <- summary.liver.jumps[which(!is.na(summary.liver.jumps$Exp.abs)),]
nodes.contrast.liver.jumps.15spe$Event <- factor(nodes.contrast.liver.jumps.15spe$events, levels=c("speciation", "duplication"))
rm(summary.liver.jump)
rm(summary.liver.jumps)

############ Testing OC using all 15 species for 4480 trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.liver.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.liver.jumps.15spe$index<-nodes.contrast.liver.jumps.15spe$index.tree
test.output.liver.processed$node<-test.output.liver.processed$From
nodes.contrast.liver.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.liver.jumps.15spe,test.output.liver.processed,by=c("index","node"))
nodes.contrast.liver.jumps.15spe.removing.jump.node <- nodes.contrast.liver.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.liver.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.liver.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.liver.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.liver.jumps.15spe.removing.jump.node$Exp.abs<-nodes.contrast.liver.jumps.15spe.removing.jump.node$pic


########### Testing OC for all nodes of 848 trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.liver.nojump<-summary.tree(Calibrated.standardized.liver.nojump)
summary.liver.nojump$Exp.abs <- abs(summary.liver.nojump$pic.liver)
summary.liver.nojumps<-summary.liver.nojump[,!(names(summary.liver.nojump) %in% c("pic.liver"))]
nodes.contrast.liver.nojumps.15spe <- summary.liver.nojump[which(!is.na(summary.liver.nojump$Exp.abs)),]
nodes.contrast.liver.nojumps.15spe$Event <- factor(nodes.contrast.liver.nojumps.15spe$events, levels=c("speciation", "duplication"))
rm(summary.liver.nojump)
rm(summary.liver.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.liver.new.fish<-test.output.liver.processed[which(test.output.liver.processed$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##8677 obs

## Step 6: Analysis on fish specific clades
test.output.liver.new.fish.mod<-unique(test.output.liver.new.fish[c(6,3,7,8,10,2,1)])  #7581 obs

########### Proportion of jump events in teleosts in the 4480 trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.trees.liver.fish<-NULL
count<-0
info.jump.trees.liver.fish<-bind_rows(lapply(Calibrated.standardized.liver.jump, tree.data.collection.fish))
info.jump.trees.liver.fish$spe.num<-info.jump.trees.liver.fish$internal.events.fish-info.jump.trees.liver.fish$dup.num
info.jump.trees.liver.fish<-info.jump.trees.liver.fish[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.liver.fish<-
  test.output.liver.new.fish  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.liver.fish)<-c("index.old","jump.dup")

test.spe.liver.fish<-
  test.output.liver.new.fish%>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.liver.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.liver.fish<-full_join(test.spe.liver.fish,test.dup.liver.fish)
info.jump.trees.liver.final.fish<-merge(info.jump.trees.liver.fish,merged.liver.fish, by = c("index.old"))
info.jump.trees.liver.final.fish[is.na(info.jump.trees.liver.final.fish)] <- 0

## Have to consider trees where atleast one fish specific duplication and speciation exists
info.jump.trees.liver.final1.fish<-info.jump.trees.liver.final.fish[which(info.jump.trees.liver.final.fish$dup.num>0 & info.jump.trees.liver.final.fish$spe.num>0),]

## Proportion of jumps
info.jump.trees.liver.final1.fish$dup.jump.prob<-(info.jump.trees.liver.final1.fish$jump.dup/(info.jump.trees.liver.final1.fish$dup.num*2))
info.jump.trees.liver.final1.fish$spe.jump.prob<-(info.jump.trees.liver.final1.fish$jump.spe/(info.jump.trees.liver.final1.fish$spe.num*2))


median(na.omit(info.jump.trees.liver.final1.fish$dup.jump.prob)) #0
median(na.omit(info.jump.trees.liver.final1.fish$spe.jump.prob)) #0.125
mean(na.omit(info.jump.trees.liver.final1.fish$dup.jump.prob)) #0.1232
mean(na.omit(info.jump.trees.liver.final1.fish$spe.jump.prob)) #0.1631

pval.jump.liver.fish <- paired.wilcox(info.jump.trees.liver.final1.fish$spe.jump.prob,info.jump.trees.liver.final1.fish$dup.jump.prob) #p-value = 2.2e-16


## For median data
dtfr.liver.fish<-NULL
dtfr.liver.fish<-data.frame(prop.jump=info.jump.trees.liver.final1.fish$spe.jump.prob,events="speciation")
dtfr.liver.fish<-rbind(dtfr.liver.fish,data.frame(prop.jump=info.jump.trees.liver.final1.fish$dup.jump.prob,events="duplication"))
dtfr.liver.fish$Exp.abs <-abs(dtfr.liver.fish$prop.jump)
dtfr.liver.fish$Event<-factor(dtfr.liver.fish$events, levels=c("speciation","duplication"))


library(tidyverse)
median.data.jump.liver.liver.fish<-NULL
median.data.jump.liver.liver.fish<- dtfr.liver.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.liver.liver.fish$pic.round<-paste0("Median = ",round(median.data.jump.liver.liver.fish$Exp.abs,4), sep="")
median.data.jump.liver.liver.fish$count<- paste0(median.data.jump.liver.liver.fish$Event,"\n\n","(n = ",median.data.jump.liver.liver.fish$freq,")")

plot4B<-boxplot.new2(dtfr.liver.fish,pval.jump.liver.fish,"Exp",median.data.jump.liver.liver.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 5 teleosts"))))+
  #ylab(expression(bold(paste("Proportions of jump events for teleosts ")))) +
  ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.liver.fish, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.liver.liver.fish$count) ## check

########### Testing OC for all fish specific nodes of 4480 trees supporting jump##################
nodes.contrast.liver.jumps.fish<-nodes.contrast.liver.jumps.15spe[which(nodes.contrast.liver.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 4480 trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.liver.jumps.fish.removing.jump.node <- NULL
nodes.contrast.liver.jumps.fish$index<-nodes.contrast.liver.jumps.fish$index.tree
test.output.liver.new.fish$node<-test.output.liver.new.fish$From
nodes.contrast.liver.jumps.fish.removing.jump.node <-anti_join(nodes.contrast.liver.jumps.fish,test.output.liver.new.fish,by=c("index","node"))
nodes.contrast.liver.jumps.fish.removing.jump.node <- nodes.contrast.liver.jumps.fish.removing.jump.node[which(!is.na(nodes.contrast.liver.jumps.fish.removing.jump.node$pic)),]
nodes.contrast.liver.jumps.fish.removing.jump.node$Event <- factor(nodes.contrast.liver.jumps.fish.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.liver.jumps.fish.removing.jump.node$Exp.abs<-nodes.contrast.liver.jumps.fish.removing.jump.node$pic

rm(test.output.liver.new.fish)

########### Testing OC for all fish specific nodes of 848 trees that do not support trait jump for teleosts ##################
nodes.contrast.liver.nojumps.fish<-nodes.contrast.liver.nojumps.15spe[which(nodes.contrast.liver.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")

################## Chunk 4: Analysis on contrasts standardized sure 3RWGD trees for traits ##########################

### Step 1a: Considering (Sure.WGD.trees) 1159 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for liver tissue for sure WGD trees 
count<-0
Calibrated.liver.3Rexp1<-lapply(Sure.WGD.trees, diagnostic.plot.test.exp.all.tissue, "TPM.liver")
Calibrated.liver.3Rexp2<-Calibrated.liver.3Rexp1[!is.na(Calibrated.liver.3Rexp1)] 
Calibrated.3R.standardized.liver<-Calibrated.liver.3Rexp2[ ! sapply(Calibrated.liver.3Rexp2, is.null) ]## finally we obtained 892/1159 trees passing the diagnostic tests
rm(Calibrated.liver.3Rexp1)
rm(Calibrated.liver.3Rexp2)

### Step 1b: Considering (Calibrated.3R.standardized.liver) 892 trees passing the diagnostic test for liver expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.3R.standardized.liver,tree.index.collect))
rm(count)

test.output.3Rexp.liver<-NULL
test.output.3Rexp.liver<-merge(index.info,test.output.liver1,by="Gene") ## 2695 observations
index.trees.jump<-unique(test.output.3Rexp.liver$tree.num) ##769 trees
Calibrated.3RWGD.liver.jump<-Calibrated.3R.standardized.liver[c(index.trees.jump)] ##769/892=86.21% Brownian trees supported trait jump model for liver expressions
Calibrated.3RWGD.liver.nojump<-Calibrated.3R.standardized.liver[-c(index.trees.jump)] ## 123/892=13.79% trees with no support for jump in trait
rm(index.info)
rm(index.trees.jump)

### Step2: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.3Rexp.liver))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.3Rexp.liver)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 2695 obs
}

test.output.3RWGD.liver.processed<-bind_cols(test.output.3Rexp.liver,processed.fout) ## 2695 obs for all the 15 species
rm(i)
rm(test.output.3Rexp.liver)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

##Since our calibrated tree shows both the FishWGD and duplication nodes, 
##we want to modify the node events of the jump output to separate out the effect of FishWGD from the SSDs in those trees.

jump.3RWGD.liver.trees.summary<-summary.tree(Calibrated.3RWGD.liver.jump)
jump.3RWGD.liver.trees.summary$Gene<-jump.3RWGD.liver.trees.summary$index.tree
jump.3RWGD.liver.trees.summary$From<-jump.3RWGD.liver.trees.summary$node
test.output.3RWGD.liver.processed.mod<-merge(test.output.3RWGD.liver.processed,jump.3RWGD.liver.trees.summary,by=c("Gene","From"))
test.output.3RWGD.liver.processed.mod<-test.output.3RWGD.liver.processed.mod[c(1,3,2,4:11,16)]
test.output.3RWGD.liver.processed.mod$node.event<-test.output.3RWGD.liver.processed.mod$events
colnames(test.output.3RWGD.liver.processed.mod)[10]<-"pic"
test.output.3RWGD.liver.processed.mod<-test.output.3RWGD.liver.processed.mod[c(1:11)] ## 2695 obs

## Step 3: Analysis on all 15 vertebrates clades
test.output.3RWGD.liver.processed.new<-unique(test.output.3RWGD.liver.processed.mod[c(6,3,7,8,10,2,1)])  #2362 obs

########### Proportion of jump events in 15 vertebrates species in the 769 3RWGD trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.3RWGD.liver<-NULL
info.jump.3RWGD.liver<-bind_rows(lapply(Calibrated.3RWGD.liver.jump, tree.data.collection.3RWGD))
info.jump.3RWGD.liver$spe.num<-info.jump.3RWGD.liver$internal.events-(info.jump.3RWGD.liver$dup.num + info.jump.3RWGD.liver$WGD.num)
info.jump.3RWGD.liver<-info.jump.3RWGD.liver[c(11,14,3,12,5,4,13)]


## Counting frquency of jump for "speciation" and "duplication" events
test.dup.3RWGD.liver<-
  test.output.3RWGD.liver.processed.mod %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.3RWGD.liver)<-c("index.old","jump.dup")

test.WGD.3RWGD.liver<-
  test.output.3RWGD.liver.processed.mod %>%
  filter(node.event=="FishWGD" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.WGD.3RWGD.liver)<-c("index.old","jump.WGD")

test.spe.3RWGD.liver<-
  test.output.3RWGD.liver.processed.mod %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.3RWGD.liver)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.3RWGD.liver1<-full_join(test.spe.3RWGD.liver,test.dup.3RWGD.liver)
merged.3RWGD.liver<-full_join(merged.3RWGD.liver1,test.WGD.3RWGD.liver)
info.jump.3RWGD.liver.final1<-merge(info.jump.3RWGD.liver,merged.3RWGD.liver, by = c("index.old"))
info.jump.3RWGD.liver.final1[is.na(info.jump.3RWGD.liver.final1)] <- 0
rm(merged.3RWGD.liver1)

## Have to consider trees where atleast one WGD and one speciation exists
info.jump.3RWGD.liver.final<-info.jump.3RWGD.liver.final1[which(info.jump.3RWGD.liver.final1$spe.num>0 & info.jump.3RWGD.liver.final1$WGD.num>0),]
rm(info.jump.3RWGD.liver.final1)

## Proportion of jumps
info.jump.3RWGD.liver.final$dup.jump.prob<-(info.jump.3RWGD.liver.final$jump.dup/(info.jump.3RWGD.liver.final$dup.num*2))
info.jump.3RWGD.liver.final$spe.jump.prob<-(info.jump.3RWGD.liver.final$jump.spe/(info.jump.3RWGD.liver.final$spe.num*2))
info.jump.3RWGD.liver.final$WGD.jump.prob<-(info.jump.3RWGD.liver.final$jump.WGD/(info.jump.3RWGD.liver.final$WGD.num*2))

## Since,many WGD trees didnot have a duplication node, it produces NA for dup.jump.prob
## To get rid of this
info.jump.3RWGD.liver.final2<-info.jump.3RWGD.liver.final[complete.cases(info.jump.3RWGD.liver.final),] ##227 obs

## Stats
median(info.jump.3RWGD.liver.final2$dup.jump.prob) #0
median(info.jump.3RWGD.liver.final2$spe.jump.prob) #0.055
median(info.jump.3RWGD.liver.final2$WGD.jump.prob) #0
mean(info.jump.3RWGD.liver.final2$dup.jump.prob) #0154
mean(info.jump.3RWGD.liver.final2$spe.jump.prob) #0.09
mean(info.jump.3RWGD.liver.final2$WGD.jump.prob) #0.046

pval.jump.3RWGD.liver.spe2dup <- paired.wilcox(info.jump.3RWGD.liver.final2$spe.jump.prob,info.jump.3RWGD.liver.final2$dup.jump.prob) #p-value = 6.01e-1
pval.jump.3RWGD.liver.spe2WGD <- paired.wilcox(info.jump.3RWGD.liver.final2$spe.jump.prob,info.jump.3RWGD.liver.final2$WGD.jump.prob) #p-value < 2.2e-16

## For median data
dtfr.3RWGD.liver.15spe<-NULL
dtfr.3RWGD.liver.15spe<-data.frame(prop.jump=info.jump.3RWGD.liver.final2$spe.jump.prob,events="speciation")
dtfr.3RWGD.liver.15spe<-rbind(dtfr.3RWGD.liver.15spe,data.frame(prop.jump=info.jump.3RWGD.liver.final2$dup.jump.prob,events="duplication"))
dtfr.3RWGD.liver.15spe<-rbind(dtfr.3RWGD.liver.15spe,data.frame(prop.jump=info.jump.3RWGD.liver.final2$WGD.jump.prob,events="FishWGD"))
dtfr.3RWGD.liver.15spe$liver.abs <-abs(dtfr.3RWGD.liver.15spe$prop.jump)
dtfr.3RWGD.liver.15spe$Event<-factor(dtfr.3RWGD.liver.15spe$events, levels=c("speciation","duplication","FishWGD"))


median.data.jump.3RWGD.liver.15spe<-NULL
median.data.jump.3RWGD.liver.15spe<- dtfr.3RWGD.liver.15spe%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.3RWGD.liver.15spe$pic.round<-paste0("Median = ",round(median.data.jump.3RWGD.liver.15spe$Exp.abs,4), sep="")
median.data.jump.3RWGD.liver.15spe$count<- paste0(median.data.jump.3RWGD.liver.15spe$Event,"\n\n","(n = ",median.data.jump.3RWGD.liver.15spe$freq,")")

plot3R1<-boxplot.new.3RWGD2(dtfr.3RWGD.liver.15spe,pval.jump.3RWGD.liver.spe2dup,pval.jump.3RWGD.liver.spe2WGD,"Exp.3RWGD",median.data.jump.3RWGD.liver.15spe) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.05, xmax = 1.7, ymin = 0.67, ymax =0.67, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.4, y = 0.75, label= pval.jump.3RWGD.liver.spe2dup, fontface = 4)+
  annotate("rect", xmin = 1.05, xmax = 2.8, ymin = 0.88, ymax =0.88, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.7, y = 0.95, label= pval.jump.3RWGD.liver.spe2WGD, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.3RWGD.liver.15spe$count) ## check

# Create Data
data.jump.3RWGD.exp <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(769,123)
)

# Compute the position of labels
data.jump.3RWGD.liver <- data.jump.3RWGD.liver %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.3RWGD.liver$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.3RWGD.liver$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.3RWGD.liver<-ggplot(data.jump.3RWGD.liver, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 769 WGD trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.3RWGD.liver.jump<-summary.tree(Calibrated.3RWGD.liver.jump)
summary.3RWGD.liver.jump$Exp.abs <- abs(summary.3RWGD.liver.jump$pic)
summary.3RWGD.liver.jumps<-summary.3RWGD.liver.jump[,!(names(summary.3RWGD.liver.jump) %in% c("pic"))]
nodes.contrast.3RWGD.liver.jumps.15spe <- summary.3RWGD.liver.jumps[which(!is.na(summary.3RWGD.liver.jumps$Exp.abs)),]
nodes.contrast.3RWGD.liver.jumps.15spe$Event <- factor(nodes.contrast.3RWGD.liver.jumps.15spe$events, levels=c("speciation", "duplication", "FishWGD"))
rm(summary.3RWGD.liver.jump)
rm(summary.3RWGD.liver.jumps)

############ Testing OC using all 15 species for 769 WGD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.3RWGD.liver.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.3RWGD.liver.jumps.15spe$index<-nodes.contrast.3RWGD.liver.jumps.15spe$index.tree
test.output.3RWGD.liver.processed.mod$node<-test.output.3RWGD.liver.processed.mod$From
nodes.contrast.3RWGD.liver.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.3RWGD.liver.jumps.15spe,test.output.3RWGD.liver.processed.mod,by=c("index","node"))
nodes.contrast.3RWGD.liver.jumps.15spe.removing.jump.node <- nodes.contrast.3RWGD.liver.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.3RWGD.liver.jumps.15spe.removing.jump.node$Exp.abs)),]
nodes.contrast.3RWGD.liver.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.3RWGD.liver.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication","FishWGD"))
nodes.contrast.3RWGD.liver.jumps.15spe.removing.jump.node$pic<-abs(nodes.contrast.3RWGD.liver.jumps.15spe.removing.jump.node$Exp.abs)

########### Testing OC for all nodes of 123 3RWGD exp trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.3RWGD.liver.nojump<-summary.tree(Calibrated.3RWGD.liver.nojump)
summary.3RWGD.liver.nojump$Exp.abs <- abs(summary.3RWGD.liver.nojump$pic.liver)
summary.3RWGD.liver.nojumps<-summary.3RWGD.liver.nojump[,!(names(summary.3RWGD.liver.nojump) %in% c("pic.liver"))]
nodes.contrast.3RWGD.liver.nojumps.15spe <- summary.3RWGD.liver.nojumps[which(!is.na(summary.3RWGD.liver.nojumps$Exp.abs)),]
nodes.contrast.3RWGD.liver.nojumps.15spe$Event <- factor(nodes.contrast.3RWGD.liver.nojumps.15spe$events, levels=c("speciation", "duplication","FishWGD"))
rm(summary.3RWGD.liver.nojump)
rm(summary.3RWGD.liver.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.3RWGD.liver.new.fish<-test.output.3RWGD.liver.processed.mod[which(test.output.3RWGD.liver.processed.mod$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##1384 obs

## Step 6: Analysis on fish specific clades
test.output.3RWGD.liver.processed.new.fish<-unique(test.output.3RWGD.liver.new.fish[c(6,3,7,8,10,2,1)])  #1022obs

########### Proportion of jump events in teleosts in the 769 WGD trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.3RWGD.fish.liver<-NULL
info.jump.3RWGD.fish.liver<-bind_rows(lapply(Calibrated.3RWGD.liver.jump, tree.data.collection.3RWGD.fish))
info.jump.3RWGD.fish.liver$spe.num<-info.jump.3RWGD.fish.liver$internal.events.fish-(info.jump.3RWGD.fish.liver$dup.num + info.jump.3RWGD.fish.liver$WGD.num)
info.jump.3RWGD.fish.liver<-info.jump.3RWGD.fish.liver[c(11,14,3,12,5,4,13)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.3RWGD.fish.liver<-
  test.output.3RWGD.liver.new.fish %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.3RWGD.fish.liver)<-c("index.old","jump.dup")

test.WGD.3RWGD.fish.liver<-
  test.output.3RWGD.liver.new.fish %>%
  filter(node.event=="FishWGD" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.WGD.3RWGD.fish.liver)<-c("index.old","jump.WGD")

test.spe.3RWGD.fish.liver<-
  test.output.3RWGD.liver.new.fish %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.3RWGD.fish.liver)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.3RWGD.liver.fish1<-full_join(test.spe.3RWGD.fish.liver,test.dup.3RWGD.fish.liver)
merged.3RWGD.liver.fish<-full_join(merged.3RWGD.liver.fish1,test.WGD.3RWGD.fish.liver)
info.jump.liver.final.3RWGD.fish<-merge(info.jump.3RWGD.fish.liver,merged.3RWGD.liver.fish, by = c("index.old"))
info.jump.liver.final.3RWGD.fish[is.na(info.jump.liver.final.3RWGD.fish)] <- 0
rm(merged.3RWGD.liver.fish1)

## Have to consider trees where atleast one WGD and one speciation exists
info.jump.liver.final1.3RWGD.fish<-info.jump.liver.final.3RWGD.fish[which(info.jump.liver.final.3RWGD.fish$spe.num>0 & info.jump.liver.final.3RWGD.fish$WGD.num>0),]
rm(info.jump.liver.final.3RWGD.fish)

## Proportion of jumps
info.jump.liver.final1.3RWGD.fish$dup.jump.prob<-(info.jump.liver.final1.3RWGD.fish$jump.dup/(info.jump.liver.final1.3RWGD.fish$dup.num*2))
info.jump.liver.final1.3RWGD.fish$spe.jump.prob<-(info.jump.liver.final1.3RWGD.fish$jump.spe/(info.jump.liver.final1.3RWGD.fish$spe.num*2))
info.jump.liver.final1.3RWGD.fish$WGD.jump.prob<-(info.jump.liver.final1.3RWGD.fish$jump.WGD/(info.jump.liver.final1.3RWGD.fish$WGD.num*2))

## Since,many WGD trees didnot have a duplication node, it produces NA for dup.jump.prob
## To get rid of this
info.jump.liver.final2.3RWGD.fish<-info.jump.liver.final1.3RWGD.fish[complete.cases(info.jump.liver.final1.3RWGD.fish),] #206 obs

## Stats
median(na.omit(info.jump.liver.final2.3RWGD.fish$dup.jump.prob)) #0
median(na.omit(info.jump.liver.final2.3RWGD.fish$spe.jump.prob)) #0.1
median(na.omit(info.jump.liver.final2.3RWGD.fish$WGD.jump.prob)) #0
mean(na.omit(info.jump.liver.final2.3RWGD.fish$dup.jump.prob)) #0.1212
mean(na.omit(info.jump.liver.final2.3RWGD.fish$spe.jump.prob)) #0.1473
mean(na.omit(info.jump.liver.final2.3RWGD.fish$WGD.jump.prob)) #0.0461

pval.jump.3RWGD.liver.fish.spe2dup <- paired.wilcox(info.jump.liver.final2.3RWGD.fish$spe.jump.prob,info.jump.liver.final2.3RWGD.fish$dup.jump.prob) #p-value = 1.08e-4
pval.jump.3RWGD.liver.fish.spe2WGD <- paired.wilcox(info.jump.liver.final2.3RWGD.fish$spe.jump.prob,info.jump.liver.final2.3RWGD.fish$WGD.jump.prob) #p-value = 2.2e-16

## For median data
dtfr.3RWGD.liver.fish<-NULL
dtfr.3RWGD.liver.fish<-data.frame(prop.jump=info.jump.liver.final2.3RWGD.fish$spe.jump.prob,events="speciation")
dtfr.3RWGD.liver.fish<-rbind(dtfr.3RWGD.liver.fish,data.frame(prop.jump=info.jump.liver.final2.3RWGD.fish$dup.jump.prob,events="duplication"))
dtfr.3RWGD.liver.fish<-rbind(dtfr.3RWGD.liver.fish,data.frame(prop.jump=info.jump.liver.final2.3RWGD.fish$WGD.jump.prob,events="FishWGD"))
dtfr.3RWGD.liver.fish$Exp.abs <-abs(dtfr.3RWGD.liver.fish$prop.jump)
dtfr.3RWGD.liver.fish$Event<-factor(dtfr.3RWGD.liver.fish$events, levels=c("speciation","duplication","FishWGD"))


median.data.jump.3RWGD.liver.fish<-NULL
median.data.jump.3RWGD.liver.fish<- dtfr.3RWGD.liver.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.3RWGD.liver.fish$pic.round<-paste0("Median = ",round(median.data.jump.3RWGD.liver.fish$liver.abs,4), sep="")
median.data.jump.3RWGD.liver.fish$count<- paste0(median.data.jump.3RWGD.liver.fish$Event,"\n","(n = ",median.data.jump.3RWGD.liver.fish$freq,")")

plot3R2<-boxplot.new.3RWGD2(dtfr.3RWGD.liver.fish,pval.jump.3RWGD.liver.fish.spe2dup,pval.jump.3RWGD.liver.fish.spe2WGD,"Exp.3RWGD",median.data.jump.3RWGD.liver.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("Jump in exp levels for 5 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.05, xmax = 1.7, ymin = 0.67, ymax =0.67, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.4, y = 0.75, label= pval.jump.3RWGD.liver.fish.spe2dup, fontface = 4)+
  annotate("rect", xmin = 1.05, xmax = 2.8, ymin = 0.88, ymax =0.88, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.7, y = 0.95, label= pval.jump.3RWGD.liver.fish.spe2WGD, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.3RWGD.liver.fish$count) ## check

########### Testing OC for all fish specific nodes of 769 WGD trees supporting jump##################
nodes.contrast.3RWGD.liver.jumps.fish<-nodes.contrast.3RWGD.liver.jumps.15spe[which(nodes.contrast.3RWGD.liver.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 769 WGD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.3RWGD.liver.jumps.fish.removing.jump.node<-nodes.contrast.3RWGD.liver.jumps.15spe.removing.jump.node[which(nodes.contrast.3RWGD.liver.jumps.15spe.removing.jump.node$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

########### Testing OC for all fish specific nodes of 123 3RWGD trees that do not support trait jump for teleosts ##################
nodes.contrast.3RWGD.liver.nojumps.fish<-nodes.contrast.3RWGD.liver.nojumps.15spe[which(nodes.contrast.3RWGD.liver.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")

################## Chunk 5: Analysis on contrasts standardized sure SSD trees for traits ##########################

### Step 1a: Considering (sure.SSD.trees) 4139 trees to passing the diagnostic test for specific trait 

## Diagnostic tests for liver tissue for sure SSD trees 
count<-0
Calibrated.liver.SSD1<-lapply(sure.SSD.trees, diagnostic.plot.test.exp.all.tissue, "TPM.liver")
Calibrated.liver.SSD2<-Calibrated.liver.SSD1[!is.na(Calibrated.liver.SSD1)] 
Calibrated.standardized.SSD.liver<-Calibrated.liver.SSD2[ ! sapply(Calibrated.liver.SSD2, is.null) ]## finally we obtained 3182/4139 trees passing the diagnostic tests
rm(Calibrated.liver.SSD1)
rm(Calibrated.liver.SSD2)

### Step 1b: Considering (Calibrated.standardized.SSD.liver) 3182 trees passing the diagnostic test for liver expressions 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.standardized.SSD.liver,tree.index.collect))
rm(count)

test.output.SSD.liver<-NULL
test.output.SSD.liver<-merge(index.info,test.output.liver1,by="Gene") ## 11048 observations
index.trees.jump<-unique(test.output.SSD.liver$tree.num) ##2641 trees
Calibrated.SSD.liver.jump<-Calibrated.standardized.SSD.liver[c(index.trees.jump)] ##2641/3182=83% Brownian trees supported trait jump model for liver expressions
Calibrated.SSD.liver.nojump<-Calibrated.standardized.SSD.liver[-c(index.trees.jump)] ## 541/3182=17%  trees with no support for jump in tau
rm(index.info)
rm(index.trees.jump)

### Step2: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.SSD.liver))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.SSD.liver)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 11048 obs
}

test.output.SSD.liver.processed<-bind_cols(test.output.SSD.liver,processed.fout) ## 11048 obs for all the 15 species
test.output.SSD.liver.processed<-test.output.SSD.liver.processed[(which(!(test.output.SSD.liver.processed$clade.name %in% c("Clupeocephala","Osteoglossocephalai")))),] ##10275 obs
rm(i)
rm(test.output.SSD.liver)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

## Step 3: Analysis on all 15 vertebrates clades
test.output.SSD.liver.processed.mod<-unique(test.output.SSD.liver.processed[c(6,3,7,8,10,2,1)])  #8358 obs

########### Proportion of jump events in 15 vertebrates species in the 2641 SSD trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.SSD.liver<-NULL
count<-0
info.jump.SSD.liver<-bind_rows(lapply(Calibrated.SSD.liver.jump, tree.data.collection))
info.jump.SSD.liver$spe.num<-info.jump.SSD.liver$internal.events-info.jump.SSD.liver$dup.num
info.jump.SSD.liver<-info.jump.SSD.liver[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.SSD.liver<-
  test.output.SSD.liver.processed %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.SSD.liver)<-c("index.old","jump.dup")

test.spe.SSD.liver<-
  test.output.SSD.liver.processed %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.SSD.liver)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.SSD.liver<-full_join(test.spe.SSD.liver,test.dup.SSD.liver)
info.jump.liver.final.SSD.all<-merge(info.jump.SSD.liver,merged.SSD.liver, by = c("index.old"))
info.jump.liver.final.SSD.all[is.na(info.jump.liver.final.SSD.all)] <- 0

## Have to consider trees where atleast one duplication and speciation exists
info.jump.liver.final1.SSD.all<-info.jump.liver.final.SSD.all[which(info.jump.liver.final.SSD.all$dup.num>0 & info.jump.liver.final.SSD.all$spe.num>0),]
rm(info.jump.liver.final.SSD.all)

## Proportion of jumps
info.jump.liver.final1.SSD.all$dup.jump.prob<-(info.jump.liver.final1.SSD.all$jump.dup/(info.jump.liver.final1.SSD.all$dup.num*2))
info.jump.liver.final1.SSD.all$spe.jump.prob<-(info.jump.liver.final1.SSD.all$jump.spe/(info.jump.liver.final1.SSD.all$spe.num*2))

## Stats
median(na.omit(info.jump.liver.final1.SSD.all$dup.jump.prob)) #0
median(na.omit(info.jump.liver.final1.SSD.all$spe.jump.prob)) #0.05
mean(na.omit(info.jump.liver.final1.SSD.all$dup.jump.prob)) #0.2173
mean(na.omit(info.jump.liver.final1.SSD.all$spe.jump.prob)) #0.0996

##Surprisingly jump in dup > jump in spe (according to the median value it should be just opposite)
pval.jump.SSD.liver <- paired.wilcox.one.sided(info.jump.liver.final1.SSD.all$spe.jump.prob,info.jump.liver.final1.SSD.all$dup.jump.prob) #p-value = 2.2e-16


########### Testing OC for all nodes of 2641 SSD trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.SSD.liver.jump<-summary.tree(Calibrated.SSD.liver.jump)
summary.SSD.liver.jump$Exp.abs <- abs(summary.SSD.liver.jump$pic.liver)
summary.SSD.liver.jumps<-summary.SSD.liver.jump[,!(names(summary.SSD.liver.jump) %in% c("pic.liver"))]
nodes.contrast.SSD.liver.jumps.15spe <- summary.SSD.liver.jumps[which(!is.na(summary.SSD.liver.jumps$Exp.abs)),]
nodes.contrast.SSD.liver.jumps.15spe$Event <- factor(nodes.contrast.SSD.liver.jumps.15spe$events, levels=c("speciation", "duplication"))
nodes.contrast.SSD.liver.jumps.15spe<-nodes.contrast.SSD.liver.jumps.15spe[(which(!(nodes.contrast.SSD.liver.jumps.15spe$label %in% c("Clupeocephala","Osteoglossocephalai")))),]
rm(summary.SSD.liver.jump)
rm(summary.SSD.liver.jumps)

############ Testing OC using all 15 species for 2641 SSD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.SSD.liver.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.SSD.liver.jumps.15spe$index<-nodes.contrast.SSD.liver.jumps.15spe$index.tree
test.output.SSD.liver.processed$node<-test.output.SSD.liver.processed$From
nodes.contrast.SSD.liver.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.SSD.liver.jumps.15spe,test.output.SSD.liver.processed,by=c("index","node"))
nodes.contrast.SSD.liver.jumps.15spe.removing.jump.node <- nodes.contrast.SSD.liver.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.SSD.liver.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.SSD.liver.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.SSD.liver.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.SSD.liver.jumps.15spe.removing.jump.node$Exp.abs<-nodes.contrast.SSD.liver.jumps.15spe.removing.jump.node$pic
nodes.contrast.SSD.liver.jumps.15spe.removing.jump.node<-nodes.contrast.SSD.liver.jumps.15spe.removing.jump.node[(which(!(nodes.contrast.SSD.liver.jumps.15spe.removing.jump.node$label %in% c("Clupeocephala","Osteoglossocephalai")))),]

########### Testing OC for all nodes of 841 SSD exp trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.SSD.liver.nojump<-summary.tree(Calibrated.SSD.liver.nojump)
summary.SSD.liver.nojump$Exp.abs <- abs(summary.SSD.liver.nojump$pic)
summary.SSD.liver.nojumps<-summary.SSD.liver.nojump[,!(names(summary.SSD.liver.nojump) %in% c("pic"))]
nodes.contrast.SSD.liver.nojumps.15spe <- summary.SSD.liver.nojumps[which(!is.na(summary.SSD.liver.nojumps$Exp.abs)),]
nodes.contrast.SSD.liver.nojumps.15spe$Event <- factor(nodes.contrast.SSD.liver.nojumps.15spe$events, levels=c("speciation", "duplication"))
nodes.contrast.SSD.liver.nojumps.15spe<-nodes.contrast.SSD.liver.nojumps.15spe[(which(!(nodes.contrast.SSD.liver.nojumps.15spe$label %in% c("Clupeocephala","Osteoglossocephalai")))),] 
rm(summary.SSD.liver.nojump)
rm(summary.SSD.liver.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.SSD.liver.new.fish<-test.output.SSD.liver.processed[which(test.output.SSD.liver.processed$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##4275 obs

## Step 6: Analysis on fish specific clades
test.output.SSD.liver.new.fish.mod<-unique(test.output.SSD.liver.new.fish[c(6,3,7,8,10,2,1)])  #3615 obs


########### Proportion of jump events in teleosts in the 2641 SSD trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.SSD.liver.trees.fish<-NULL
count<-0
info.jump.SSD.liver.trees.fish<-bind_rows(lapply(Calibrated.SSD.liver.jump, tree.data.collection.fish))
info.jump.SSD.liver.trees.fish$spe.num<-info.jump.SSD.liver.trees.fish$internal.events.fish-info.jump.SSD.liver.trees.fish$dup.num
info.jump.SSD.liver.trees.fish<-info.jump.SSD.liver.trees.fish[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.SSD.liver.fish<-
  test.output.SSD.liver.new.fish  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.SSD.liver.fish)<-c("index.old","jump.dup")


test.spe.SSD.liver.fish<-
  test.output.SSD.liver.new.fish %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.SSD.liver.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.SSD.liver.fish<-full_join(test.spe.SSD.liver.fish,test.dup.SSD.liver.fish)
info.jump.trees.SSD.liver.final.fish<-merge(info.jump.SSD.liver.trees.fish,merged.SSD.liver.fish, by = c("index.old"))
info.jump.trees.SSD.liver.final.fish[is.na(info.jump.trees.SSD.liver.final.fish)] <- 0

## Have to consider trees where atleast one fish specific duplication and speciation exists
info.jump.trees.SSD.liver.final1.fish<-info.jump.trees.SSD.liver.final.fish[which(info.jump.trees.SSD.liver.final.fish$dup.num>0 & info.jump.trees.SSD.liver.final.fish$spe.num>0),]

## Proportion of jumps
info.jump.trees.SSD.liver.final1.fish$dup.jump.prob<-(info.jump.trees.SSD.liver.final1.fish$jump.dup/(info.jump.trees.SSD.liver.final1.fish$dup.num*2))
info.jump.trees.SSD.liver.final1.fish$spe.jump.prob<-(info.jump.trees.SSD.liver.final1.fish$jump.spe/(info.jump.trees.SSD.liver.final1.fish$spe.num*2))


median(na.omit(info.jump.trees.SSD.liver.final1.fish$dup.jump.prob)) #0.036
median(na.omit(info.jump.trees.SSD.liver.final1.fish$spe.jump.prob)) #0.125
mean(na.omit(info.jump.trees.SSD.liver.final1.fish$dup.jump.prob)) #0.1973
mean(na.omit(info.jump.trees.SSD.liver.final1.fish$spe.jump.prob)) #0.1537

pval.jump.SSD.liver.fish <- paired.wilcox.one.sided(info.jump.trees.SSD.liver.final1.fish$spe.jump.prob,info.jump.trees.SSD.liver.final1.fish$dup.jump.prob) #p-value = 1.84e-4
########### Testing OC for all fish specific nodes of 2641 SSD trees supporting jump##################
nodes.contrast.SSD.liver.jumps.fish<-nodes.contrast.SSD.liver.jumps.15spe[which(nodes.contrast.SSD.liver.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] #10684 obs

############ Testing OC for teleosts in 2641 SSD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 

nodes.contrast.SSD.liver.jumps.fish.removing.jump.node <- NULL
nodes.contrast.SSD.liver.jumps.fish$index<-nodes.contrast.SSD.liver.jumps.fish$index.tree
test.output.SSD.liver.new.fish$node<-test.output.SSD.liver.new.fish$From
nodes.contrast.SSD.liver.jumps.fish.removing.jump.node <-anti_join(nodes.contrast.SSD.liver.jumps.fish,test.output.SSD.liver.new.fish,by=c("index","node"))
nodes.contrast.SSD.liver.jumps.fish.removing.jump.node <- nodes.contrast.SSD.liver.jumps.fish.removing.jump.node[which(!is.na(nodes.contrast.SSD.liver.jumps.fish.removing.jump.node$pic)),]
nodes.contrast.SSD.liver.jumps.fish.removing.jump.node$Event <- factor(nodes.contrast.SSD.liver.jumps.fish.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.SSD.liver.jumps.fish.removing.jump.node$Exp.abs<-nodes.contrast.SSD.liver.jumps.fish.removing.jump.node$pic
nodes.contrast.SSD.liver.jumps.fish.removing.jump.node<-nodes.contrast.SSD.liver.jumps.fish.removing.jump.node[(which(!(nodes.contrast.SSD.liver.jumps.fish.removing.jump.node$label %in% c("Clupeocephala","Osteoglossocephalai")))),]

########### Testing OC for all fish specific nodes of 841 SSD exp trees that do not support trait jump for teleosts ##################
nodes.contrast.SSD.liver.nojumps.fish<-nodes.contrast.SSD.liver.nojumps.15spe[which(nodes.contrast.SSD.liver.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] #1886 obs

#save.image("norm_exp_15spe_4tips_new.RData")


