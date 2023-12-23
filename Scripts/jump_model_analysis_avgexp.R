

### This code is written to work with "trait jump model" for average expression

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

## Analysis (with 70% posterior probability of trait jump)

######################### Chunk2: reading output of trait jump #############################

test.output.exp<-read.csv(paste0(folder1,"To_pablo/Avg_exp/trait_jump_output_avg_exp.csv", sep=""),sep=",",header=T) ##159980 obs

## Only considering trees with posterior probability of jump more than 70%
test.output1<-test.output.exp[which(test.output.exp$JumpProbability>=0.7),] ##33854 obs


########################## Chunk 3: Analysis on all contrasts standardized trees for traits ###############################

### Step 1: Considering (Calibrated.all.standardized.exp) 4603 tree data passing the diagnostic test for average expression levels 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.all.standardized.exp,tree.index.collect))
rm(count)

test.output.exp.all<-NULL
test.output.exp.all<-merge(index.info,test.output1,by="Gene") ## 17106 observations
index.trees.jump<-unique(test.output.exp.all$tree.num) ##2604 trees
Calibrated.standardized.exp.jump<-Calibrated.all.standardized.exp[c(index.trees.jump)] ##2604/4603=56.57% Brownian trees supported trait jump model for average expression levels
Calibrated.standardized.exp.nojump<-Calibrated.all.standardized.exp[-c(index.trees.jump)] ## 1999 trees with no support for jump in tau
rm(index.info)
rm(index.trees.jump)

### Step2: Calculating PICs of required trait for trees to match the index number of the trees on which levolution was run
Count<-0
trees.all.interest<-lapply(trees.all.interest,contrast.calc,trait="avg.exp")
#trees.all.interest<-lapply(trees.all.interest,contrast.calc,trait="TPM.testis")
Count<-0
Calibrated.standardized.exp.jump<-lapply(Calibrated.standardized.exp.jump,contrast.calc,trait="avg.exp") 
#Calibrated.standardized.exp.jump<-lapply(Calibrated.standardized.exp.jump,contrast.calc,trait="TPM.testis") 
rm(Count)

### Step3: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.exp.all))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.exp.all)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 17106 obs
}

test.output.exp.processed<-bind_cols(test.output.exp.all,processed.fout) ## 17106 obs for all the 15 species
rm(i)
rm(test.output.exp.all)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

## Is OC supported for exp for vertebrate clades supporting jump?
pic.spe.all<-test.output.exp.processed$pic[which(test.output.exp.processed$node.event=="speciation")]
pic.dup.all<-test.output.exp.processed$pic[which(test.output.exp.processed$node.event=="duplication")]
median(pic.spe.all) ##0.06
median(pic.dup.all) ##0.074
Pval.all.jump.all<-two.tailed.wilcox(pic.spe.all,pic.dup.all) # < 2.2e-16

## Step 4: Analysis on all 15 vertebrates clades
test.output.exp.new.mod<-unique(test.output.exp.processed[c(6,3,7,8,10,2,1)])  #17135 obs

########### Proportion of jump events in 15 vertebrates species in the 2604 trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.trees<-NULL
count<-0
info.jump.trees<-bind_rows(lapply(Calibrated.standardized.exp.jump, tree.data.collection))
info.jump.trees$spe.num<-info.jump.trees$internal.events-info.jump.trees$dup.num
info.jump.trees<-info.jump.trees[c(11,12,3:5)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.all<-
  test.output.exp.processed  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.all)<-c("index.old","jump.dup")

test.spe.all<-
  test.output.exp.processed%>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.all)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.all<-full_join(test.spe.all,test.dup.all)
info.jump.trees.final.all<-merge(info.jump.trees,merged.all, by = c("index.old"))
info.jump.trees.final.all[is.na(info.jump.trees.final.all)] <- 0

## Have to consider trees where atleast one duplication and speciation exists
info.jump.trees.final1.all<-info.jump.trees.final.all[which(info.jump.trees.final.all$dup.num>0 & info.jump.trees.final.all$spe.num>0),]
rm(info.jump.trees.final.all)

## Proportion of jumps
info.jump.trees.final1.all$dup.jump.prob<-(info.jump.trees.final1.all$jump.dup/(info.jump.trees.final1.all$dup.num*2))
info.jump.trees.final1.all$spe.jump.prob<-(info.jump.trees.final1.all$jump.spe/(info.jump.trees.final1.all$spe.num*2))

## Stats
median(na.omit(info.jump.trees.final1.all$dup.jump.prob)) #0.2
median(na.omit(info.jump.trees.final1.all$spe.jump.prob)) #0.07
mean(na.omit(info.jump.trees.final1.all$dup.jump.prob)) #0.3059
mean(na.omit(info.jump.trees.final1.all$spe.jump.prob)) #0.1610

pval.jump.all <- paired.wilcox(info.jump.trees.final1.all$spe.jump.prob,info.jump.trees.final1.all$dup.jump.prob) #p-value = 2.2e-16

## For median data
dtfr.all<-NULL
dtfr.all<-data.frame(prop.jump=info.jump.trees.final1.all$spe.jump.prob,events="speciation")
dtfr.all<-rbind(dtfr.all,data.frame(prop.jump=info.jump.trees.final1.all$dup.jump.prob,events="duplication"))
dtfr.all$Exp.abs <-abs(dtfr.all$prop.jump)
dtfr.all$Event<-factor(dtfr.all$events, levels=c("speciation","duplication"))

library(tidyverse)
median.data.jump.exp<-NULL
median.data.jump.exp<- dtfr.all%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.exp$pic.round<-paste0("Median = ",round(median.data.jump.exp$Exp.abs,4), sep="")
median.data.jump.exp$count<- paste0(median.data.jump.exp$Event,"\n\n","(n = ",median.data.jump.exp$freq,")")

plot4A<-boxplot.new2(dtfr.all,pval.jump.all,"Exp",median.data.jump.exp) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.all, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.exp$count) ## check

# Create Data
data.jump.exp <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(2604,1999)
)

# Compute the position of labels
data.jump.exp <- data.jump.exp %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.exp$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.exp$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.exp<-ggplot(data.jump.exp, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
  #scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 2604 trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.exp.jump<-summary.tree(Calibrated.standardized.exp.jump)
summary.exp.jump$Exp.abs <- abs(summary.exp.jump$pic_Exp)
summary.exp.jumps<-summary.exp.jump[,!(names(summary.exp.jump) %in% c("pic_Exp"))]
nodes.contrast.exp.jumps.15spe <- summary.exp.jumps[which(!is.na(summary.exp.jumps$Exp.abs)),]
nodes.contrast.exp.jumps.15spe$Event <- factor(nodes.contrast.exp.jumps.15spe$events, levels=c("speciation", "duplication"))
rm(summary.exp.jump)
rm(summary.exp.jumps)

############ Testing OC using all 15 species for 2604 trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
library(dplyr)
nodes.contrast.exp.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.exp.jumps.15spe$index<-nodes.contrast.exp.jumps.15spe$index.tree
test.output.exp.processed$node<-test.output.exp.processed$From
nodes.contrast.exp.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.exp.jumps.15spe,test.output.exp.processed,by=c("index","node"))
nodes.contrast.exp.jumps.15spe.removing.jump.node <- nodes.contrast.exp.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.exp.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.exp.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.exp.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.exp.jumps.15spe.removing.jump.node$Exp.abs<-nodes.contrast.exp.jumps.15spe.removing.jump.node$pic

pic.spe.removing.jump.15species.2604trees<-speciation.contrast(nodes.contrast.exp.jumps.15spe.removing.jump.node,"pic") #(n=30534)
pic.dup.removing.jump.15species.2604trees<-duplication.contrast(nodes.contrast.exp.jumps.15spe.removing.jump.node,"pic") #(n=3834)


########### Testing OC for all nodes of 1999 trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.exp.nojump<-summary.tree(Calibrated.standardized.exp.nojump)
summary.exp.nojump$Exp.abs <- abs(summary.exp.nojump$pic_Exp)
summary.exp.nojumps<-summary.exp.nojump[,!(names(summary.exp.nojump) %in% c("pic_Exp"))]
nodes.contrast.exp.nojumps.15spe <- summary.exp.nojump[which(!is.na(summary.exp.nojump$Exp.abs)),]
nodes.contrast.exp.nojumps.15spe$Event <- factor(nodes.contrast.exp.nojumps.15spe$events, levels=c("speciation", "duplication"))
rm(summary.exp.nojump)
rm(summary.exp.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.exp.new.fish<-test.output.exp.processed[which(test.output.exp.processed$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##7101 obs

## Is OC supported for exp for fish clades supporting jump?
pic.spe.fish<-test.output.exp.new.fish$pic[which(test.output.exp.new.fish$node.event=="speciation")]
pic.dup.fish<-test.output.exp.new.fish$pic[which(test.output.exp.new.fish$node.event=="duplication")]
median(pic.spe.fish) ##0.064
median(pic.dup.fish) ##0.100
Pval.all.jump.fish<-two.tailed.wilcox(pic.spe.fish,pic.dup.fish) # < 2.21e-16

## Step 6: Analysis on fish specific clades
test.output.exp.new.fish.mod<-unique(test.output.exp.new.fish[c(6,3,7,8,10,2,1)])  #7101 obs
rm(test.output.exp.new.fish)

########### Proportion of jump events in teleosts in the 2604 trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.trees.fish<-NULL
count<-0
info.jump.trees.fish<-bind_rows(lapply(Calibrated.standardized.exp.jump, tree.data.collection.fish))
info.jump.trees.fish$spe.num<-info.jump.trees.fish$internal.events.fish-info.jump.trees.fish$dup.num
info.jump.trees.fish<-info.jump.trees.fish[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.fish<-
  test.output.exp.new.fish  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.fish)<-c("index.old","jump.dup")

test.spe.fish<-
  test.output.exp.new.fish%>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.fish<-full_join(test.spe.fish,test.dup.fish)
info.jump.trees.final.fish<-merge(info.jump.trees.fish,merged.fish, by = c("index.old"))
info.jump.trees.final.fish[is.na(info.jump.trees.final.fish)] <- 0


## Have to consider trees where atleast one fish specific duplication and speciation exists
info.jump.trees.final1.fish<-info.jump.trees.final.fish[which(info.jump.trees.final.fish$dup.num>0 & info.jump.trees.final.fish$spe.num>0),]

## Proportion of jumps
info.jump.trees.final1.fish$dup.jump.prob<-(info.jump.trees.final1.fish$jump.dup/(info.jump.trees.final1.fish$dup.num*2))
info.jump.trees.final1.fish$spe.jump.prob<-(info.jump.trees.final1.fish$jump.spe/(info.jump.trees.final1.fish$spe.num*2))


median(na.omit(info.jump.trees.final1.fish$dup.jump.prob)) #0.1666
median(na.omit(info.jump.trees.final1.fish$spe.jump.prob)) #0.1364
mean(na.omit(info.jump.trees.final1.fish$dup.jump.prob)) #0.29
mean(na.omit(info.jump.trees.final1.fish$spe.jump.prob)) #0.25

pval.jump.fish <- paired.wilcox(info.jump.trees.final1.fish$spe.jump.prob,info.jump.trees.final1.fish$dup.jump.prob) #p-value = 2.24e-04
wilcox.test(info.jump.trees.final1.fish$spe.jump.prob,info.jump.trees.final1.fish$dup.jump.prob,paired = T, alternative = "less")$p.value ## one sided p-value = 1.11e-04


## For median data
dtfr.fish<-NULL
dtfr.fish<-data.frame(prop.jump=info.jump.trees.final1.fish$spe.jump.prob,events="speciation")
dtfr.fish<-rbind(dtfr.fish,data.frame(prop.jump=info.jump.trees.final1.fish$dup.jump.prob,events="duplication"))
dtfr.fish$Exp.abs <-abs(dtfr.fish$prop.jump)
dtfr.fish$Event<-factor(dtfr.fish$events, levels=c("speciation","duplication"))


library(tidyverse)
median.data.jump.exp.fish<-NULL
median.data.jump.exp.fish<- dtfr.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.exp.fish$pic.round<-paste0("Median = ",round(median.data.jump.exp.fish$Exp.abs,4), sep="")
median.data.jump.exp.fish$count<- paste0(median.data.jump.exp.fish$Event,"\n\n","(n = ",median.data.jump.exp.fish$freq,")")

plot4B<-boxplot.new2(dtfr.fish,pval.jump.fish,"Exp",median.data.jump.exp.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 5 teleosts"))))+
  #ylab(expression(bold(paste("Proportions of jump events for teleosts ")))) +
  ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.fish, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.exp.fish$count) ## check

########### Testing OC for all fish specific nodes of 2604 trees supporting jump##################
nodes.contrast.exp.jumps.fish<-nodes.contrast.exp.jumps.15spe[which(nodes.contrast.exp.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 2604 trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
library(dplyr)
nodes.contrast.exp.jumps.fish.removing.jump.node <- NULL
nodes.contrast.exp.jumps.fish$index<-nodes.contrast.exp.jumps.fish$index.tree
test.output.exp.new.fish$node<-test.output.exp.new.fish$From
nodes.contrast.exp.jumps.fish.removing.jump.node <-anti_join(nodes.contrast.exp.jumps.fish,test.output.exp.new.fish,by=c("index","node"))
nodes.contrast.exp.jumps.fish.removing.jump.node <- nodes.contrast.exp.jumps.fish.removing.jump.node[which(!is.na(nodes.contrast.exp.jumps.fish.removing.jump.node$pic)),]
nodes.contrast.exp.jumps.fish.removing.jump.node$Event <- factor(nodes.contrast.exp.jumps.fish.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.exp.jumps.fish.removing.jump.node$Exp.abs<-nodes.contrast.exp.jumps.fish.removing.jump.node$pic

pic.spe.removing.jump.fish.2604trees<-speciation.contrast(nodes.contrast.exp.jumps.fish.removing.jump.node,"pic") #(n=9024)
pic.dup.removing.jump.fish.2604trees<-duplication.contrast(nodes.contrast.exp.jumps.fish.removing.jump.node,"pic") #(n=2185)


########### Testing OC for all fish specific nodes of 1999 trees that do not support trait jump for teleosts ##################
nodes.contrast.exp.nojumps.fish<-nodes.contrast.exp.nojumps.15spe[which(nodes.contrast.exp.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")

################## Chunk 4: Analysis on contrasts standardized sure 3RWGD trees for traits ##########################

### Step 1: Considering (Calibrated.3R.exp) 774 3RWGD trees data passing the diagnostic test for average expression levels 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.3R.exp,tree.index.collect))
rm(count)

test.output.exp.all<-NULL
test.output.exp.all<-merge(index.info,test.output1,by="Gene") ## 2300 observations
index.trees.jump<-unique(test.output.exp.all$tree.num) ##398 trees
Calibrated.3RWGD.exp.jump<-Calibrated.3R.exp[c(index.trees.jump)] ##398/774=51.42% Brownian trees supported trait jump model for average expression levels
Calibrated.3RWGD.exp.nojump<-Calibrated.3R.exp[-c(index.trees.jump)] ## 376 trees with no support for jump in tau
rm(index.info)
rm(index.trees.jump)

### Step2: Calculating PICs of required trait for trees to match the index number of the trees on which levolution was run
Count<-0
trees.all.interest<-lapply(trees.all.interest,contrast.calc,trait="avg.exp")

Count<-0
Calibrated.3RWGD.exp.jump<-lapply(Calibrated.3RWGD.exp.jump,contrast.calc,trait="avg.exp") 
rm(Count)

### Step3: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.exp.all))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.exp.all)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 2300 obs
}

test.output.3RWGD.exp.processed<-bind_cols(test.output.exp.all,processed.fout) ## 2300 obs for all the 15 species
rm(i)
rm(test.output.exp.all)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

##Since our calibrated tree shows both the FishWGD and duplication nodes, 
##we want to modify the node events of the jump output to separate out the effect of FishWGD from the SSDs in those trees.

jump.3RWGD.exptrees.summary<-summary.tree(Calibrated.3RWGD.exp.jump)
jump.3RWGD.exptrees.summary$Gene<-jump.3RWGD.exptrees.summary$index.tree
jump.3RWGD.exptrees.summary$From<-jump.3RWGD.exptrees.summary$node
test.output.3RWGD.exp.processed.mod<-merge(test.output.3RWGD.exp.processed,jump.3RWGD.exptrees.summary,by=c("Gene","From"))
test.output.3RWGD.exp.processed.mod<-test.output.3RWGD.exp.processed.mod[c(1,3,2,4:11,16)]
test.output.3RWGD.exp.processed.mod$node.event<-test.output.3RWGD.exp.processed.mod$events
colnames(test.output.3RWGD.exp.processed.mod)[10]<-"pic"
test.output.3RWGD.exp.processed.mod<-test.output.3RWGD.exp.processed.mod[c(1:11)] ## 2300 obs

## Is OC supported for vertebrate clades of 3RWGD exp trees that are supporting jump?
## Answer is no
pic.spe.3RWGD.exp<-test.output.3RWGD.exp.processed.mod$pic[which(test.output.3RWGD.exp.processed.mod$node.event=="speciation")]
pic.WGD.3RWGD.exp<-test.output.3RWGD.exp.processed.mod$pic[which(test.output.3RWGD.exp.processed.mod$node.event=="FishWGD")]
pic.SSD.3RWGD.exp<-test.output.3RWGD.exp.processed.mod$pic[which(test.output.3RWGD.exp.processed.mod$node.event=="duplication")]

## Step 4: Analysis on all 15 vertebrates clades
test.output.3RWGD.exp.processed.exp.new.mod<-unique(test.output.3RWGD.exp.processed.mod[c(6,3,7,8,10,2,1)])  #1799 obs

########### Proportion of jump events in 15 vertebrates species in the 398 3RWGD trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.3RWGD.trees<-NULL
info.jump.3RWGD.trees<-bind_rows(lapply(Calibrated.3RWGD.exp.jump, tree.data.collection.3RWGD))
info.jump.3RWGD.trees$spe.num<-info.jump.3RWGD.trees$internal.events-(info.jump.3RWGD.trees$dup.num + info.jump.3RWGD.trees$WGD.num)
info.jump.3RWGD.trees<-info.jump.3RWGD.trees[c(11,14,3,12,5,4,13)]


## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.3RWGD.all<-
  test.output.3RWGD.exp.processed.mod %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.3RWGD.all)<-c("index.old","jump.dup")

test.WGD.3RWGD.all<-
  test.output.3RWGD.exp.processed.mod %>%
  filter(node.event=="FishWGD" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.WGD.3RWGD.all)<-c("index.old","jump.WGD")

test.spe.3RWGD.all<-
  test.output.3RWGD.exp.processed.mod%>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.3RWGD.all)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.3RWGD.all1<-full_join(test.spe.3RWGD.all,test.dup.3RWGD.all)
merged.3RWGD.all<-full_join(merged.3RWGD.all1,test.WGD.3RWGD.all)
info.jump.trees.final.3RWGD.all<-merge(info.jump.3RWGD.trees,merged.3RWGD.all, by = c("index.old"))
info.jump.trees.final.3RWGD.all[is.na(info.jump.trees.final.3RWGD.all)] <- 0
rm(merged.3RWGD.all1)

## Have to consider trees where atleast one WGD and one speciation exists
info.jump.trees.final1.3RWGD.all<-info.jump.trees.final.3RWGD.all[which(info.jump.trees.final.3RWGD.all$spe.num>0 & info.jump.trees.final.3RWGD.all$WGD.num>0),]
rm(info.jump.trees.final.3RWGD.all)

## Proportion of jumps
info.jump.trees.final1.3RWGD.all$dup.jump.prob<-(info.jump.trees.final1.3RWGD.all$jump.dup/(info.jump.trees.final1.3RWGD.all$dup.num*2))
info.jump.trees.final1.3RWGD.all$spe.jump.prob<-(info.jump.trees.final1.3RWGD.all$jump.spe/(info.jump.trees.final1.3RWGD.all$spe.num*2))
info.jump.trees.final1.3RWGD.all$WGD.jump.prob<-(info.jump.trees.final1.3RWGD.all$jump.WGD/(info.jump.trees.final1.3RWGD.all$WGD.num*2))

## Since,many WGD trees didnot have a duplication node, it produces NA for dup.jump.prob
## To get rid of this
info.jump.trees.final2.3RWGD.all<-info.jump.trees.final1.3RWGD.all[complete.cases(info.jump.trees.final1.3RWGD.all),]

## Stats
median(info.jump.trees.final2.3RWGD.all$dup.jump.prob) #0.25
median(info.jump.trees.final2.3RWGD.all$spe.jump.prob) #0.045
median(info.jump.trees.final2.3RWGD.all$WGD.jump.prob) #0
mean(info.jump.trees.final2.3RWGD.all$dup.jump.prob) #0.294
mean(info.jump.trees.final2.3RWGD.all$spe.jump.prob) #0.1447
mean(info.jump.trees.final2.3RWGD.all$WGD.jump.prob) #0.1492

pval.jump.3RWGD.exp.all.spe2dup <- paired.wilcox(info.jump.trees.final2.3RWGD.all$spe.jump.prob,info.jump.trees.final2.3RWGD.all$dup.jump.prob) #p-value = 7e-10
pval.jump.3RWGD.exp.all.spe2WGD <- paired.wilcox(info.jump.trees.final2.3RWGD.all$spe.jump.prob,info.jump.trees.final2.3RWGD.all$WGD.jump.prob) #p-value = 2.08e-3

## For median data
dtfr.3RWGD.exp.15spe<-NULL
dtfr.3RWGD.exp.15spe<-data.frame(prop.jump=info.jump.trees.final2.3RWGD.all$spe.jump.prob,events="speciation")
dtfr.3RWGD.exp.15spe<-rbind(dtfr.3RWGD.exp.15spe,data.frame(prop.jump=info.jump.trees.final2.3RWGD.all$dup.jump.prob,events="duplication"))
dtfr.3RWGD.exp.15spe<-rbind(dtfr.3RWGD.exp.15spe,data.frame(prop.jump=info.jump.trees.final2.3RWGD.all$WGD.jump.prob,events="FishWGD"))
dtfr.3RWGD.exp.15spe$Exp.abs <-abs(dtfr.3RWGD.exp.15spe$prop.jump)
dtfr.3RWGD.exp.15spe$Event<-factor(dtfr.3RWGD.exp.15spe$events, levels=c("speciation","duplication","FishWGD"))

library(tidyverse)
median.data.jump.3RWGD.exp.15spe<-NULL
median.data.jump.3RWGD.exp.15spe<- dtfr.3RWGD.exp.15spe%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.3RWGD.exp.15spe$pic.round<-paste0("Median = ",round(median.data.jump.3RWGD.exp.15spe$Exp.abs,4), sep="")
median.data.jump.3RWGD.exp.15spe$count<- paste0(median.data.jump.3RWGD.exp.15spe$Event,"\n\n","(n = ",median.data.jump.3RWGD.exp.15spe$freq,")")

plot3R1<-boxplot.new.3RWGD2(dtfr.3RWGD.exp.15spe,pval.jump.3RWGD.exp.all.spe2dup,pval.jump.3RWGD.exp.all.spe2WGD,"Exp.3RWGD",median.data.jump.3RWGD.exp.15spe) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.05, xmax = 1.7, ymin = 0.67, ymax =0.67, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.4, y = 0.75, label= pval.jump.3RWGD.exp.all.spe2dup, fontface = 4)+
  annotate("rect", xmin = 1.05, xmax = 2.8, ymin = 0.88, ymax =0.88, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.7, y = 0.95, label= pval.jump.3RWGD.exp.all.spe2WGD, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.3RWGD.exp.15spe$count) ## check

# Create Data
data.jump.3RWGD.exp <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(398,376)
)

# Compute the position of labels
data.jump.3RWGD.exp <- data.jump.3RWGD.exp %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.3RWGD.exp$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.3RWGD.exp$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.3RWGD.exp<-ggplot(data.jump.3RWGD.exp, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 398 WGD trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.3RWGD.exp.jump<-summary.tree(Calibrated.3RWGD.exp.jump)
summary.3RWGD.exp.jump$Exp.abs <- abs(summary.3RWGD.exp.jump$pic_Exp)
summary.3RWGD.exp.jumps<-summary.3RWGD.exp.jump[,!(names(summary.3RWGD.exp.jump) %in% c("pic_Exp"))]
nodes.contrast.3RWGD.exp.jumps.15spe <- summary.3RWGD.exp.jumps[which(!is.na(summary.3RWGD.exp.jumps$Exp.abs)),]
nodes.contrast.3RWGD.exp.jumps.15spe$Event <- factor(nodes.contrast.3RWGD.exp.jumps.15spe$events, levels=c("speciation", "duplication", "FishWGD"))
rm(summary.3RWGD.exp.jump)
rm(summary.3RWGD.exp.jumps)

############ Testing OC using all 15 species for 398 WGD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.3RWGD.exp.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.3RWGD.exp.jumps.15spe$index<-nodes.contrast.3RWGD.exp.jumps.15spe$index.tree
test.output.3RWGD.exp.processed.mod$node<-test.output.3RWGD.exp.processed.mod$From
nodes.contrast.3RWGD.exp.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.3RWGD.exp.jumps.15spe,test.output.3RWGD.exp.processed.mod,by=c("index","node"))
nodes.contrast.3RWGD.exp.jumps.15spe.removing.jump.node <- nodes.contrast.3RWGD.exp.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.3RWGD.exp.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.3RWGD.exp.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.3RWGD.exp.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication","FishWGD"))
nodes.contrast.3RWGD.exp.jumps.15spe.removing.jump.node$Exp.abs<-abs(nodes.contrast.3RWGD.exp.jumps.15spe.removing.jump.node$pic)

pic.spe.removing.jump.15species.3RWGD.exp.trees<-speciation.contrast(nodes.contrast.3RWGD.exp.jumps.15spe.removing.jump.node,"pic") 
pic.WGD.removing.jump.15species.3RWGD.exp.trees<-duplication.contrast.3R(nodes.contrast.3RWGD.exp.jumps.15spe.removing.jump.node,"pic")
pic.dup.removing.jump.15species.3RWGD.exp.trees<-duplication.contrast(nodes.contrast.3RWGD.exp.jumps.15spe.removing.jump.node,"pic")

########### Testing OC for all nodes of 376 3RWGD exp trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.3RWGD.exp.nojump<-summary.tree(Calibrated.3RWGD.exp.nojump)
summary.3RWGD.exp.nojump$Exp.abs <- abs(summary.3RWGD.exp.nojump$pic_Exp)
summary.3RWGD.exp.nojumps<-summary.3RWGD.exp.nojump[,!(names(summary.3RWGD.exp.nojump) %in% c("pic_Exp"))]
nodes.contrast.3RWGD.exp.nojumps.15spe <- summary.3RWGD.exp.nojumps[which(!is.na(summary.3RWGD.exp.nojumps$Exp.abs)),]
nodes.contrast.3RWGD.exp.nojumps.15spe$Event <- factor(nodes.contrast.3RWGD.exp.nojumps.15spe$events, levels=c("speciation", "duplication","FishWGD"))
rm(summary.3RWGD.exp.nojump)
rm(summary.3RWGD.exp.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.3RWGD.exp.new.fish<-test.output.3RWGD.exp.processed.mod[which(test.output.3RWGD.exp.processed.mod$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##7101 obs

## Step 6: Analysis on fish specific clades
test.output.3RWGD.exp.processed.new.fish.mod<-unique(test.output.3RWGD.exp.new.fish[c(6,3,7,8,10,2,1)])  #875obs
rm(test.output.3RWGD.exp.new.fish)

########### Proportion of jump events in teleosts in the 398 WGD trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.3RWGD.fish.trees<-NULL
info.jump.3RWGD.fish.trees<-bind_rows(lapply(Calibrated.3RWGD.exp.jump, tree.data.collection.3RWGD.fish))
info.jump.3RWGD.fish.trees$spe.num<-info.jump.3RWGD.fish.trees$internal.events.fish-(info.jump.3RWGD.fish.trees$dup.num + info.jump.3RWGD.fish.trees$WGD.num)
info.jump.3RWGD.fish.trees<-info.jump.3RWGD.fish.trees[c(11,14,3,12,5,4,13)]

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.3RWGD.fish<-
  test.output.3RWGD.exp.new.fish %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.3RWGD.fish)<-c("index.old","jump.dup")

test.WGD.3RWGD.fish<-
  test.output.3RWGD.exp.new.fish %>%
  filter(node.event=="FishWGD" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.WGD.3RWGD.fish)<-c("index.old","jump.WGD")

test.spe.3RWGD.fish<-
  test.output.3RWGD.exp.new.fish%>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.3RWGD.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.3RWGD.fish1<-full_join(test.spe.3RWGD.fish,test.dup.3RWGD.fish)
merged.3RWGD.fish<-full_join(merged.3RWGD.fish1,test.WGD.3RWGD.fish)
info.jump.trees.final.3RWGD.fish<-merge(info.jump.3RWGD.fish.trees,merged.3RWGD.fish, by = c("index.old"))
info.jump.trees.final.3RWGD.fish[is.na(info.jump.trees.final.3RWGD.fish)] <- 0
rm(merged.3RWGD.fish1)

## Have to consider trees where atleast one WGD and one speciation exists
info.jump.trees.final1.3RWGD.fish<-info.jump.trees.final.3RWGD.fish[which(info.jump.trees.final.3RWGD.fish$spe.num>0 & info.jump.trees.final.3RWGD.fish$WGD.num>0),]
rm(info.jump.trees.final.3RWGD.fish)

## Proportion of jumps
info.jump.trees.final1.3RWGD.fish$dup.jump.prob<-(info.jump.trees.final1.3RWGD.fish$jump.dup/(info.jump.trees.final1.3RWGD.fish$dup.num*2))
info.jump.trees.final1.3RWGD.fish$spe.jump.prob<-(info.jump.trees.final1.3RWGD.fish$jump.spe/(info.jump.trees.final1.3RWGD.fish$spe.num*2))
info.jump.trees.final1.3RWGD.fish$WGD.jump.prob<-(info.jump.trees.final1.3RWGD.fish$jump.WGD/(info.jump.trees.final1.3RWGD.fish$WGD.num*2))

## Since,many WGD trees didnot have a duplication node, it produces NA for dup.jump.prob
## To get rid of this
info.jump.trees.final2.3RWGD.fish<-info.jump.trees.final1.3RWGD.fish[complete.cases(info.jump.trees.final1.3RWGD.fish),]

## Stats
median(na.omit(info.jump.trees.final2.3RWGD.fish$dup.jump.prob)) #0.25
median(na.omit(info.jump.trees.final2.3RWGD.fish$spe.jump.prob)) #0.1899
median(na.omit(info.jump.trees.final2.3RWGD.fish$WGD.jump.prob)) #0
mean(na.omit(info.jump.trees.final2.3RWGD.fish$dup.jump.prob)) #0.295
mean(na.omit(info.jump.trees.final2.3RWGD.fish$spe.jump.prob)) #0.278
mean(na.omit(info.jump.trees.final2.3RWGD.fish$WGD.jump.prob)) #0.240

pval.jump.3RWGD.exp.fish.spe2dup <- paired.wilcox(info.jump.trees.final2.3RWGD.fish$spe.jump.prob,info.jump.trees.final2.3RWGD.fish$dup.jump.prob) #p-value = 6.92e-1
pval.jump.3RWGD.exp.fish.spe2WGD <- paired.wilcox(info.jump.trees.final2.3RWGD.fish$spe.jump.prob,info.jump.trees.final2.3RWGD.fish$WGD.jump.prob) #p-value = 1.26e-1

## For median data
dtfr.3RWGD.exp.fish<-NULL
dtfr.3RWGD.exp.fish<-data.frame(prop.jump=info.jump.trees.final2.3RWGD.fish$spe.jump.prob,events="speciation")
dtfr.3RWGD.exp.fish<-rbind(dtfr.3RWGD.exp.fish,data.frame(prop.jump=info.jump.trees.final2.3RWGD.fish$dup.jump.prob,events="duplication"))
dtfr.3RWGD.exp.fish<-rbind(dtfr.3RWGD.exp.fish,data.frame(prop.jump=info.jump.trees.final2.3RWGD.fish$WGD.jump.prob,events="FishWGD"))
dtfr.3RWGD.exp.fish$Exp.abs <-abs(dtfr.3RWGD.exp.fish$prop.jump)
dtfr.3RWGD.exp.fish$Event<-factor(dtfr.3RWGD.exp.fish$events, levels=c("speciation","duplication","FishWGD"))

library(tidyverse)
median.data.jump.3RWGD.exp.fish<-NULL
median.data.jump.3RWGD.exp.fish<- dtfr.3RWGD.exp.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.3RWGD.exp.fish$pic.round<-paste0("Median = ",round(median.data.jump.3RWGD.exp.fish$Exp.abs,4), sep="")
median.data.jump.3RWGD.exp.fish$count<- paste0(median.data.jump.3RWGD.exp.fish$Event,"\n","(n = ",median.data.jump.3RWGD.exp.fish$freq,")")

plot3R2<-boxplot.new.3RWGD2(dtfr.3RWGD.exp.fish,pval.jump.3RWGD.exp.fish.spe2duppval.jump.3RWGD.exp.fish.spe2WGD,"Exp.3RWGD",median.data.jump.3RWGD.exp.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("Jump in exp levels for 5 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.05, xmax = 1.7, ymin = 0.67, ymax =0.67, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.4, y = 0.75, label= pval.jump.3RWGD.exp.fish.spe2dup, fontface = 4)+
  annotate("rect", xmin = 1.05, xmax = 2.8, ymin = 0.88, ymax =0.88, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.7, y = 0.95, label= pval.jump.3RWGD.exp.fish.spe2WGD, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.3RWGD.exp.fish$count) ## check

########### Testing OC for all fish specific nodes of 398 WGD trees supporting jump##################
nodes.contrast.3RWGD.exp.jumps.fish<-nodes.contrast.3RWGD.exp.jumps.15spe[which(nodes.contrast.3RWGD.exp.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 398 WGD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
nodes.contrast.3RWGD.exp.jumps.fish.removing.jump.node<-nodes.contrast.3RWGD.exp.jumps.15spe.removing.jump.node[which(nodes.contrast.3RWGD.exp.jumps.15spe.removing.jump.node$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

########### Testing OC for all fish specific nodes of 376 3RWGD trees that do not support trait jump for teleosts ##################
nodes.contrast.3RWGD.exp.nojumps.fish<-nodes.contrast.3RWGD.exp.nojumps.15spe[which(nodes.contrast.3RWGD.exp.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae","Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")

################## Chunk 5: Analysis on contrasts standardized sure SSD trees for traits ##########################

### Step 1: Considering (Calibrated.SSD.exp) 2717 SSD tree data passing the diagnostic test for average expression levels 
### and matching their index with the output to identify standardized trees with and without trait jump

index.info<-NULL
count<-0
index.info<-bind_rows(lapply(Calibrated.SSD.exp,tree.index.collect))
rm(count)

test.output.exp.all<-NULL
test.output.exp.all<-merge(index.info,test.output1,by="Gene") ## 10197 observations
index.trees.jump<-unique(test.output.exp.all$tree.num) ##1593 trees
Calibrated.SSD.exp.jump<-Calibrated.SSD.exp[c(index.trees.jump)] ##1593/2717=58.63% Brownian trees supported trait jump model for average expression levels
Calibrated.SSD.exp.nojump<-Calibrated.SSD.exp[-c(index.trees.jump)] ## 1124 trees with no support for jump in tau
rm(index.info)
rm(index.trees.jump)

### Step2: Calculating PICs of required trait for trees to match the index number of the trees on which levolution was run
Count<-0
trees.all.interest<-lapply(trees.all.interest,contrast.calc,trait="avg.exp")

Count<-0
Calibrated.SSD.exp.jump<-lapply(Calibrated.SSD.exp.jump,contrast.calc,trait="avg.exp") 
rm(Count)

### Step3: Collecting jump node events info along with the corresponding pic values
Frequency<-0
processed.fout<-NULL
for(i in 1:nrow(test.output.exp.all))
{
  print(i)
  processed.output<-NULL
  processed.output<-process.jump.output2(test.output.exp.all)
  processed.fout<-bind_rows(processed.fout,processed.output) ## 10197 obs
}

test.output.SSD.exp.processed<-bind_cols(test.output.exp.all,processed.fout) ## 10197 obs for all the 15 species
test.output.SSD.exp.processed<-test.output.SSD.exp.processed[(which(!(test.output.SSD.exp.processed$clade.name %in% c("Clupeocephala","Osteoglossocephalai")))),] ##9667 obs
rm(i)
rm(test.output.exp.all)
rm(processed.fout)
rm(processed.output)
rm(Frequency)

## Is OC supported for vertebrate clades of SSD exp trees that are supporting jump?
pic.spe.all<-test.output.SSD.exp.processed$pic[which(test.output.SSD.exp.processed$node.event=="speciation")]
pic.dup.all<-test.output.SSD.exp.processed$pic[which(test.output.SSD.exp.processed$node.event=="duplication")]
median(pic.spe.all) ##0.061
median(pic.dup.all) ##0.088
Pval.all.jump.all<-two.tailed.wilcox(pic.spe.all,pic.dup.all) # < 2.2e-16

## Step 4: Analysis on all 15 vertebrates clades
test.output.SSD.exp.processed.exp.new.mod<-unique(test.output.SSD.exp.processed[c(6,3,7,8,10,2,1)])  #7102 obs

########### Proportion of jump events in 15 vertebrates species in the 1593 SSD trees supporting jump ##################

##Numbers of speciation and duplication events for trees with jumps
info.jump.SSD.trees<-NULL
count<-0
info.jump.SSD.trees<-bind_rows(lapply(Calibrated.SSD.brain.jump, tree.data.collection))
info.jump.SSD.trees$spe.num<-info.jump.SSD.trees$internal.events-info.jump.SSD.trees$dup.num
info.jump.SSD.trees<-info.jump.SSD.trees[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.SSD.all<-
  test.output.SSD.exp.processed %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.SSD.all)<-c("index.old","jump.dup")

test.spe.SSD.all<-
  test.output.SSD.exp.processed%>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.SSD.all)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.SSD.all<-full_join(test.spe.SSD.all,test.dup.SSD.all)
info.jump.trees.final.SSD.all<-merge(info.jump.SSD.trees,merged.SSD.all, by = c("index.old"))
info.jump.trees.final.SSD.all[is.na(info.jump.trees.final.SSD.all)] <- 0

## Have to consider trees where atleast one duplication and speciation exists
info.jump.trees.final1.SSD.all<-info.jump.trees.final.SSD.all[which(info.jump.trees.final.SSD.all$dup.num>0 & info.jump.trees.final.SSD.all$spe.num>0),]
rm(info.jump.trees.final.SSD.all)

## Proportion of jumps
info.jump.trees.final1.SSD.all$dup.jump.prob<-(info.jump.trees.final1.SSD.all$jump.dup/(info.jump.trees.final1.SSD.all$dup.num*2))
info.jump.trees.final1.SSD.all$spe.jump.prob<-(info.jump.trees.final1.SSD.all$jump.spe/(info.jump.trees.final1.SSD.all$spe.num*2))

## Stats
median(na.omit(info.jump.trees.final1.SSD.all$dup.jump.prob)) #0.25
median(na.omit(info.jump.trees.final1.SSD.all$spe.jump.prob)) #0.062
mean(na.omit(info.jump.trees.final1.SSD.all$dup.jump.prob)) #0.355
mean(na.omit(info.jump.trees.final1.SSD.all$spe.jump.prob)) #0.148

pval.jump.SSD.all <- paired.wilcox(info.jump.trees.final1.SSD.all$spe.jump.prob,info.jump.trees.final1.SSD.all$dup.jump.prob) #p-value = 2.2e-16

## For median data
dtfr.SSD.15spe<-NULL
dtfr.SSD.15spe<-data.frame(prop.jump=info.jump.trees.final1.SSD.all$spe.jump.prob,events="speciation")
dtfr.SSD.15spe<-rbind(dtfr.all,data.frame(prop.jump=info.jump.trees.final1.SSD.all$dup.jump.prob,events="duplication"))
dtfr.SSD.15spe$Exp.abs <-abs(dtfr.SSD.15spe$prop.jump)
dtfr.SSD.15spe$Event<-factor(dtfr.SSD.15spe$events, levels=c("speciation","duplication"))

library(tidyverse)
median.data.jump.SSD.exp<-NULL
median.data.jump.SSD.exp<- dtfr.SSD.15spe%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.SSD.exp$pic.round<-paste0("Median = ",round(median.data.jump.SSD.exp$Exp.abs,4), sep="")
median.data.jump.SSD.exp$count<- paste0(median.data.jump.SSD.exp$Event,"\n\n","(n = ",median.data.jump.SSD.exp$freq,")")

plot1<-boxplot.new2(dtfr.SSD.15spe,pval.jump.SSD.all,"Exp",median.data.jump.SSD.exp) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 15 species"))))+
  ylab(expression(bold(paste("Proportions of jump events")))) +
  #ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.SSD.all, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.SSD.exp$count) ## check

# Create Data
data.jump.SSD.exp <- data.frame(
  group=c("Supporting trait jump model","No support for trait jump model"),
  value=c(1593,1124)
)

# Compute the position of labels
data.jump.SSD.exp <- data.jump.SSD.exp %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data.jump.SSD.exp$value) *100) %>%
  mutate(percent=paste0(prop = round((value / sum(data.jump.SSD.exp$value) *100),0),"%","\n(",value," trees)", sep="")) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie.chart.jump.SSD.exp<-ggplot(data.jump.SSD.exp, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position=c(1.2, 0.8)) +
  geom_text(aes(y = ypos, label = percent), color = "black", size=4,fontface = "bold") +
  scale_fill_brewer(palette = "Greens")
#scale_fill_manual(values=c("#c5d9ed", "#facfd2"))

########### Testing OC for all nodes of 1593 SSD trees supporting jump for 15 species ##################

## Summarizing contrasts for normalized exp levels 
summary.SSD.exp.jump<-summary.tree(Calibrated.SSD.exp.jump)
summary.SSD.exp.jump$Exp.abs <- abs(summary.SSD.exp.jump$pic_Exp)
summary.SSD.exp.jumps<-summary.SSD.exp.jump[,!(names(summary.SSD.exp.jump) %in% c("pic_Exp"))]
nodes.contrast.SSD.exp.jumps.15spe <- summary.SSD.exp.jumps[which(!is.na(summary.SSD.exp.jumps$Exp.abs)),]
nodes.contrast.SSD.exp.jumps.15spe$Event <- factor(nodes.contrast.SSD.exp.jumps.15spe$events, levels=c("speciation", "duplication"))
nodes.contrast.SSD.exp.jumps.15spe<-nodes.contrast.SSD.exp.jumps.15spe[(which(!(nodes.contrast.SSD.exp.jumps.15spe$label %in% c("Clupeocephala","Osteoglossocephalai")))),]
rm(summary.SSD.exp.jump)
rm(summary.SSD.exp.jumps)

############ Testing OC using all 15 species for 1593 SSD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
library(dplyr)
nodes.contrast.SSD.exp.jumps.15spe.removing.jump.node <- NULL
nodes.contrast.SSD.exp.jumps.15spe$index<-nodes.contrast.SSD.exp.jumps.15spe$index.tree
test.output.SSD.exp.processed$node<-test.output.SSD.exp.processed$From
nodes.contrast.SSD.exp.jumps.15spe.removing.jump.node <-anti_join(nodes.contrast.SSD.exp.jumps.15spe,test.output.SSD.exp.processed,by=c("index","node"))
nodes.contrast.SSD.exp.jumps.15spe.removing.jump.node <- nodes.contrast.SSD.exp.jumps.15spe.removing.jump.node[which(!is.na(nodes.contrast.SSD.exp.jumps.15spe.removing.jump.node$pic)),]
nodes.contrast.SSD.exp.jumps.15spe.removing.jump.node$Event <- factor(nodes.contrast.SSD.exp.jumps.15spe.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.SSD.exp.jumps.15spe.removing.jump.node$Exp.abs<-nodes.contrast.SSD.exp.jumps.15spe.removing.jump.node$pic
nodes.contrast.SSD.exp.jumps.15spe.removing.jump.node<-nodes.contrast.SSD.exp.jumps.15spe.removing.jump.node[(which(!(nodes.contrast.SSD.exp.jumps.15spe.removing.jump.node$label %in% c("Clupeocephala","Osteoglossocephalai")))),]

pic.spe.removing.jump.15species.SSD.exp.trees<-speciation.contrast(nodes.contrast.SSD.exp.jumps.15spe.removing.jump.node,"pic") #(n=30534)
pic.dup.removing.jump.15species.SSD.exp.trees<-duplication.contrast(nodes.contrast.SSD.exp.jumps.15spe.removing.jump.node,"pic") #(n=3834)


########### Testing OC for all nodes of 1124 SSD exp trees that do not support trait jump for 15 species ##################
## Summarizing contrasts for normalized exp levels 
summary.SSD.exp.nojump<-summary.tree(Calibrated.SSD.exp.nojump)
summary.SSD.exp.nojump$Exp.abs <- abs(summary.SSD.exp.nojump$pic_Exp)
summary.SSD.exp.nojumps<-summary.SSD.exp.nojump[,!(names(summary.SSD.exp.nojump) %in% c("pic_Exp"))]
nodes.contrast.SSD.exp.nojumps.15spe <- summary.SSD.exp.nojumps[which(!is.na(summary.SSD.exp.nojumps$Exp.abs)),]
nodes.contrast.SSD.exp.nojumps.15spe$Event <- factor(nodes.contrast.SSD.exp.nojumps.15spe$events, levels=c("speciation", "duplication"))
nodes.contrast.SSD.exp.nojumps.15spe<-nodes.contrast.SSD.exp.nojumps.15spe[(which(!(nodes.contrast.SSD.exp.nojumps.15spe$label %in% c("Clupeocephala","Osteoglossocephalai")))),] 
rm(summary.SSD.exp.nojump)
rm(summary.SSD.exp.nojumps)

### Step 5: Identifying teleosts fish specific clades
test.output.SSD.exp.new.fish<-test.output.SSD.exp.processed[which(test.output.SSD.exp.processed$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,] ##7101 obs

## Step 6: Analysis on fish specific clades
test.output.SSD.exp.new.fish.mod<-unique(test.output.SSD.exp.new.fish[c(6,3,7,8,10,2,1)])  #2254 obs
rm(test.output.SSD.exp.new.fish)

########### Proportion of jump events in teleosts in the 1593 SSD trees supporting jump ##################

##Numbers of fish specific speciation and duplication events for trees with jumps
info.jump.SSD.exp.trees.fish<-NULL
count<-0
info.jump.SSD.exp.trees.fish<-bind_rows(lapply(Calibrated.SSD.exp.jump, tree.data.collection.fish))
info.jump.SSD.exp.trees.fish$spe.num<-info.jump.SSD.exp.trees.fish$internal.events.fish-info.jump.SSD.exp.trees.fish$dup.num
info.jump.SSD.exp.trees.fish<-info.jump.SSD.exp.trees.fish[c(11,12,3:5)]
rm(count)

## Counting frquency of jump  for "speciation" and "duplication" events
test.dup.SSD.exp.fish<-
  test.output.SSD.exp.new.fish  %>%
  filter(node.event=="duplication" ) %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.dup.SSD.exp.fish)<-c("index.old","jump.dup")


test.spe.SSD.exp.fish<-
  test.output.SSD.exp.new.fish %>%
  filter(node.event=="speciation") %>%
  .$Gene %>% 
  table() %>% 
  as.data.frame()
colnames(test.spe.SSD.exp.fish)<-c("index.old","jump.spe")  

## Merging datasets to get actual statististics of trait jump
merged.SSD.exp.fish<-full_join(test.spe.SSD.exp.fish,test.dup.SSD.exp.fish)
info.jump.trees.SSD.exp.final.fish<-merge(info.jump.SSD.exp.trees.fish,merged.SSD.exp.fish, by = c("index.old"))
info.jump.trees.SSD.exp.final.fish[is.na(info.jump.trees.SSD.exp.final.fish)] <- 0


## Have to consider trees where atleast one fish specific duplication and speciation exists
info.jump.trees.SSD.exp.final1.fish<-info.jump.trees.SSD.exp.final.fish[which(info.jump.trees.SSD.exp.final.fish$dup.num>0 & info.jump.trees.SSD.exp.final.fish$spe.num>0),]

## Proportion of jumps
info.jump.trees.SSD.exp.final1.fish$dup.jump.prob<-(info.jump.trees.SSD.exp.final1.fish$jump.dup/(info.jump.trees.SSD.exp.final1.fish$dup.num*2))
info.jump.trees.SSD.exp.final1.fish$spe.jump.prob<-(info.jump.trees.SSD.exp.final1.fish$jump.spe/(info.jump.trees.SSD.exp.final1.fish$spe.num*2))


median(na.omit(info.jump.trees.SSD.exp.final1.fish$dup.jump.prob)) #0.33
median(na.omit(info.jump.trees.SSD.exp.final1.fish$spe.jump.prob)) #0.125
mean(na.omit(info.jump.trees.SSD.exp.final1.fish$dup.jump.prob)) #0.3572
mean(na.omit(info.jump.trees.SSD.exp.final1.fish$spe.jump.prob)) #0.1936

pval.jump.SSD.exp.fish <- paired.wilcox(info.jump.trees.SSD.exp.final1.fish$spe.jump.prob,info.jump.trees.SSD.exp.final1.fish$dup.jump.prob) #p-value < 2.2e-16

## For median data
dtfr.SSD.exp.fish<-NULL
dtfr.SSD.exp.fish<-data.frame(prop.jump=info.jump.trees.SSD.exp.final1.fish$spe.jump.prob,events="speciation")
dtfr.SSD.exp.fish<-rbind(dtfr.SSD.exp.fish,data.frame(prop.jump=info.jump.trees.SSD.exp.final1.fish$dup.jump.prob,events="duplication"))
dtfr.SSD.exp.fish$Exp.abs <-abs(dtfr.SSD.exp.fish$prop.jump)
dtfr.SSD.exp.fish$Event<-factor(dtfr.SSD.exp.fish$events, levels=c("speciation","duplication"))


library(tidyverse)
median.data.jump.SSD.exp.fish<-NULL
median.data.jump.SSD.exp.fish<- dtfr.SSD.exp.fish%>% 
  group_by(Event) %>% 
  summarise(Mean=mean(prop.jump),
            Exp.abs=median(prop.jump),
            freq = n())
#median <- aggregate(pic ~ Event, Frame, median)
median.data.jump.SSD.exp.fish$pic.round<-paste0("Median = ",round(median.data.jump.SSD.exp.fish$Exp.abs,4), sep="")
median.data.jump.SSD.exp.fish$count<- paste0(median.data.jump.SSD.exp.fish$Event,"\n\n","(n = ",median.data.jump.SSD.exp.fish$freq,")")

plot2<-boxplot.new2(dtfr.SSD.exp.fish,pval.jump.SSD.exp.fish,"Exp",median.data.jump.SSD.exp.fish) +  
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(expression(bold(paste("For 5 teleosts"))))+
  #ylab(expression(bold(paste("Proportions of jump events for teleosts ")))) +
  ylab(" ")+
  theme(plot.title = element_text(color="black", size=12, face="bold"))+
  annotate("rect", xmin = 1.2, xmax = 1.7, ymin = 0.78, ymax =0.78, alpha=0.8,colour = "blue")+
  annotate("text", x = 1.45, y = 0.85, label= pval.jump.SSD.exp.fish, fontface = 4)+
  scale_x_discrete(labels=median.data.jump.SSD.exp.fish$count) ## check

########### Testing OC for all fish specific nodes of 1593 SSD trees supporting jump##################
nodes.contrast.SSD.exp.jumps.fish<-nodes.contrast.SSD.exp.jumps.15spe[which(nodes.contrast.SSD.exp.jumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

############ Testing OC for teleosts in 1593 SSD trees by removing nodes supporting trait jump ###############
##Removing nodes supporting trait jump 
library(dplyr)
nodes.contrast.SSD.exp.jumps.fish.removing.jump.node <- NULL
nodes.contrast.SSD.exp.jumps.fish$index<-nodes.contrast.SSD.exp.jumps.fish$index.tree
test.output.SSD.exp.new.fish$node<-test.output.SSD.exp.new.fish$From
nodes.contrast.SSD.exp.jumps.fish.removing.jump.node <-anti_join(nodes.contrast.SSD.exp.jumps.fish,test.output.SSD.exp.new.fish,by=c("index","node"))
nodes.contrast.SSD.exp.jumps.fish.removing.jump.node <- nodes.contrast.SSD.exp.jumps.fish.removing.jump.node[which(!is.na(nodes.contrast.SSD.exp.jumps.fish.removing.jump.node$pic)),]
nodes.contrast.SSD.exp.jumps.fish.removing.jump.node$Event <- factor(nodes.contrast.SSD.exp.jumps.fish.removing.jump.node$Event, levels=c("speciation", "duplication"))
nodes.contrast.SSD.exp.jumps.fish.removing.jump.node$Exp.abs<-nodes.contrast.SSD.exp.jumps.fish.removing.jump.node$pic
nodes.contrast.SSD.exp.jumps.fish.removing.jump.node<-nodes.contrast.SSD.exp.jumps.fish.removing.jump.node[(which(!(nodes.contrast.SSD.exp.jumps.fish.removing.jump.node$label %in% c("Clupeocephala","Osteoglossocephalai")))),]

########### Testing OC for all fish specific nodes of 1124 SSD exp trees that do not support trait jump for teleosts ##################
nodes.contrast.SSD.exp.nojumps.fish<-nodes.contrast.SSD.exp.nojumps.15spe[which(nodes.contrast.SSD.exp.nojumps.15spe$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

#save.image("norm_exp_15spe_4tips_new.RData")


###################### Working principles of trait jump model #########################

#### Here we aim to identify the duplicates where the jump took place, and the other copy 

## Asymmetric vs. symmetric branches for average expression levels
## Matching indices of trees with dataframe of jumps in expression levels
info.jump.trees$Gene<-info.jump.trees$index.old
df.exp1<-merge(info.jump.trees,test.output.exp.processed,by=c("Gene"))
df.exp1<-df.exp1[c(1,7:8,12,13)]
df.exp1<-df.exp1[order(df.exp1$tree.num,df.exp1$From),]

## Only considering duplicates, we identified the asymmetric and symmetric branches of trait jumps
df.exp2<-df.exp1[which(df.exp1$node.event=="duplication"),] #3486 obs
jump.symmetric<-df.exp2[duplicated(df.exp2),] #938*2=1876 obs
#jump.symmetric<-jump.symmetric[order(jump.symmetric$tree.num,jump.symmetric$From),]
jump.asymmetric<-anti_join(df.exp2,jump.symmetric,by=c("tree.num","From")) #1610 obs

## Only considering species-specific duplicates in vertebrates
#jump.asymmetric.teleosts<-jump.asymmetric[which(jump.asymmetric$clade.name %in% c("Astyanax.mexicanus","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Esox.lucius", "Danio.rerio")),]

## To obtain the specific branch of jump
jump.asymmetric.new<-merge(jump.asymmetric,test.output.exp.processed,by=c("tree.num","Gene","From","node.event","clade.name")) #455 obs
jump.asymmetric.new<-jump.asymmetric.new[c(1,2,3,6,4,5,7)]
jump.symmetric.new<-merge(jump.symmetric,test.output.exp.processed,by=c("tree.num","Gene","From","node.event","clade.name"))
jump.symmetric.new<-jump.symmetric.new[c(1,2,3,6,4,5,7)]

#rm(jump.asymmetric.teleosts)
rm(df.exp1)
rm(df.exp2)

################# For asymmetric trait jump ###############
### First we check how many number of duplicates are there for the nodes,where one of them has experienced a jump in their expression levels 
df.spe.specific.duplicates<-data.frame(tree=NA,duplicate.no=NA)
for(i in 1:nrow(jump.asymmetric.new))
{
  nd<-NULL
  duplicates<-NULL
  tree.ind<-jump.asymmetric.new$tree.num[i]
  duplicates<-tips(Calibrated.all.standardized.exp[[tree.ind]]@phylo,jump.asymmetric.new$From[i])
  df.spe.specific.duplicates<-rbind(df.spe.specific.duplicates,data.frame(tree=tree.ind,duplicate.no=length(duplicates)))
  
}
df.spe.specific.duplicates<-df.spe.specific.duplicates[-1,] ##1610 obs
length(df.spe.specific.duplicates$tree[which(df.spe.specific.duplicates$duplicate.no==2)]) ##995/1610=61.8%
length(df.spe.specific.duplicates$tree[which(df.spe.specific.duplicates$duplicate.no==3)]) ##200/1610=12.42%
length(df.spe.specific.duplicates$tree[which(df.spe.specific.duplicates$duplicate.no>3)]) ##18/1610=25.77%


### Analysing two duplicates

df.specific.duplicate.data<-data.frame(tree=NA,jump.copy=NA,other.copy=NA,Tau.jump=NA,Tau.other=NA,Exp.jump=NA,Exp.other=NA,Max.exp.jump.tissue=NA,Max.exp.other.tissue=NA,label=NA)
for(i in 1:nrow(jump.asymmetric.new))
{
  nd<-NULL
  duplicates<-NULL
  label<-NULL
  tree.ind<-jump.asymmetric.new$tree.num[i]
  duplicates<-tips(Calibrated.all.standardized.exp[[tree.ind]]@phylo,jump.asymmetric.new$From[i])
  outgroups<-ancestor(Calibrated.all.standardized.exp[[tree.ind]]@phylo,jump.asymmetric.new$From[i])
  if(length(duplicates)==2)
  {
    nd<-as.numeric(jump.asymmetric.new$To[i]) ##Node jump
    jump.copy<-Calibrated.all.standardized.exp[[tree.ind]]@data$label[Calibrated.all.standardized.exp[[tree.ind]]@data$node[nd]] ##Duplicate with jump
    other.copy<-duplicates[!(duplicates %in% jump.copy)] ##Duplicate without jump
    nd.other<-as.numeric(Calibrated.all.standardized.exp[[tree.ind]]@data$node[which(Calibrated.all.standardized.exp[[tree.ind]]@data$label %in% other.copy)]) ##Node duplicate without jump
    Tau.jump<-Calibrated.all.standardized.exp[[tree.ind]]@data$Tau[which(Calibrated.all.standardized.exp[[tree.ind]]@data$label==jump.copy)] ###Tau of duplicate with jump
    Tau.other<-Calibrated.all.standardized.exp[[tree.ind]]@data$Tau[which(Calibrated.all.standardized.exp[[tree.ind]]@data$label==other.copy)] ###Tau of other copy
    Exp.jump<-Calibrated.all.standardized.exp[[tree.ind]]@data$avg.exp[which(Calibrated.all.standardized.exp[[tree.ind]]@data$label==jump.copy)]###Exp level of duplicate with jump
    Exp.other<-Calibrated.all.standardized.exp[[tree.ind]]@data$avg.exp[which(Calibrated.all.standardized.exp[[tree.ind]]@data$label==other.copy)]###Exp level of other copy
    Max.Exp.jump<-Calibrated.all.standardized.exp[[tree.ind]]@data$max.exp[which(Calibrated.all.standardized.exp[[tree.ind]]@data$label==jump.copy)]###Max Exp level of duplicate with jump
    tissue.maxexp.jump1<-colnames(Calibrated.all.standardized.exp[[tree.ind]]@data[,10:15][which(Calibrated.all.standardized.exp[[tree.ind]]@data[nd,][,10:15]==Max.Exp.jump)]) ###Tissue of Max Exp level of duplicate with jump
    tissue.maxexp.jump<-unlist(strsplit(tissue.maxexp.jump1,"TPM."))[2]
    Max.Exp.other<-Calibrated.all.standardized.exp[[tree.ind]]@data$max.exp[which(Calibrated.all.standardized.exp[[tree.ind]]@data$label==other.copy)] ###Max Exp level of other copy
    tissue.maxexp.other1<-colnames(Calibrated.all.standardized.exp[[tree.ind]]@data[,10:15][which(Calibrated.all.standardized.exp[[tree.ind]]@data[nd.other,][,10:15]==Max.Exp.other)]) ###Tissue of Max Exp level of other copy
    tissue.maxexp.other<-unlist(strsplit(tissue.maxexp.other1,"TPM."))[2]
    label<-ifelse(tissue.maxexp.jump %in% tissue.maxexp.other,"same","other" )
    df.specific.duplicate.data<-rbind(df.specific.duplicate.data,data.frame(tree=tree.ind,
                                                                                            jump.copy=jump.copy,
                                                                                            other.copy=other.copy,
                                                                                            Tau.jump=Tau.jump,
                                                                                            Tau.other=Tau.other,
                                                                                            Exp.jump=Exp.jump,
                                                                                            Exp.other=Exp.other,
                                                                                            Max.exp.jump.tissue=tissue.maxexp.jump,
                                                                                            Max.exp.other.tissue=tissue.maxexp.other,
                                                                                            label=label))
  }  
}
df.specific.duplicate.data<-df.specific.duplicate.data[-1,] ##995 obs

### Analysing output
nrow(df.specific.duplicate.data[which(df.specific.duplicate.data$Tau.jump >= 0.8 & df.specific.duplicate.data$Tau.other >= 0.8),])#101 cases (10.15%) where both copies are tissue-specific
nrow(df.specific.duplicate.data[which(df.specific.duplicate.data$Tau.jump >= 0.8 & df.specific.duplicate.data$Tau.other <= 0.3),])#131 cases (21.3%) where the jumped copy are tissue-specific but the other is not
nrow(df.specific.duplicate.data[which(df.specific.duplicate.data$Tau.jump <= 0.4 & df.specific.duplicate.data$Tau.other >= 0.8),])#35 cases (14.14%) where the jumped copy are non-tissue-specific but the other is tissue-specific
nrow(df.specific.duplicate.data[which(df.specific.duplicate.data$Tau.jump <= 0.3 & df.specific.duplicate.data$Tau.other <= 0.3),])#52 cases (16.87%) where both copies are non-tissue-specific
##This means that (177+101)/403=69% cases where the jumped copy is tissue-specific

## When both the copies are tissue specific
TF.both<-df.specific.duplicate.data[which(df.specific.duplicate.data$Tau.jump >= 0.8 & df.specific.duplicate.data$Tau.other >= 0.8),]

length(TF.both$label[which(TF.both$label %in% "same")]) #59 obs
length(TF.both$label[which(TF.both$label %in% "other")]) #42 obs






################# For symmetric trait jump ###############
### First we check how many number of duplicates are there for the nodes,where one of them has experienced a jump in their expression levels 
df.spe.specific.duplicates.sym<-data.frame(tree=NA,duplicate.no=NA)
for(i in 1:nrow(jump.symmetric.new))
{
  nd<-NULL
  duplicates<-NULL
  tree.ind<-jump.symmetric.new$tree.num[i]
  duplicates<-tips(Calibrated.all.standardized.exp[[tree.ind]]@phylo,jump.symmetric.new$From[i])
  df.spe.specific.duplicates.sym<-rbind(df.spe.specific.duplicates.sym,data.frame(tree=tree.ind,duplicate.no=length(duplicates)))
  
}
df.spe.specific.duplicates.sym<-df.spe.specific.duplicates.sym[-1,] ##1876 obs
length(df.spe.specific.duplicates.sym$tree[which(df.spe.specific.duplicates.sym$duplicate.no==2)]) ##500/1876=26.65%
length(df.spe.specific.duplicates.sym$tree[which(df.spe.specific.duplicates.sym$duplicate.no==3)]) ##30/1876=1.6%
length(df.spe.specific.duplicates.sym$tree[which(df.spe.specific.duplicates.sym$duplicate.no>3)]) ##1346/1876=71.75%


### Analysing two duplicates

df.specific.duplicate.data.sym<-data.frame(tree=NA,jump.copy=NA,other.copy=NA,Tau.jump=NA,Tau.other=NA,Exp.jump=NA,Exp.other=NA,Max.exp.jump.tissue=NA,Max.exp.other.tissue=NA,label=NA)
for(i in 1:nrow(jump.symmetric.new))
{
  nd<-NULL
  duplicates<-NULL
  label<-NULL
  tree.ind<-jump.symmetric.new$tree.num[i]
  duplicates<-tips(Calibrated.all.standardized.exp[[tree.ind]]@phylo,jump.symmetric.new$From[i])
  if(length(duplicates)==2)
  {
    nd<-as.numeric(jump.symmetric.new$To[i]) ##Node jump
    jump.copy<-Calibrated.all.standardized.exp[[tree.ind]]@data$label[Calibrated.all.standardized.exp[[tree.ind]]@data$node[nd]] ##Duplicate with jump
    other.copy<-duplicates[!(duplicates %in% jump.copy)] ##Duplicate without jump
    nd.other<-as.numeric(Calibrated.all.standardized.exp[[tree.ind]]@data$node[which(Calibrated.all.standardized.exp[[tree.ind]]@data$label %in% other.copy)]) ##Node duplicate without jump
    Tau.jump<-Calibrated.all.standardized.exp[[tree.ind]]@data$Tau[which(Calibrated.all.standardized.exp[[tree.ind]]@data$label==jump.copy)] ###Tau of duplicate with jump
    Tau.other<-Calibrated.all.standardized.exp[[tree.ind]]@data$Tau[which(Calibrated.all.standardized.exp[[tree.ind]]@data$label==other.copy)] ###Tau of other copy
    Exp.jump<-Calibrated.all.standardized.exp[[tree.ind]]@data$avg.exp[which(Calibrated.all.standardized.exp[[tree.ind]]@data$label==jump.copy)]###Exp level of duplicate with jump
    Exp.other<-Calibrated.all.standardized.exp[[tree.ind]]@data$avg.exp[which(Calibrated.all.standardized.exp[[tree.ind]]@data$label==other.copy)]###Exp level of other copy
    Max.Exp.jump<-Calibrated.all.standardized.exp[[tree.ind]]@data$max.exp[which(Calibrated.all.standardized.exp[[tree.ind]]@data$label==jump.copy)]###Max Exp level of duplicate with jump
    tissue.maxexp.jump1<-colnames(Calibrated.all.standardized.exp[[tree.ind]]@data[,10:15][which(Calibrated.all.standardized.exp[[tree.ind]]@data[nd,][,10:15]==Max.Exp.jump)]) ###Tissue of Max Exp level of duplicate with jump
    tissue.maxexp.jump<-unlist(strsplit(tissue.maxexp.jump1,"TPM."))[2]
    Max.Exp.other<-Calibrated.all.standardized.exp[[tree.ind]]@data$max.exp[which(Calibrated.all.standardized.exp[[tree.ind]]@data$label==other.copy)] ###Max Exp level of other copy
    tissue.maxexp.other1<-colnames(Calibrated.all.standardized.exp[[tree.ind]]@data[,10:15][which(Calibrated.all.standardized.exp[[tree.ind]]@data[nd.other,][,10:15]==Max.Exp.other)]) ###Tissue of Max Exp level of other copy
    tissue.maxexp.other<-unlist(strsplit(tissue.maxexp.other1,"TPM."))[2]
    label<-ifelse(tissue.maxexp.jump %in% tissue.maxexp.other,"same","other" )
    df.specific.duplicate.data.sym<-rbind(df.specific.duplicate.data.sym,data.frame(tree=tree.ind,
                                                                            jump.copy=jump.copy,
                                                                            other.copy=other.copy,
                                                                            Tau.jump=Tau.jump,
                                                                            Tau.other=Tau.other,
                                                                            Exp.jump=Exp.jump,
                                                                            Exp.other=Exp.other,
                                                                            Max.exp.jump.tissue=tissue.maxexp.jump,
                                                                            Max.exp.other.tissue=tissue.maxexp.other,
                                                                            label=label))
  }  
}
df.specific.duplicate.data.sym<-df.specific.duplicate.data.sym[-1,] ##500 obs

### Analysing output
nrow(df.specific.duplicate.data.sym[which(df.specific.duplicate.data.sym$Tau.jump >= 0.8 & df.specific.duplicate.data.sym$Tau.other >= 0.8),])#68/500 cases (10.15%) where both copies are tissue-specific
nrow(df.specific.duplicate.data.sym[which(df.specific.duplicate.data.sym$Tau.jump >= 0.8 & df.specific.duplicate.data.sym$Tau.other <= 0.3),])#25 cases (21.3%) where the jumped copy are tissue-specific but the other is not
nrow(df.specific.duplicate.data.sym[which(df.specific.duplicate.data.sym$Tau.jump <= 0.4 & df.specific.duplicate.data.sym$Tau.other >= 0.8),])#33 cases (14.14%) where the jumped copy are non-tissue-specific but the other is tissue-specific
nrow(df.specific.duplicate.data.sym[which(df.specific.duplicate.data.sym$Tau.jump <= 0.3 & df.specific.duplicate.data.sym$Tau.other <= 0.3),])#20 cases (16.87%) where both copies are non-tissue-specific
##This means that (177+101)/403=69% cases where the jumped copy is tissue-specific

## When both the copies are tissue specific
TF.both.sym<-df.specific.duplicate.data.sym[which(df.specific.duplicate.data.sym$Tau.jump >= 0.8 & df.specific.duplicate.data.sym$Tau.other >= 0.8),]

length(TF.both.sym$label[which(TF.both.sym$label %in% "same")]) #42 obs
length(TF.both.sym$label[which(TF.both.sym$label %in% "other")]) #26 obs
