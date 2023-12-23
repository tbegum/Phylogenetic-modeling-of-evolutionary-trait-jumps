
### This script is written to analyse the results from Selectome based on lrt
###If positive selection promotes trait evolution through rapid jumps following duplications, we expect to find more positive selections in the branches that also support jumps in their traits after duplication

##### This part is for 15 vertebrates ########

############  Step1:To match index of selectome tree with the trees supporting jump for traits##############
## Reading output from selectome
positive.selection.table1<-read.table(paste0(folder4,"mapped.selectome.results.tsv",sep=""), sep = "\t",header = T) #50632 obs
selectome.output.all<-positive.selection.table1
selectome.output.all$Gene<-sub("genetree.", "", selectome.output.all$tree_file)
selectome.output.all$Gene<-as.integer(selectome.output.all$Gene)
selectome.output.all$From<-selectome.output.all$parent
selectome.output.all$To<-selectome.output.all$node # 50632 obs
selectome.output.all<-selectome.output.all[,c(-4)] ##label is removed as it denotes the child branch name instead of parent node

############  Calculating lr for events using average expression levels as trait##################

####### Group1: Analyzing trusted SSD trees supporting jumps in avg exp levels, we aim to see the differences in LR for branches supporting jumps as well as positive selections ########
##Step2: Analyses on branches supporting both the trait jumps and positive selections

index.info.exp.jump<-NULL
count<-0
index.info.exp.jump<-bind_rows(lapply(Calibrated.SSD.exp,tree.index.collect))

## Matching calibrated jump trees index
test.output.exp.processed.jump.calibrated<-test.output.exp.processed[which(test.output.exp.processed$Gene %in% index.info.exp.jump$Gene),] ## 10197 obs

## Matching with branches supporting jump and positive selections
selectome.output.jump.avg.exp<-merge(selectome.output.all,test.output.exp.processed.jump.calibrated,by=c("Gene","From","To")) ## 997 obs
selectome.output.jump.avg.exp$Event <- factor(selectome.output.jump.avg.exp$node.event, levels=c("speciation", "duplication"))
#selectome.output.avgexp.jump.trees.only.ps<-merge(selectome.output.avgexp.jump.trees.only.ps1,summary.exp.jump,by=c("Gene","From"))

## LRT Analyses
ps.exp.jump.spe<-selectome.output.jump.avg.exp$lrt[selectome.output.jump.avg.exp$Event=="speciation"] #1.74e-07
ps.exp.jump.dup<-selectome.output.jump.avg.exp$lrt[selectome.output.jump.avg.exp$Event=="duplication"] #0.4448
pval.jump.lrt.avg.exp<-two.tailed.wilcox(ps.exp.jump.spe,ps.exp.jump.dup)


## Getting node.events from summary data
summary.exp.jump<-summary.tree(Calibrated.SSD.exp)
summary.exp.jump$Gene<-summary.exp.jump$index.tree
summary.exp.jump$node.event<-summary.exp.jump$events
summary.exp.jump$From<-summary.exp.jump$node
selectome.output.jump.avg.exp1<-merge(selectome.output.jump.avg.exp,summary.exp.jump,by=c("Gene","From"))
selectome.output.jump.avg.exp1$Event <- factor(selectome.output.jump.avg.exp1$node.event, levels=c("speciation", "duplication"))


## PIC Analyses
ps.exp.jump.spe.pic<-abs(selectome.output.jump.avg.exp1$pic_Exp[selectome.output.jump.avg.exp1$Event=="speciation"]) #0.036
ps.exp.jump.dup.pic<-abs(selectome.output.jump.avg.exp1$pic_Exp[selectome.output.jump.avg.exp1$Event=="duplication"]) #0.025
pval.jump.lrt.avg.exp.pic<-two.tailed.wilcox(ps.exp.jump.spe.pic,ps.exp.jump.dup.pic) #2.09 × 10-01


#dodge <- position_dodge(width = 0.41)
# plot.selectome.exp <- ggplot(selectome.output.jump.avg.exp,aes(x=Event, y=lrt, fill=Event)) + 
#   guides(colour = guide_legend(override.aes = list(shape = 16))) +
#   geom_boxplot(width=0.5,position = dodge, outlier.colour=NA,notch = T) +
#   #geom_text(data = med, aes(label = lrt.new),fontface = 2,hjust=0.5,vjust =-1)+
#   xlab( NULL ) +
#   ylab(expression(bold(paste("LRT value")))) +
#   theme_classic()+
#   theme(legend.title=element_blank(),legend.position="none") +
#   theme(axis.text=element_text(size=10,face="bold")) +
#   theme(legend.text=element_text(size=10,face="bold"))+
#   annotate("rect", xmin = 1.1, xmax = 1.8, ymin = 0.18, ymax =0.18, alpha=0.8,colour = "blue")+
#   annotate("text", x = 1.5, y = 0.2, label= pval.jump.lrt.avg.exp, fontface = 4)+
#   ggtitle(expression(bold(paste("LRT test for avg expression levels"))))
#   
#  
########### Group2: Analyzing trusted SSD trees supporting jumps in avg exp levels, we aim to see the differences in LRT for branches not supporting jumps but positive selections #############
##Step2: Considering trusted SSD trees with trait jumps where branches just support positive selections

## Matching trusted SSD trees indices those support jump 
selectome.output.jumptrees.avg.exp<-selectome.output.all[which(selectome.output.all$Gene %in% test.output.exp.processed.jump.calibrated$Gene),] ## 10829 obs

## Finding out branches supporting positive selection but no jump 
selectome.output.avgexp.jump.trees.only.ps1<-anti_join(selectome.output.jumptrees.avg.exp,test.output.exp.processed.jump.calibrated,by=c("Gene","From","To")) ## 9832 obs 

## Getting node.events from summary data
summary.exp.jump<-summary.tree(Calibrated.SSD.exp)
summary.exp.jump$Gene<-summary.exp.jump$index.tree
summary.exp.jump$node.event<-summary.exp.jump$events
summary.exp.jump$From<-summary.exp.jump$node

selectome.output.avgexp.jump.trees.only.ps<-merge(selectome.output.avgexp.jump.trees.only.ps1,summary.exp.jump,by=c("Gene","From"))
selectome.output.avgexp.jump.trees.only.ps$Event <- factor(selectome.output.avgexp.jump.trees.only.ps$node.event, levels=c("speciation", "duplication"))

## Analyses
ps.spe.otherthanjump.exp<-selectome.output.avgexp.jump.trees.only.ps$lrt[selectome.output.avgexp.jump.trees.only.ps$Event=="speciation"] #1.177e-07
ps.dup.otherthanjump.exp<-selectome.output.avgexp.jump.trees.only.ps$lrt[selectome.output.avgexp.jump.trees.only.ps$Event=="duplication"] #4.36e-07
pval.lrt.avg.exp.otherthanjump<-two.tailed.wilcox(ps.spe.otherthanjump.exp,ps.dup.otherthanjump.exp)


###### Group3: Analyzing trees not supporting jumps in avg exp levels, we aim to see the differences in LRT for branches just supporting positive selections #####
##Step3: Considering trees without trait jumps where branches supporting positive selections
index.info.exp.nojump<-NULL
count<-0
index.info.exp.nojump<-bind_rows(lapply(Calibrated.SSD.exp.nojump,tree.index.collect))

## Finding out the calibrated gene index for trees without supporting jump in avg exp levels but support positive selections
selectome.output.calibrated.trees.nojump<-selectome.output.all[which(selectome.output.all$Gene %in% index.info.exp.nojump$Gene),] ## 7452 obs

## Getting node.events from summary data
summary.exp.nojump<-summary.tree(Calibrated.SSD.exp.nojump)
summary.exp.nojump$Gene<-summary.exp.nojump$index.tree
summary.exp.nojump$node.event<-summary.exp.nojump$events
summary.exp.nojump$From<-summary.exp.nojump$node

selectome.output.calibrated.trees.nojump1<-merge(selectome.output.calibrated.trees.nojump,summary.exp.nojump,by=c("Gene","From"))
selectome.output.calibrated.trees.nojump1$Event<-factor(selectome.output.calibrated.trees.nojump1$node.event,levels = c("speciation","duplication"))

## Analyses
ps.exp.nojump.spe<-selectome.output.calibrated.trees.nojump1$lrt[selectome.output.calibrated.trees.nojump1$node.event=="speciation"] #1.41e-07
ps.exp.nojump.dup<-selectome.output.calibrated.trees.nojump1$lrt[selectome.output.calibrated.trees.nojump1$node.event=="duplication"] #0.6411
pval.lrt.avg.exp.nojump<-two.tailed.wilcox(ps.exp.nojump.spe,ps.exp.nojump.dup)
#save.image("norm_exp_15spe_4tips_new.RData")


############  Calculating lrt for events using tau as the trait ##################

##Step2: Analyses on branches supporting both the trait jumps and positive selections for trusted SSD trees

index.info.tau.jump<-NULL
count<-0
index.info.tau.jump<-bind_rows(lapply(Calibrated.SSD.tau,tree.index.collect))

## Matching calibrated jump trees index
test.output.tau.processed.jump.calibrated<-test.output.tau.processed[which(test.output.tau.processed$Gene %in% index.info.tau.jump$Gene),] ## 5021 obs

###### Group1: Analyzing trusted SSD trees supporting jumps in tau, we aim to see the differences in LRT values for branches supporting jumps as well as positive selections #########
##Step2: Analyses on branches supporting both the trait jumps and positive selections

## Matching with branches supporting jump and positive selections
selectome.output.jump.tau<-merge(selectome.output.all,test.output.tau.processed.jump.calibrated,by=c("Gene","From","To")) ## 503 obs
selectome.output.jump.tau$Event <-factor(selectome.output.jump.tau$node.event,levels=c("speciation","duplication"))
#selectome.output.jump.tau<-merge(selectome.output.all,test.output.tau.processed,by=c("Gene","From","To"))


## LRT analyses
ps.tau.jump.spe<-selectome.output.jump.tau$lrt[selectome.output.jump.tau$Event=="speciation"] #2.10e-07
ps.tau.jump.dup<-selectome.output.jump.tau$lrt[selectome.output.jump.tau$Event=="duplication"] #1.317
pval.jump.lrt.tau<-two.tailed.wilcox(ps.tau.jump.spe,ps.tau.jump.dup)


## Getting node.events from summary data
summary.tau.jump<-summary.tree(Calibrated.SSD.tau)
summary.tau.jump$Gene<-summary.tau.jump$index.tree
summary.tau.jump$node.event<-summary.tau.jump$events
summary.tau.jump$From<-summary.tau.jump$node
selectome.output.jump.tau1<-merge(selectome.output.jump.tau,summary.tau.jump,by=c("Gene","From"))

## PIC Analyses
ps.tau.jump.spe.pic<-abs(selectome.output.jump.tau1$pic_Tau[selectome.output.jump.tau1$Event=="speciation"]) #0.0048
ps.tau.jump.dup.pic<-abs(selectome.output.jump.tau1$pic_Tau[selectome.output.jump.tau1$Event=="duplication"]) #0.0042
pval.jump.lrt.tau.pic<-two.tailed.wilcox(ps.tau.jump.spe.pic,ps.tau.jump.dup.pic) #3.78 × 10-01


# plot.selectome.tau <- ggplot(data=selectome.output.jump.tau, aes(x=lrt, group=node.event, fill=node.event)) +
#   guides(colour = guide_legend(override.aes = list(shape = 16))) +
#   geom_density(adjust=1.5, alpha=.8) +
#   #coord_cartesian(xlim=c(0, 15)) +
#   coord_cartesian(ylim=c(0, 1)) +
#   theme_classic()+
#   xlab("lrt")+
#   ylab("density")+
#   theme(legend.title=element_blank(),legend.position=c(0.9,0.9)) +
#   #theme(legend.position="none")+
#   theme(axis.text=element_text(size=9,face="bold")) +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
#   theme(legend.text=element_text(size=10,face="bold")) +
#   scale_fill_manual(values=c("#F8766D","#00BFC4"))+
#   annotate("text", x = 20, y = 0.7, label=pval.jump.lrt.tau, fontface = 4)+
#   ggtitle(expression(bold(paste("LRT test on tau, ",tau))))


######### Group2: Analyzing trees supporting jumps in tau, we aim to see the differences in LRT for branches not supporting jumps but positive selections #########
##Step2: Considering trees with trait jumps where branches just support positive selections

## Matching trees supporting jump 
selectome.output.jumptrees.tau<-selectome.output.all[which(selectome.output.all$Gene %in% test.output.tau.processed.jump.calibrated$Gene),] ## 8148 obs

## Finding out branches supporting positive selection but no jump 
selectome.output.tau.jump.trees.only.ps1<-anti_join(selectome.output.jumptrees.tau,test.output.exp.processed.jump.calibrated,by=c("Gene","From","To")) ## 7536 obs 

## Getting node.events from summary data
summary.tau.jump<-summary.tree(Calibrated.SSD.tau)
summary.tau.jump$Gene<-summary.tau.jump$index.tree
summary.tau.jump$node.event<-summary.tau.jump$events
summary.tau.jump$From<-summary.tau.jump$node

selectome.output.tau.jump.trees.only.ps<-merge(selectome.output.tau.jump.trees.only.ps1,summary.tau.jump,by=c("Gene","From"))
selectome.output.tau.jump.trees.only.ps$Event <- factor(selectome.output.tau.jump.trees.only.ps$node.event, levels=c("speciation", "duplication"))

## Analyses
ps.spe.otherthanjump.tau<-selectome.output.tau.jump.trees.only.ps$lrt[selectome.output.tau.jump.trees.only.ps$Event=="speciation"] #1.29e-07
ps.dup.otherthanjump.tau<-selectome.output.tau.jump.trees.only.ps$lrt[selectome.output.tau.jump.trees.only.ps$Event=="duplication"] #9.06e-07
pval.lrt.tau.otherthanjump<-two.tailed.wilcox(ps.spe.otherthanjump.tau,ps.dup.otherthanjump.tau)

###### Group3: Analyzing trusted SSD trees not supporting jumps in tau, we aim to see the differences in LRT for branches just supporting positive selections #####
##Step3: Considering trees without trait jumps where branches supporting positive selections
index.info.tau.nojump<-NULL
count<-0
index.info.tau.nojump<-bind_rows(lapply(Calibrated.SSD.tau.nojump,tree.index.collect))

## Finding out the calibrated gene index for trees without supporting jump in avg exp levels but support positive selections
selectome.output.calibrated.trees.nojump.tau<-selectome.output.all[which(selectome.output.all$Gene %in% index.info.tau.nojump$Gene),] ## 5551 obs

## Getting node.events from summary data
summary.tau.nojump<-summary.tree(Calibrated.SSD.tau.nojump)
summary.tau.nojump$Gene<-summary.tau.nojump$index.tree
summary.tau.nojump$node.event<-summary.tau.nojump$events
summary.tau.nojump$From<-summary.tau.nojump$node

selectome.output.calibrated.trees.tau.nojump1<-merge(selectome.output.calibrated.trees.nojump.tau,summary.tau.nojump,by=c("Gene","From"))
selectome.output.calibrated.trees.tau.nojump1$Event<-factor(selectome.output.calibrated.trees.tau.nojump1$node.event,levels = c("speciation","duplication"))

## Analyses
ps.tau.nojump.spe<-selectome.output.calibrated.trees.tau.nojump1$lrt[selectome.output.calibrated.trees.tau.nojump1$Event=="speciation"] #1.25e-07
ps.tau.nojump.dup<-selectome.output.calibrated.trees.tau.nojump1$lrt[selectome.output.calibrated.trees.tau.nojump1$Event=="duplication"] #6.72e-07
pval.lrt.tau.nojump<-two.tailed.wilcox(ps.tau.nojump.spe,ps.tau.nojump.dup)

#save.image("norm_exp_15spe_4tips_new.RData")

################### This part is for 5 teleost fishes ###################

#############  Calculating lrt for events using average expression levels as trait ##################

####### Group1: Analyzing trusted SSD trees supporting jumps in avg exp levels, we aim to see the differences in LRT for branches supporting jumps as well as positive selections in teleosts ########
##Step2: Analyses on branches supporting both the trait jumps and positive selections in teleosts
selectome.output.jump.avg.exp.fish<-selectome.output.jump.avg.exp1[which(selectome.output.jump.avg.exp1$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Clupeocephala","Osteoglossocephalai","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

## Analyses
ps.exp.jump.spe.fish<-selectome.output.jump.avg.exp.fish$lrt[selectome.output.jump.avg.exp.fish$Event=="speciation"] #1.4e-07
ps.exp.jump.dup.fish<-selectome.output.jump.avg.exp.fish$lrt[selectome.output.jump.avg.exp.fish$Event=="duplication"] #0.724
pval.jump.lrt.avg.exp.fish<-two.tailed.wilcox(ps.exp.jump.spe.fish,ps.exp.jump.dup.fish) #1.68 × 10-02


## PIC Analyses
ps.exp.jump.spe.pic.fish<-abs(selectome.output.jump.avg.exp.fish$pic_Exp[selectome.output.jump.avg.exp.fish$Event=="speciation"]) #0.036
ps.exp.jump.dup.pic.fish<-abs(selectome.output.jump.avg.exp.fish$pic_Exp[selectome.output.jump.avg.exp.fish$Event=="duplication"]) #0.025
pval.jump.lrt.avg.exp.pic.fish<-two.tailed.wilcox(ps.exp.jump.spe.pic.fish,ps.exp.jump.dup.pic.fish) #6.76 × 10-01


########### Group2: Analyzing trusted SSD trees supporting jumps in avg exp levels, we aim to see the differences in LRT for branches not supporting jumps but positive selections in teleosts #############
##Step2: Considering trusted SSD trees with trait jumps where branches just support positive selections in teleosts

selectome.output.avgexp.jump.trees.only.ps.fish<-selectome.output.avgexp.jump.trees.only.ps[which(selectome.output.avgexp.jump.trees.only.ps$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Clupeocephala","Osteoglossocephalai","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

## Analyses
ps.spe.otherthanjump.exp.fish<-selectome.output.avgexp.jump.trees.only.ps.fish$lrt[selectome.output.avgexp.jump.trees.only.ps.fish$Event=="speciation"] #5.89e-07
ps.dup.otherthanjump.exp.fish<-selectome.output.avgexp.jump.trees.only.ps.fish$lrt[selectome.output.avgexp.jump.trees.only.ps.fish$Event=="duplication"] #0.05
pval.lrt.avg.exp.otherthanjump.fish<-two.tailed.wilcox(ps.spe.otherthanjump.exp.fish,ps.dup.otherthanjump.exp.fish) #4.58 × 10-02

###### Group3: Analyzing trusted SSD trees not supporting jumps in average exp levels, we aim to see the differences in LRT for branches just supporting positive selections in teleosts #####
##Step3: Considering trees without trait jumps where branches supporting positive selections in teleosts

selectome.output.calibrated.trees.exp.nojump1.fish<-selectome.output.calibrated.trees.nojump1[which(selectome.output.calibrated.trees.nojump1$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Clupeocephala","Osteoglossocephalai","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

## Analyses
ps.exp.nojump.spe.fish<-selectome.output.calibrated.trees.exp.nojump1.fish$lrt[selectome.output.calibrated.trees.exp.nojump1.fish$Event=="speciation"] #8.61e-07
ps.exp.nojump.dup.fish<-selectome.output.calibrated.trees.exp.nojump1.fish$lrt[selectome.output.calibrated.trees.exp.nojump1.fish$Event=="duplication"] #0.016
pval.lrt.exp.nojump.fish<-two.tailed.wilcox(ps.exp.nojump.spe.fish,ps.exp.nojump.dup.fish) #3.41 × 10-01

#############  Calculating lrt for events using tau as trait ##################

####### Group1: Analyzing trusted SSD trees supporting jumps in tau, we aim to see the differences in LRT for branches supporting jumps as well as positive selections in teleosts ########
##Step2: Analyses on branches supporting both the trait jumps and positive selections in teleosts
selectome.output.jump.tau.fish<-selectome.output.jump.tau1[which(selectome.output.jump.tau1$clade.name %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Clupeocephala","Osteoglossocephalai","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

## Analyses
ps.tau.jump.spe.fish<-selectome.output.jump.tau.fish$lrt[selectome.output.jump.tau.fish$Event=="speciation"] #3.11e-07
ps.tau.jump.dup.fish<-selectome.output.jump.tau.fish$lrt[selectome.output.jump.tau.fish$Event=="duplication"] #2.16
pval.jump.lrt.tau.fish<-two.tailed.wilcox(ps.tau.jump.spe.fish,ps.tau.jump.dup.fish) #9.58 × 10-03


## PIC Analyses
ps.tau.jump.spe.pic.fish<-abs(selectome.output.jump.tau.fish$pic_Tau[selectome.output.jump.tau.fish$Event=="speciation"]) #0.0036
ps.tau.jump.dup.pic.fish<-abs(selectome.output.jump.tau.fish$pic_Tau[selectome.output.jump.tau.fish$Event=="duplication"]) #0.013
pval.jump.lrt.tau.pic.fish<-two.tailed.wilcox(ps.tau.jump.spe.pic.fish,ps.tau.jump.dup.pic.fish) #6.51 × 10-01


########### Group2: Analyzing trusted SSD trees supporting jumps in tau, we aim to see the differences in LRT for branches not supporting jumps but positive selections in teleosts #############
##Step2: Considering trusted SSD trees with trait jumps where branches just support positive selections in teleosts

selectome.output.tau.jump.trees.only.ps.fish<-selectome.output.tau.jump.trees.only.ps[which(selectome.output.tau.jump.trees.only.ps$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Clupeocephala","Osteoglossocephalai","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

## Analyses
ps.spe.otherthanjump.tau.fish<-selectome.output.tau.jump.trees.only.ps.fish$lrt[selectome.output.tau.jump.trees.only.ps.fish$Event=="speciation"] #3.98e-07
ps.dup.otherthanjump.tau.fish<-selectome.output.tau.jump.trees.only.ps.fish$lrt[selectome.output.tau.jump.trees.only.ps.fish$Event=="duplication"] #8.40e-07
pval.lrt.tau.otherthanjump.fish<-two.tailed.wilcox(ps.spe.otherthanjump.tau.fish,ps.dup.otherthanjump.tau.fish) #2.7 × 10-01

###### Group3: Analyzing trusted SSD trees not supporting jumps in tau, we aim to see the differences in LRT for branches just supporting positive selections in teleosts #####
##Step3: Considering trees without trait jumps where branches supporting positive selections in teleosts

selectome.output.calibrated.trees.tau.nojump1.fish<-selectome.output.calibrated.trees.tau.nojump1[which(selectome.output.calibrated.trees.tau.nojump1$label %in% c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio","Acanthomorphata","Clupeocephala","Osteoglossocephalai","Oryzias.latipes","Characiphysae","Euteleosteomorpha","Characoidei","Oryzias","Atherinomorphae","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus","Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae","Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")) ,]

## Analyses
ps.tau.nojump.spe.fish<-selectome.output.calibrated.trees.tau.nojump1.fish$lrt[selectome.output.calibrated.trees.tau.nojump1.fish$Event=="speciation"] #6.37e-07
ps.tau.nojump.dup.fish<-selectome.output.calibrated.trees.tau.nojump1.fish$lrt[selectome.output.calibrated.trees.tau.nojump1.fish$Event=="duplication"] #0.215
pval.lrt.tau.nojump.fish<-two.tailed.wilcox(ps.tau.nojump.spe.fish,ps.tau.nojump.dup.fish) #7.71 × 10-02

#save.image("norm_exp_15spe_4tips_new.RData")


