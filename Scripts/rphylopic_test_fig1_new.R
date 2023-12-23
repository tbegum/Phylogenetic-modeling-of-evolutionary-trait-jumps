
## This code is written to use silhouette of rphylopic images
## We modified the previous code as many functions of "rphylopic" is no longer exists

## Setting working directory
setwd("~/Desktop/nwork/Empirical/Manuscript_3/Codes_Manuscipt3/")
#install.packages(c("devtools","rphylopic","ggplot2","ggtree","phytools","deeptime"))

## Required packages
library(devtools)
library(rphylopic)
library(ggplot2)
library(ggtree)
library(phytools)
library(deeptime)

##Drawing the phylogeny
##loading previously stored species tree
load("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/norm_exp_15spe_4tips.rda")

## Reading one tree with only speciation events 
tr<-trees.all.speciation.new[[1]]@phylo
dd<-trees.all.speciation.new[[1]]@data
tr$tip.label<-trees.all.speciation.new[[1]]@data$S[1:length(trees.all.speciation.new[[1]]@phylo$tip.label)]
tr$tip.label[3]<-"Oryzias.latipes"
tr$tip.label[13]<-"Mustela.putorius"
tr$tip.label<-gsub("\\.", " ",tr$tip.label)

## Creating dataframe
d <- ggimage::phylopic_uid(tr$tip.label)
d1<-d
d1$color<-c(10,10,10,10,10,4,4,4,4,4,4,4,4,4,4)
d1$color1<-c("gr1","gr2","gr3","gr4","gr5","gr6","gr7","gr8","gr9","gr10","gr11","gr12","gr13","gr14","gr15")
d1$color2<-c("#000000","#000000","#000000","#000000","#000000","#00868B","#00868B","#00868B","#00868B","#00868B","#00868B","#00868B","#00868B","#00868B","#00868B")
#d1$color2<-c("#e6ccff","#e6ccff","#e6ccff","#e6ccff","#e6ccff","#CCCCFF","#CCCCFF","#CCCCFF","#CCCCFF","#CCCCFF","#CCCCFF","#CCCCFF","#CCCCFF","#CCCCFF","#CCCCFF")
## To maintain the uid of phylopic 


## We used rphylopic v.1.1.1 (https://CRAN.R-project.org/package=rphylopic) to get silhoutte images

## For pike
pike_id_all <- get_uuid(name = "Esox lucius",n=3)  # list images
pike_id <- get_phylopic(uuid = "f8369dec-bdf6-432b-a0c4-41ee5d75286d",height = 256)  # get individual image id
d1$uid[1]<-pike_id_all[2]

## For tilapia
tilapia_id_all<-get_uuid(name = "Oreochromis niloticus", n=2) # list images
tilapia_id <- get_phylopic(uuid = "84c7e672-2593-44a6-a807-cffbd3156cc5",height = 256)  # get individual image id
d1$uid[2]<-tilapia_id_all[1]

## For medaka
medaka_id_all<-get_uuid(name = "Oryzias latipes", n=1)
medaka_id <- get_phylopic(uuid = "b7ae3bbf-c0a1-43bf-80ac-babf90554073",height = 256)  # get individual image id
d1$uid[3]<- medaka_id_all[1]

## For zebrafish
zfish_id_all <- get_uuid(name = "Danio rerio", n=3) # list images
zfish_id <- get_phylopic(uuid = "b7ae3bbf-c0a1-43bf-80ac-babf90554073",height = 256)  # get individual image id
d1$uid[4]<-zfish_id_all[2]

## For cavefish
cavefish_id_all <- get_uuid(name = "Astyanax mexicanus", n=2) # list images
cavefish_id <- get_phylopic(uuid =  "56971d2a-799c-4993-9ffc-7455417efd90",height = 256)  # get individual image id
d1$uid[5]<-cavefish_id_all[1]

## For chicken
chicken_id_all <- get_uuid(name = "Gallus gallus", n=3) # list images
chicken_id <- get_phylopic(uuid = "2de1c95c-7e1f-429b-9c08-17f0a27d176f",height = 256)  # get individual image id
d1$uid[6]<-chicken_id_all[1]

## For monodelphis
Monodelphis_id_all <-get_uuid(name =  "Monodelphis domestica", n=1) # list images
Monodelphis_id <- get_phylopic(uuid =  "dba15a95-cccd-46d5-9efe-9529f51857c0",height = 256)  # get individual image id
d1$uid[7]<-Monodelphis_id_all[1]

## For monodelphis
Macac_id_all <-get_uuid(name =  "Macaca mulatta", n=1) # list images
Macac_id <- get_phylopic(uuid =  "eedde61f-3402-4f7c-9350-49b74f5e1dba",height = 256)  # get individual image id
d1$uid[8]<-Macac_id_all[1]

## For human
human_id_all <- get_uuid(name = "Homo sapiens", n = 20)# list images
human_id<-get_phylopic(uuid =  "c089caae-43ef-4e4e-bf26-973dd4cb65c5",height = 256)
#save_phylopic(img = human_id, path = "/Users/admin/Desktop/nwork/Empirical/Manuscript_3/phylopic_image/human1.png",
#              width = 500, height = 500)
d1$uid[9]<-human_id_all[12]

## For rabbit
rabbit_id_all <- get_uuid(name ="Oryctolagus cuniculus",  n = 2)# list images
rabbit_id <- get_phylopic(uuid =  "d71fcd12-90f6-40c4-b918-0dac0ba3e809",height = 256)  # get individual image id
d1$uid[10]<-rabbit_id_all[2]

## For mouse
mouse_id_all <- get_uuid(name = "Mus musculus", n = 4)# list images
mouse_id <- get_phylopic(uuid =  "c4572239-3b7c-4259-8049-3f44e6d52e6f",height = 256)  # get individual image id
d1$uid[11]<-mouse_id_all[4]

## For rat
rat_id_all <- get_uuid(name =  "Rattus norvegicus", n = 4)# list images
rat_id <- get_phylopic(uuid =  "53f8f372-adc0-4247-bc66-0075e241fc95",height = 256)  # get individual image id
d1$uid[12]<-rat_id_all[2]

## For Mustela (## Since we did not find the rphylopic image of "Mustela putorius furo", we modified the name to ""Mustela putorius")
mustela_id_all <- get_uuid(name = "Mustela putorius", n = 4)# list images
mustela_id <- get_phylopic(uuid =  "0070ddbf-fdcd-4a7b-97c3-4670504dc06f",height = 256)  # get individual image id
d1$uid[13]<-mustela_id_all[1]
#d1$name[13]<-"Mustela putorius"
#rownames(d1)[13]<-"Mustela putorius furo"

## For dog
dog_id_all <- get_uuid(name = "Canis lupus familiaris", n = 14) # list images
dog_id <- get_phylopic(uuid =  "6f3ebbc6-be53-4216-b45b-946f7984669b",height = 256)  # get individual image id
d1$uid[14]<-dog_id_all[5]

## For cow
cow_id_all <- get_uuid(name = "Bos taurus", n = 12) # list images
cow_id <- get_phylopic(uuid =  "dc5c561e-e030-444d-ba22-3d427b60e58a",height = 256)  # get individual image id
d1$uid[15]<-cow_id_all[1]

#d1$img<-c("pike_id","tilapia_id","medaka_id","zfish_id","cavefish_id","chicken_id","Monodelphis_id","Macac_id","human_id","rabbit_id","mouse_id","rat_id","mustela_id", "dog_id","cow_id")

##Silhoutte image for teleosts
fish <- get_uuid(name = "teleostei", n = 500) # list images 
fish_img <- get_phylopic(uuid =  "0b4e39a1-3bca-4b58-8805-ebaae22cd860")  # get individual image id
fish_img1<-recolor_phylopic(img = fish_img, color= "darkgrey",remove_background = T)
save_phylopic(img = fish_img1, path = "/Users/admin/Desktop/nwork/Empirical/Manuscript_3/phylopic_image/teleost.png",
                           width = 500, height = 500)
fish_id <- fish[2]

##Silhoutte image for tetrapods
tetrapods <- get_uuid(name ="tetrapoda", n = 500) # list images 
tetrapod_img <- get_phylopic(uuid =  "24c49c6a-0052-4004-b571-b0ba2cb1d073")  # get individual image id
save_phylopic(img = tetrapod_img, path = "/Users/admin/Desktop/nwork/Empirical/Manuscript_3/phylopic_image/tetrapods.png",
                                         width = 500, height = 500)
tetrapod_id <- tetrapods[15]

## To get the coordinate of nodes in a phylogeny for inserting rphylopic images with backgrounds
xx<-ggtree::fortify(tr)
trr<-ggtree(xx,branch.length="none",layout="fan", size=0.7) 
coord<-ggplot2::ggplot_build(trr)[["data"]][[2]]

## Remodifying the name of "Mustela putorius" to the original one
tr$tip.label[13]<-"Mustela putorius furo"
d1$name[13]<-"Mustela putorius furo"
#save.image("species_uid.rda")

### Plot
#load("species_uid.rda")
setwd("/Users/admin/Desktop/nwork/Empirical/Manuscript_3/Plot_Manuscipt3/")
#png("Fig1_new_rphylo.png",width=400,height=400)
png("Fig1_rphylo_latest3_300dpi",width=400,height=400)
##After setting resolution at 300 dpi, the figure is saved as "Fig1_rphylo_latest.png", "Fig1_rphylo_latest2.png" and "Fig1_rphylo_latest3.png" with different sizes

    ggtree(tr,branch.length="none",layout="fan", size=1) %<+% d1 + 
    #revts(ggtree(tr,branch.length="none",layout="fan",size=1)) %<+% d1 + 
    geom_hilight(node=17, fill="#e4d192")+
    geom_hilight(node=21, fill="#f5efe6") +
    geom_phylopic(aes(uuid=uid,colour=color),size=0.3,alpha=0.9)+
    geom_tiplab(aes(label=label),color='black', angle=0,hjust=0.5, align=T,size=3,fontface=4,vjust=4) +
    theme(legend.position = "none")+
    theme(plot.margin=margin(10,10,10,10))+
    #scale_x_continuous(breaks = seq(‐500, 0, 100),
    #                   labels = seq(500, 0, ‐100)) +
    #coord_geo_polar(dat = "periods") +
    add_phylopic(img=fish_img,x = 4.0,y = 2.62,alpha = 0.6, ysize = 0.25,color = "#47B5FF")+
    add_phylopic(img=tetrapod_img,x = 1,y = 7.42,alpha = 0.7, ysize = 0.15,color = "black")+
    #geom_nodelab(aes(image="../phylopic_image/teleost.png", subset= node == 17), geom="image", align=0,size=0.25,hjust=0.5, alpha=0.5,colour="black")+
    #geom_nodelab(aes(image=tetrapod_id,subset = node == 21), geom="phylopic", align=0,size=0.15,hjust=0.5, alpha=0.5, colour="#47B5FF")+
    geom_nodelab(aes(label="Teleosts",subset = node == 17), hjust=0.5,vjust=2.5, size=3.5, angle=0, fontface=2)+
    geom_nodelab(aes(label="Tetrapods",subset = node == 21), hjust=0.5,vjust=2.5, size=3.5,angle=0,fontface=2)
 
dev.off()

  
  
    


  