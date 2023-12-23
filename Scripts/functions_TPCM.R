############## These functions are used to process and analyze our data ################# 

########### Function written to parse gene trees from Genomicus v 95.01 and to label the nodes with "speciation" and "duplication" events ###########

Parse.trees.Genomicus.95<-function(line)
{
  treetxt <- textConnection(line)
  tree <- treeio::read.nhx(treetxt)
  close(treetxt)
  
  tree@data$label <- c(tree@phylo$tip.label, tree@phylo$node.label)
  serial_no <- nrow(tree@data)
  ntips <- length(tree@phylo$tip.label) 
  internal_nodes <- (ntips+1):(ntips+tree@phylo$Nnode) 
  
  # Adding a column of events of interest
  tree@data$events <- NA
  tree@data$events[ tree@data$D == "Y" ] <- "duplication"
  tree@data$events[ tree@data$D == 0 & tree@data$node > ntips ] <- "speciation"
  tree@data$events <- factor(tree@data$events, levels=c( "speciation", "duplication"))
  return(tree)
}

## Dropping tips without the 15 species IDs (for which tau can be calculated) from a tree. 
## column @phylo$tip.label containing Ensembl.Gene.ID for the species
## From PhyloFish we have only 5 species common between Genomicus & PhyloFishs
## The 5 fish species and 10 outgroup tetrapoda we kept are - (1) zebrafish (^ENSDARG), (2) northern pike (^ENSELUG), 
## (3) cavefish (^ENSAMXG), (4) medaka (^ENSORLG), (5) tilapia (^ENSONIG),(6) human (^ENSG), (7) mouse (^ENSMUSG),
## (8) rat (^ENSRNOG), (9) cow (^ENSBTAG), (10) macaque (^ENSMMUG), (11) chicken (^ENSGALG), (12) dog (^ENSCAFG),
## (13) rabbit (^ENSOCUG), (14) ferret (^ENSMPUG), and (15) opossum (^ENSMODG)

trees_with_specific_IDs <- function(tree) 
{
  ##Collecting tip labels (Ensembl.ID) of trees
  treep <- tree@phylo$tip.label
  
  ##Tips with species data to be retained
  spec <- c("ENSDARG","ENSELUG","ENSAMXG","ENSORLG00000","ENSG000","ENSMUSG","ENSONIG","ENSRNOG","ENSBTAG","ENSMMUG","ENSGALG","ENSCAFG","ENSOCUG","ENSMPUG","ENSMODG")
  search.string <- paste0(spec, collapse = "|")
  Species_data_tip <- treep[grepl(search.string,treep)] 
  
  ## We considered trees with atleast 5 tips
  if(length(Species_data_tip) >= 15)
  {
    tips_to_remove <- treep[!(treep %in% Species_data_tip)]
    pruned.tree<-treeio::drop.tip(tree, tips_to_remove)
    return(pruned.tree)
  }
  else{return(NA)}
}

##This function is written to collect the root node number for a tree
root.node<-function(tree)
{
  ##Collecting gene tree and their data and ntips (number of tips)
  gene.tree<-tree@phylo
  ntips<-as.numeric(length(gene.tree$tip.label))
  
  ## Identifyting root node 
  root.node<-ntips + 1
  return(root.node)
}  

##This function is written to collect the internal nodes for a tree
internal.nodes<-function(tree)
{
  ##Collecting data tibble
  gene.data<-tree@data
  
  ## Identifyting internal nodes with events 
  internal.nodes<-gene.data$node[which(!is.na(gene.data$events))]
  return(internal.nodes)
}  


## This function adds branch length and node depth to our dataframe of tree@data for further use in building time calibration matrix
tree.nodedepth <- function(tree)
{
  gene.tree <- tree@phylo
  
  ##Computing branch length for the trees
  gene.tree <- ape::compute.brlen(gene.tree,1 )
  gene.tree.data <- tree@data
  gene.tree.nodes<-gene.tree.data$node
  
  ##Computing node depth and returned it to the tree into "@data" slot
  nodedepth <- ape::node.depth(gene.tree, method = 1)
  tree@phylo$edge.length<-gene.tree$edge.length
  tree@data$depth_node<-nodedepth
  return(tree)
}

## This function adds node height to our dataframe of tree@data for further use in building time calibration matrix
tree.height <- function(tree)
{
  gene_tree <- tree@phylo
  gene_data <- tree@data
  
  ## Identifying nodes
  gene_tree_nodes<-gene_data$node
  height_node<-NULL
  
  ##Computing node height and returned it to the tree into "@data" slot
  height_node<-ape::node.height(gene_tree)
  tree@data$nodeheight<-height_node
  return(tree)
}


## Dropping tips without tau data of species and considering trees with atleast 4 tips
drop.leaf <- function(tree) 
{
  ##Reading the tips labels and data table
  tree.tip <- tree@phylo$tip.label
  gene.treedata <- tree@data
  
  ##Identify the tips without expression data
  tips.to.remove <- gene.treedata$label[which((is.na(gene.treedata$Tau)) & (is.na(gene.treedata$events)))]
  
  ##Identify the tips with expression data
  tips.with.Tau <- gene.treedata$label[which((!is.na(gene.treedata$Tau)) & (is.na(gene.treedata$events)))]
  
  ##Returning the pruned tree if atleast 4 tips have Tau data
  if(length(tips.with.Tau) >= 4)
  {
    pruned.tree<-treeio::drop.tip(tree, tips.to.remove)
    return(pruned.tree)
  }
  else{return(NA)}
}

##Checking trees with all duplication events 
all.duplication<- function(tree)
{
  ## Collecting total internal node data, no. of tips,speciation and duplication nodes
  gene.tree <- tree@phylo
  gene.treedata <- tree@data
  ntips <- length(gene.tree$tip.label) 
  Internal.nodedata <- nrow(gene.treedata)-ntips
  Duplication.node <- length(gene.treedata$events[which(gene.treedata$events == "duplication")])
  
  ## Now checking for trees with all duplication events
  ## If does not match the criteria returns NA
  if(Internal.nodedata - Duplication.node == 0) 
  {
    return(tree)
  }
  else{
    return(NA)
  }
}

##Checking trees with all speciation events 
all.speciation<- function(tree)
{
  ## Collecting total internal node data, no. of tips,speciation and duplication nodes
  gene.tree <- tree@phylo
  gene.treedata <- tree@data
  ntips <- length(gene.tree$tip.label) 
  Internal.nodedata <- nrow(gene.treedata)-ntips
  Speciation.node <- length(gene.treedata$events[which(gene.treedata$events == "speciation")])
  
  ## Now checking for trees with all speciation events
  ## If does not match the criteria returns NA
  if(Internal.nodedata - Speciation.node == 0) 
  {
    return(tree)
  }
  else{
    return(NA)
  }
}
##Checking trees with all events are duplications events or speciations; we remove them from our analysis and kept trees with atleast one speciation and one duplication node
filter.trees<- function(tree)
{
  ## Collecting total internal node data, no. of tips,speciation and duplication nodes
  gene.tree <- tree@phylo
  gene.treedata <- tree@data
  ntips <- length(gene.tree$tip.label) 
  Internal.nodedata <- nrow(gene.treedata)-ntips
  Duplication.node <- length(gene.treedata$events[ gene.treedata$D == "Y" & gene.treedata$node > ntips ])
  Speciation.node <- length(gene.treedata$events[ gene.treedata$D == 0 & gene.treedata$node > ntips ])
  
  ## Now checking for trees with all speciation and events
  ## If does not match the criteria returns NA
  ifelse((Duplication.node >= 1 & Speciation.node >= 1),return(tree),return(NA))
}

## This function is written to modify the labels of all duplication nodes 
## This function helps us to calibrate speciation nodes properly with different time points
modify.label <- function(tree)
{
  gene.tree <- tree@phylo
  gene.treedata <- tree@data
  
  ## Collecting internal clade labels of duplication nodes to modify them
  Internal.cladelabel <- gene.treedata$label[which(!is.na(gene.treedata$events))]
  dup.node<-gene.treedata$node[which(gene.treedata$events=="duplication")]
  dupnode.datalabel<-gene.treedata$label[dup.node]
  tree@data$label[dup.node] <- paste(dupnode.datalabel,"_d",sep = "")
  return(tree)
}

## This function is written to modify the labels of all "Cluepeocephala" duplication nodes 
modify.label2 <- function(tree)
{
  ## Collecting total internal node data, no. of tips, speciation and duplication nodes
  gene.tree <- tree@phylo
  gene.treedata <- tree@data
  Internal.cladelabel <- gene.treedata$label
  Internal.cladeevent <- gene.treedata$events
  Clupeocephala.dupnodes <- gene.treedata$node[which(Internal.cladelabel=="Clupeocephala" & Internal.cladeevent == "duplication")]
  tree@data$label[Clupeocephala.dupnodes] <- "Clupeocephala_d"
  return(tree)
}

## This function is written to remodify the changed labels of all duplication nodes 
remodify.label <- function(tree)
{
  gene.tree <- tree@phylo
  gene.treedata <- tree@data
  
  
  ## Collecting internal clade labels of duplication nodes to modify them
  Internal.cladelabel <- gene.treedata$label[which(!is.na(gene.treedata$events))]
  Internal.node<-gene.treedata$node[which(!is.na(gene.treedata$events))]
  New.internal.cladelabel<-sapply(Internal.cladelabel,
                                   function(x) 
                                   {x <- unlist(strsplit(toString(x), split='_d', fixed=TRUE))[1]})
  tree@data$label[Internal.node] <-as.character(New.internal.cladelabel)
  return(tree)
}

## This function is written to remodify the previously changed labels of all duplication nodes 
remodify.label2 <- function(tree)
{
  gene.tree <- tree@phylo
  gene.treedata <- tree@data
  
  ##Collecting tiplabels
  tiplabel<-gene.tree$tip.label
  node<-gene.treedata$node
  
  
  ## Collecting internal clade labels of duplication nodes to modify them
  Internal.cladelabel <- gene.treedata$label[which(!is.na(gene.treedata$events))]
  Internal_node<-gene.treedata$node[which(!is.na(gene.treedata$events))]
  New.internal.cladelabel<-as.character(lapply(Internal.cladelabel,
                                                function(x) 
                                                {x <- unlist(strsplit(toString(x), split='_', fixed=TRUE))[1]}))
  
  New.cladelabel<-c(tiplabel,New.internal.cladelabel)
  tree@data$label<- New.cladelabel
  return(tree)
}

## This function adds node height to our dataframe of tree@data for further use in building time calibration matrix
tree.height <- function(tree)
{
  gene.tree <- tree@phylo
  tree@data$nodeheight<-NULL
  
  if (class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  
  ## Identifying nodes
  genetree.nodes<-tree@data$node
  height.node<-NULL
  
  ##Computing node height and returned it to the tree into "@data" slot
  height.node<-ape::node.height(gene.tree)
  tree@data$nodeheight<-height.node
  return(tree)
}

## This function calibrates our gene trees based on the calibration time points
## Maintaining the trees original topology
tree.calibrate <- function(tree, timeframe, model="correlated")
{ 
  
  gene.tree <- tree@phylo
  gene.treedata <- tree@data
  
  ## Identifyting nodes and their depth (those are not representing tips)
  genetree.nodes <- gene.treedata$node[which(!is.na(gene.treedata$events))] 
  node.depth <- gene.treedata$depth_node[genetree.nodes]
  
  ## Merging tree@data based on "label"
  tree.data<-merge(gene.treedata,timeframe,by=c("label"))
  if( nrow(tree.data) > 0 )
  {
    tree.data<- tree.data[order(tree.data$depth_node, decreasing = T),]
    
    ## Buiding calibration matrix
    calibration_matrix <- data.frame(
      node = tree.data$node,
      age.min = tree.data$Mya,
      age.max = tree.data$Mya,
      soft.bounds = NA) 
    ## Time calibrating trees
    calibrate.trees <- try(ape::chronos(gene.tree, calibration = calibration_matrix,  model = model))
    
    ##Trees those are not time calibrated can not be used further
    ##To avoid error due to non calibrated tree we did the following
    if( "phylo" %in% class(calibrate.trees))
    {
      class(calibrate.trees) <- "phylo"
      tree@phylo <- calibrate.trees
    }
    
    else{tree = NA}
    return(tree)
  }
}

## This function adds age of the nodes to the "@data" slot
node.age <-function (tree)
{
  ## Tree and "@data" slots
  gene.tree<-tree@phylo
  gene.data<-tree@data
  
  Internal.label<-gene.data$label[which(!is.na(gene.data$events))]
  Internal.node<-data.frame(node=gene.data$node[which(!is.na(gene.data$events))])
  
  ## Computing node age and adding the values to the "@data"slot
  age.count<-ape::branching.times(gene.tree)
  names(age.count)<-Internal.node
  tree@data$node.age <- c(rep(0, length(gene.tree$tip.label)), as.numeric(age.count))
  return(tree)
} 

### This function helps to identify the list of trees where root event is a speciation event and there is atleast one speciation node before and after a duplication node 
## Therefore the post duplication speciation node should have a nodeheight higher than the nodeheight of the duplication node  and
## The post duplication speciation node should have a node age lower than the age of the duplication node 
special.speciation.criteria.new<-function(tree)
{
  count<<-count+1
  ##Collecting gene data frame
  gene.tree<-tree@phylo
  gene.data<-tree@data
  flag<-NULL
  
  ## Identifyting root node 
  rootnode<-root.node(tree)
  event.root.node<-as.character(gene.data$events[which(gene.data$node==rootnode)]) 
  ntips<-as.numeric(length(gene.tree$tip.label))
  Internal.node<-internal.nodes(tree)
  
  ## Now we need to check if the root node is speciation
  flag<-ifelse(event.root.node %in%"speciation",0, 1)
  
  ## We also need the list of node paths from root to the leaf
  List<-phytools::node.paths(gene.tree,rootnode)
  
  ## If the root event start with speciation, we need to check the length of the List is atleast 3 + 
  ## we also need to check if there are atleast one speciation before and after the duplication
  if(flag==0)
  {
    for(i in 1:length(List))
    {
      path<-NULL
      pathn<-NULL
      path.event<-NULL
      path.age<-NULL
      path<-append(rootnode,as.vector(List[[i]])) ## path should include root node 
      pathn<-path[-length(path)] ## path should exclude tip
      path.event<-as.character(gene.data$events[pathn])
      names(path.event)<-round(gene.data$node.age[pathn],1)
      if((length(path.event) >=3 ) & ("duplication" %in% path.event) & ("speciation" %in% path.event))
      {
        old.speciation<-as.numeric(names(path.event)[which(path.event=="speciation")][[1]])
        old.duplication<-as.numeric(names(path.event)[which(path.event=="duplication")][[1]])
        all.speciation<-as.numeric(names(path.event)[which(path.event=="speciation")])
        young.speciation<-all.speciation[which(all.speciation < old.duplication)]
        
        if((old.speciation > old.duplication) & (length(young.speciation)>0))
        {
          #print(count)
          return(tree)
        }
        else{return(NA)}  
      }
    }
  }
  if(flag==1){return(NA)}
}

## This function helps to make a tibble to find the phylogenetc independent contrast (PIC) of a pre-duplication speciation node and a post-duplication speciation node 
## It also helps to assess the effect of duplication period,predup & postdup speciation age, number of internal duplications between preduplication and postduplication nodes
pic.spe.dup.pair.old<- function(tree,trait)      ## trait can be Tau or Exp
{
  count<<-count+1
  print(count)
  
  ##Collecting gene tree and their data and ntips (number of tips)
  gene.tree<-tree@phylo
  gene.data<-tree@data
  ntips<-as.numeric(length(gene.tree$tip.label))
  
  ## Identifyting internal nodes and root node 
  Internal.node<-gene.data$node[!is.na(gene.data$events)]
  root.node<-ntips + 1
  
  ## Declaring dataframe
  paired.df<-data.frame(olddup.age=NA,
                        predup.speciation=NA,
                        predup.speciation.age=NA,
                        predup.speciation.node.depth=NA,
                        duplication.node=NA,
                        postdup.speciation=NA,
                        postdup.speciation.age=NA,
                        postdup.speciation.node.depth=NA,
                        no.of.duplication=NA)
  
  
  ## Next we plan to get a list of node path to get all possible predup and postdup speciation nodes 
  ## Then we will individually check each path of the list to find if there is a duplication events
  list<-phytools::node.paths(gene.tree,root.node)
  for(i in 1:length(list))
  {
    ## reading the elements of the list as a dataframe
    element<-NULL 
    dupnode.old<-NULL
    dup.index<-NULL
    spenode.older<-NULL
    spenode.older.age<-NULL
    spenode.older.depth<-NULL
    spenode.older.index<-NULL
    spenode.younger<-NULL
    spenode.younger.age<-NULL
    spenode.younger.depth<-NULL
    spenode.younger.index<-NULL
    dupnode.old.age<-NULL
    element<-data.frame(node=root.node)
    element<-rbind(element,data.frame(node=as.numeric(list[[i]])))
    element$events<-gene.data$events[element$node]
    element$age<-gene.data$node.age[element$node]
    element$depth<-gene.data$depth_node[element$node]
    element$index<-as.numeric(row.names(element))
    
    ## Identifying the oldest duplication node (among all duplication nodes) in a path to identify preduplication speciation node 
    dupnode.all<-as.vector(element$node[which(element$events %in% "duplication")])
    total.dupnodes<-as.numeric(length(dupnode.all))
    speciation.all<-as.vector(element$node[which(element$events %in% "speciation")])
    names(dupnode.all)<-round(element$age[which(element$node %in% dupnode.all)],1)
    names(speciation.all)<-round(element$age[which(element$node %in% speciation.all)],1)
    spenode.older1<-as.numeric(speciation.all[1])
    spenode.older1.age<-round(element$age[which(element$node %in% spenode.older1)],1)
    
    if(length(spenode.older1) > 0)
    {
      old.duplication<-NULL
      old.duplication<-as.numeric(dupnode.all)[which(as.numeric(names(dupnode.all)) < spenode.older1.age)][1]
      dup.index<-element$index[which(element$node==old.duplication)]
      old.duplication.age<-element$age[dup.index]
    
      ## Collecting pre-duplication speciation information 
      dupnode.old<-as.numeric(old.duplication)
      spenode.older<-element$node[dup.index-1]
      spenode.older.age<-element$age[dup.index-1]
      spenode.older.depth<-element$depth[dup.index-1]
        
      ## Considering the oldest duplication node 
      if(length(dup.index)> 0)
      {
        for(j in (dup.index+1):(nrow(element)))
        {
          new.event<-element$events[j] 
          if(new.event %in% "speciation")
          {
            spenode.younger<-element$node[j]
            spenode.younger.age<-element$age[j]
            spenode.younger.depth<-element$depth[j]
            spenode.younger.index<-element$index[j]
            break
          }
        }
      }
    }
    ## Collection and pairing up preduplication speciation node & postduplication speciation node
    if((nrow(element)>=4) & (length(spenode.older) > 0) & (length(spenode.younger > 0)))
    {
      if((!is.na(dupnode.old)) & (spenode.younger > ntips))
      {
        paired.df<-rbind(paired.df,data.frame(olddup.age=old.duplication.age,
                                              predup.speciation=spenode.older,
                                              predup.speciation.age=spenode.older.age,
                                              predup.speciation.node.depth=spenode.older.depth,
                                              duplication.node=dupnode.old,
                                              postdup.speciation=spenode.younger,
                                              #postdup_speciation=dupnode_old,
                                              postdup.speciation.age=spenode.younger.age,
                                              postdup.speciation.node.depth=spenode.younger.depth,
                                              no.of.duplication= as.numeric(spenode.younger.index - dup.index)))
      }
    }
  }
  
  ## Removing the row with NA value
  paired.df<-paired.df[-1,]
  
  ## Removing duplicated entries 
  paired.df<- paired.df[!duplicated(paired.df),]
  
  ##Removing redundant rows
  if((trait=="pic_Exp") & (nrow(paired.df) >= 1))
  {
    
    ##Returning the tibble
    return(tibble(
      tree.num = count,
      gene = digest(tree),
      olddup.age = paired.df$olddup.age,
      predup.spenode = paired.df$predup.speciation,
      predup.spe.age = paired.df$predup.speciation.age,
      predup.spe.label=gene.data$label[paired.df$predup.speciation],
      predup.spe.pic = abs(gene.data$pic_Exp[paired.df$predup.speciation]),
      predup.spe.node.depth = paired.df$predup.speciation.node.depth,
      predup.spe.var = gene.data$variance[paired.df$predup.speciation],
      duplication.pic = abs(gene.data$pic_Exp[paired.df$duplication.node]),
      dupnode=paired.df$duplication.node,
      duplication.label=gene.data$label[paired.df$duplication.node],
      postdup.spenode = paired.df$postdup.speciation,
      postdup.spe.age = paired.df$postdup.speciation.age,
      postdup.spe.label=gene.data$label[paired.df$postdup.speciation],
      postdup.spe.pic = abs(gene.data$pic_Exp[paired.df$postdup.speciation]),
      postdup.spe.node.depth = paired.df$postdup.speciation.node.depth,
      postdup.spe.var = gene.data$variance[paired.df$postdup.speciation],
      duplication.num = paired.df$no.of.duplication))
  }
  if((trait=="pic_Exp_rand") & (nrow(paired.df) >= 1))
  {
    ##Returning the tibble
    return(tibble(
      tree.num = count,
      gene = digest(tree),
      olddup.age = paired.df$olddup.age,
      predup.spenode = paired.df$predup.speciation,
      predup.spe.age = paired.df$predup.speciation.age,
      predup.spe.label=gene.data$label[paired.df$predup.speciation],
      predup.spe.pic = abs(gene.data$pic[paired.df$predup.speciation]),
      predup.spe.node.depth = paired.df$predup.speciation.node.depth,
      predup.spe.var = gene.data$variance[paired.df$predup.speciation],
      duplication.pic = abs(gene.data$pic[paired.df$duplication.node]),
      dupnode=paired.df$duplication.node,
      duplication.label=gene.data$label[paired.df$duplication.node],
      postdup.spenode = paired.df$postdup.speciation,
      postdup.spe.age = paired.df$postdup.speciation.age,
      postdup.spe.label=gene.data$label[paired.df$postdup.speciation],
      postdup.spe.pic = abs(gene.data$pic[paired.df$postdup.speciation]),
      postdup.spe.node.depth = paired.df$postdup.speciation.node.depth,
      postdup.spe.var = gene.data$variance[paired.df$postdup.speciation],
      duplication.num = paired.df$no.of.duplication))
  }
  if((trait=="pic_Tau") & (nrow(paired.df) >= 1))
  {
    ##Returning the tibble
    return(tibble(
      tree.num = count,
      gene = digest(tree),
      olddup.age = paired.df$olddup.age,
      predup.spenode = paired.df$predup.speciation,
      predup.spe.age = paired.df$predup.speciation.age,
      predup.spe.label=gene.data$label[paired.df$predup.speciation],
      predup.spe.pic = abs(gene.data$pic_Tau[paired.df$predup.speciation]),
      predup.spe.node.depth = paired.df$predup.speciation.node.depth,
      predup.spe.var = gene.data$variance[paired.df$predup.speciation],
      duplication.pic = abs(gene.data$pic_Tau[paired.df$duplication.node]),
      dupnode=paired.df$duplication.node,
      duplication.label=gene.data$label[paired.df$duplication.node],
      postdup.spenode = paired.df$postdup.speciation,
      postdup.spe.age = paired.df$postdup.speciation.age,
      postdup.spe.label=gene.data$label[paired.df$postdup.speciation],
      postdup.spe.pic = abs(gene.data$pic_Tau[paired.df$postdup.speciation]),
      postdup.spe.node.depth = paired.df$postdup.speciation.node.depth,
      postdup.spe.var = gene.data$variance[paired.df$postdup.speciation],
      duplication.num = paired.df$no.of.duplication))
  }
}

### This function helps to identify the list of 3RWGD trees (where root event is a speciation event) + atleast one speciation node before and after a FishWGD node 
## Therefore the post duplication speciation node should have a nodeheight higher than the nodeheight of the previous FishWGD node  and
## The post duplication speciation node should have a node age lower than the age of the previous FishWGD node
## We consider oldest FishWGD among many in this case
special.speciation.criteria.3RWGD.old <- function(tree)
{
  count<<-count+1
  print(count)
  ##Collecting gene data frame
  gene.tree<-tree@phylo
  gene.data<-tree@data
  flag<-NULL
  
  ## Identifyting root node 
  rootnode<-root.node(tree)
  event.root.node<-as.character(gene.data$events[which(gene.data$node==rootnode)]) 
  ntips<-as.numeric(length(gene.tree$tip.label))
  Internal.node<-internal.nodes(tree)
  
  ## Now we need to check if the root node is speciation
  flag<-ifelse(event.root.node %in%"speciation",0, 1)
  
  ## Root event may be a speciation/duplication event
  #flag<-0
  
  ## We also need the list of node paths from root to the leaf
  List<-phytools::node.paths(gene.tree,rootnode)
  
  ## Then we need to check the length of the events in a path is atleast 3 + 
  ## if there are atleast one speciation before and after the whole genome duplication
  if(flag==0)
  {
    for(i in 1:length(List))
    {
      path<-NULL
      pathn<-NULL
      path.event<-NULL
      path.age<-NULL
      path<-append(rootnode,as.vector(List[[i]])) ## path should include root node 
      pathn<-path[-length(path)] ## path should exclude tip
      path.event<-as.character(gene.data$events[pathn])
      names(path.event)<-round(gene.data$node.age[pathn],1)
      if((length(path.event) >=3 ) & ("FishWGD" %in% path.event) & ("speciation" %in% path.event))
      {
        old.speciation<-as.numeric(names(path.event)[which(path.event=="speciation")][[1]])
        old.duplication<-as.numeric(names(path.event)[which(path.event=="FishWGD")][[1]])
        all.speciation<-as.numeric(names(path.event)[which(path.event=="speciation")])
        young.speciation<-all.speciation[which(all.speciation < old.duplication)]
        
        if((old.speciation > old.duplication) & (length(young.speciation)>0))
        {
          #print(count)
          return(tree)
        }
        else{return(NA)}  
      }
    }
  }
  if(flag==1){return(NA)}
}

## This function helps to make a tibble to find the phylogenetc independent contrast (PIC) of a pre-duplication speciation node and a post-duplication speciation node for 3RWGD trees
## This considers old WGD event
## It also helps to assess the effect of duplication period,predup & postdup speciation age, number of internal duplications between preduplication and postduplication nodes
pic.speciation.pair.3RWGD<- function(tree,trait)  ## trait can be Tau or Exp
{
  count<<-count+1
  print(count)
  
  ##Collecting gene tree and their data and ntips (number of tips)
  gene.tree<-tree@phylo
  gene.data<-tree@data
  ntips<-as.numeric(length(gene.tree$tip.label))
  
  ## Identifyting internal nodes and root node 
  Internal.node<-gene.data$node[!is.na(gene.data$events)]
  root.node<-ntips + 1
  
  ## Declaring dataframe
  paired.df<-data.frame(olddup.age=NA,
                        predup.speciation=NA,
                        predup.speciation.age=NA,
                        predup.speciation.node.depth=NA,
                        duplication.node=NA,
                        postdup.speciation=NA,
                        postdup.speciation.age=NA,
                        postdup.speciation.node.depth=NA,
                        no.of.duplication=NA)
  
  
  ## Next we plan to get a list of node path to get all possible predup and postdup speciation nodes 
  ## Then we will individually check each path of the list to find if there is a duplication events
  list<-phytools::node.paths(gene.tree,root.node)
  for(i in 1:length(list))
  {
    ## reading the elements of the list as a dataframe
    element<-NULL 
    dupnode.old<-NULL
    dupnode.young<-NULL
    dup.index<-NULL
    dup.index.young<-NULL
    spenode.before.duplication<-NULL
    spenode.older<-NULL
    spenode.older.age<-NULL
    spenode.older.depth<-NULL
    spenode.younger<-NULL
    spenode.younger.age<-NULL
    spenode.younger.depth<-NULL
    spenode.younger.index<-NULL
    dupnode.old.age<-NULL
    element<-data.frame(node=root.node)
    element<-rbind(element,data.frame(node=as.numeric(list[[i]])))
    element$events<-gene.data$events[element$node]
    element$age<-gene.data$node.age[element$node]
    element$depth<-gene.data$depth_node[element$node]
    element$index<-as.numeric(row.names(element))
    
    ## Identifying the oldest duplication node (among all duplication nodes) in a path to identify preduplication speciation node 
    dupnode.all<-as.vector(element$node[which(element$events %in% "FishWGD" | element$events %in% "duplication")])
    total.dupnodes<-length(dupnode.all)
    
    ## Collecting the oldest duplication node 
    dupnode.old<-as.numeric(element$node[which(element$events %in% "FishWGD")][1])
    dup.index<-element$index[which(element$node==dupnode.old)]
    dupnode.old.age<-element$age[dup.index]
    #spenode.older<-element$node[dup.index-1]
    
    spenode.before.duplication<-as.numeric(element$node[which(element$events %in% "speciation" & element$index < dup.index)])
    spenode.older<-spenode.before.duplication[length(spenode.before.duplication)]
    spenode.older.age<-element$age[which(element$node==spenode.older)]
    spenode.older.depth<-element$depth[which(element$node==spenode.older)]
    
    ## Identifying the next speciation node and number of duplication in between
    if(length(dup.index)> 0)
    {
      for(j in (dup.index+1):(nrow(element)))
      {
        new.event<-element$events[j] 
        if(new.event %in% "speciation")
        {
          spenode.younger<-element$node[j]
          spenode.younger.age<-element$age[j]
          spenode.younger.depth<-element$depth[j]
          spenode.younger.index<-element$index[j]
          break
        }
      }
    }
    
    ## Collection and pairing up preduplication speciation node & postduplication speciation node
    if((nrow(element)>=4) & (length(spenode.older) > 0) & (length(spenode.younger > 0)))
    {
      if((!is.na(dupnode.old))  & (spenode.younger > ntips))
      {
        paired.df<-rbind(paired.df,data.frame(olddup.age=dupnode.old.age,
                                              predup.speciation=spenode.older,
                                              predup.speciation.age=spenode.older.age,
                                              predup.speciation.node.depth=spenode.older.depth,
                                              duplication.node=dupnode.old,
                                              postdup.speciation=spenode.younger,
                                              postdup.speciation.age=spenode.younger.age,
                                              postdup.speciation.node.depth=spenode.younger.depth,
                                              no.of.duplication= as.numeric(spenode.younger.index - dup.index)))
      }
    }
  }
  
  ## Removing the row with NA value
  paired.df<-paired.df[-1,]
  
  ##Removing redundant rows
  if((trait=="pic_Exp") & (nrow(paired.df) >= 1))
  {
    paired.df<- paired.df[!duplicated(paired.df),]
    
    ##Returning the tibble
    return(tibble(
      tree.num = count,
      gene = digest(tree),
      olddup.age = paired.df$olddup.age,
      predup.spenode = paired.df$predup.speciation,
      predup.spe.age = paired.df$predup.speciation.age,
      predup.spe.label=gene.data$label[paired.df$predup.speciation],
      predup.spe.pic = abs(gene.data$pic_Exp[paired.df$predup.speciation]),
      predup.spe.node.depth = paired.df$predup.speciation.node.depth,
      predup.spe.var = gene.data$variance[paired.df$predup.speciation],
      duplication.pic = abs(gene.data$pic_Exp[paired.df$duplication.node]),
      dupnode=paired.df$duplication.node,
      duplication.label=gene.data$label[paired.df$duplication.node],
      postdup.spenode = paired.df$postdup.speciation,
      postdup.spe.age = paired.df$postdup.speciation.age,
      postdup.spe.label=gene.data$label[paired.df$postdup.speciation],
      postdup.spe.pic = abs(gene.data$pic_Exp[paired.df$postdup.speciation]),
      postdup.spe.node.depth = paired.df$postdup.speciation.node.depth,
      postdup.spe.var = gene.data$variance[paired.df$postdup.speciation],
      duplication.num = paired.df$no.of.duplication))
  }
  if((trait=="pic_Exp_rand") & (nrow(paired.df) >= 1))
  {
    paired.df<- paired.df[!duplicated(paired.df),]
    
    ##Returning the tibble
    return(tibble(
      tree.num = count,
      gene = digest(tree),
      olddup.age = paired.df$olddup.age,
      predup.spenode = paired.df$predup.speciation,
      predup.spe.age = paired.df$predup.speciation.age,
      predup.spe.label=gene.data$label[paired.df$predup.speciation],
      predup.spe.pic = abs(gene.data$pic[paired.df$predup.speciation]),
      predup.spe.node.depth = paired.df$predup.speciation.node.depth,
      predup.spe.var = gene.data$variance[paired.df$predup.speciation],
      duplication.pic = abs(gene.data$pic[paired.df$duplication.node]),
      dupnode=paired.df$duplication.node,
      duplication.label=gene.data$label[paired.df$duplication.node],
      postdup.spenode = paired.df$postdup.speciation,
      postdup.spe.age = paired.df$postdup.speciation.age,
      postdup.spe.label=gene.data$label[paired.df$postdup.speciation],
      postdup.spe.pic = abs(gene.data$pic[paired.df$postdup.speciation]),
      postdup.spe.node.depth = paired.df$postdup.speciation.node.depth,
      postdup.spe.var = gene.data$variance[paired.df$postdup.speciation],
      duplication.num = paired.df$no.of.duplication))
  }
  if((trait=="pic_Tau") & (nrow(paired.df) >= 1))
  {
    paired.df<- paired.df[!duplicated(paired.df),]
    
    ##Returning the tibble
    return(tibble(
      tree.num = count,
      gene = digest(tree),
      olddup.age = paired.df$olddup.age,
      predup.spenode = paired.df$predup.speciation,
      predup.spe.age = paired.df$predup.speciation.age,
      predup.spe.label=gene.data$label[paired.df$predup.speciation],
      predup.spe.pic = abs(gene.data$pic_Tau[paired.df$predup.speciation]),
      predup.spe.node.depth = paired.df$predup.speciation.node.depth,
      predup.spe.var = gene.data$variance[paired.df$predup.speciation],
      duplication.pic = abs(gene.data$pic_Tau[paired.df$duplication.node]),
      dupnode=paired.df$duplication.node,
      duplication.label=gene.data$label[paired.df$duplication.node],
      postdup.spenode = paired.df$postdup.speciation,
      postdup.spe.age = paired.df$postdup.speciation.age,
      postdup.spe.label=gene.data$label[paired.df$postdup.speciation],
      postdup.spe.pic = abs(gene.data$pic_Tau[paired.df$postdup.speciation]),
      postdup.spe.node.depth = paired.df$postdup.speciation.node.depth,
      postdup.spe.var = gene.data$variance[paired.df$postdup.speciation],
      duplication.num = paired.df$no.of.duplication))
  }
}

### This function helps to identify the list of 3RWGD trees (where root event is a speciation event) + atleast one speciation node before and after a FishWGD node 
## Therefore the post duplication speciation node should have a nodeheight higher than the nodeheight of the previous FishWGD node  and
## The post duplication speciation node should have a node age lower than the age of the previous FishWGD node
## We consider youngest FishWGD among many in this case
special.speciation.criteria.3RWGD.young <- function(tree)
{
  count<<-count+1
  print(count)
  
  ##Collecting gene data frame
  gene.tree<-tree@phylo
  gene.data<-tree@data
  flag<-NULL
  
  ## Identifyting root node 
  rootnode<-root.node(tree)
  event.root.node<-as.character(gene.data$events[which(gene.data$node==rootnode)]) 
  ntips<-as.numeric(length(gene.tree$tip.label))
  Internal.node<-internal.nodes(tree)
  
  ## Now we need to check if the root node is speciation
  flag<-ifelse(event.root.node %in%"speciation",0, 1)
  
  ## We also need the list of node paths from root to the leaf
  List<-phytools::node.paths(gene.tree,rootnode)
  
  ## Then we need to check the length of the events in a path is atleast 3 + 
  ## if there are atleast one speciation before and after the whole genome duplication
  if(flag==0)
  {
    for(i in 1:length(List))
    {
      path<-NULL
      pathn<-NULL
      path.event<-NULL
      path.age<-NULL
      path<-append(rootnode,as.vector(List[[i]])) ## path should include root node 
      pathn<-path[-length(path)] ## path should exclude tip
      path.event<-as.character(gene.data$events[pathn])
      names(path.event)<-round(gene.data$node.age[pathn],1)
      if((length(path.event) >=3 ) & ("FishWGD" %in% path.event) & ("speciation" %in% path.event))
      {
        all.speciation<-as.numeric(names(path.event)[which(path.event=="speciation")])
        all.WGD<-as.numeric(names(path.event)[which(path.event=="FishWGD")])
        if(length(all.WGD) > 0)
        {
          young.WGD<-all.WGD[length(all.WGD)]
          spenodes.before.WGD<-all.speciation[which(all.speciation > young.WGD)]
          spenode.older<-spenodes.before.WGD[length(spenodes.before.WGD)]
          young.speciation<-all.speciation[which(all.speciation < young.WGD)][1]
        
        if((spenode.older > young.WGD) & (length(young.speciation)>0))
        {
          #print(count)
          return(tree)
        }}
        else{return(NA)}  
      }
    }
  }
  if(flag==1){return(NA)}
}

## This function helps to make a tibble to find the phylogenetc independent contrast (PIC) of a pre-duplication speciation node and a post-duplication speciation node for 3RWGD trees
## This considers young FishWGD event
## It also helps to assess the effect of duplication period,predup & postdup speciation age, number of internal duplications between preduplication and postduplication nodes
pic.speciation.pair.3RWGD.young<- function(tree,trait)  ## trait can be Tau or Exp
{
  count<<-count+1
  print(count)
  
  ##Collecting gene tree and their data and ntips (number of tips)
  gene.tree<-tree@phylo
  gene.data<-tree@data
  ntips<-as.numeric(length(gene.tree$tip.label))
  
  ## Identifyting internal nodes and root node 
  Internal.node<-gene.data$node[!is.na(gene.data$events)]
  root.node<-ntips + 1
  
  ## Declaring dataframe
  paired.df<-data.frame(youngdup.age=NA,
                        predup.speciation=NA,
                        predup.speciation.age=NA,
                        predup.speciation.node.depth=NA,
                        duplication.node=NA,
                        postdup.speciation=NA,
                        postdup.speciation.age=NA,
                        postdup.speciation.node.depth=NA,
                        no.of.duplication=NA)
  
  
  ## Next we plan to get a list of node path to get all possible predup and postdup speciation nodes 
  ## Then we will individually check each path of the list to find if there is a duplication events
  list<-phytools::node.paths(gene.tree,root.node)
  for(i in 1:length(list))
  {
    ## reading the elements of the list as a dataframe
    element<-NULL
    dupnode.WGD<-NULL
    dupnode.young<-NULL
    dup.index<-NULL
    dup.index.young<-NULL
    spenode.before.duplication<-NULL
    spenode.older<-NULL
    spenode.older.age<-NULL
    spenode.older.depth<-NULL
    spenode.younger<-NULL
    spenode.younger.age<-NULL
    spenode.younger.depth<-NULL
    spenode.younger.index<-NULL
    dupnode.old.age<-NULL
    element<-data.frame(node=root.node)
    element<-rbind(element,data.frame(node=as.numeric(list[[i]])))
    element$events<-gene.data$events[element$node]
    element$age<-gene.data$node.age[element$node]
    element$depth<-gene.data$depth_node[element$node]
    element$index<-as.numeric(row.names(element))
    
    ## Identifying all duplication nodes 
    dupnode.all<-as.vector(element$node[which(element$events %in% "FishWGD" | element$events %in% "duplication")])
    total.dupnodes<-length(dupnode.all)
    
    ## Collecting the youngest WGD node 
    dupnode.WGD<-as.numeric(element$node[which(element$events %in% "FishWGD")])
    dupnode.young<-dupnode.WGD[length(dupnode.WGD)]
    dup.index<-element$index[which(element$node==dupnode.young)]
    dupnode.young.age<-element$age[dup.index] 

    ## Collecting old speciation node
    spenode.before.duplication<-as.numeric(element$node[which(element$events %in% "speciation" & element$index < dup.index)])
    spenode.older<-spenode.before.duplication[length(spenode.before.duplication)]
    spenode.older.age<-element$age[which(element$node==spenode.older)]
    spenode.older.depth<-element$depth[which(element$node==spenode.older)]
    
    ## Identifying the next speciation node and number of duplication in between
    if(length(dup.index)> 0)
    {
      for(j in (dup.index+1):(nrow(element)))
      {
        new.event<-element$events[j] 
        if(new.event %in% "speciation")
        {
          spenode.younger<-element$node[j]
          spenode.younger.age<-element$age[j]
          spenode.younger.depth<-element$depth[j]
          spenode.younger.index<-element$index[j]
          break
        }
      }
    }
    
    ## Collection and pairing up preduplication speciation node & postduplication speciation node
    if((nrow(element)>=4) & (length(spenode.older) > 0) & (length(spenode.younger > 0)))
    {
      #if((!is.na(dupnode.old)) & (!is.na(dupnode.young)) & (spenode.younger > ntips))
      if((!is.na(dupnode.young))  & (spenode.younger > ntips))
      {
        paired.df<-rbind(paired.df,data.frame(youngdup.age=dupnode.young.age,
                                              predup.speciation=spenode.older,
                                              predup.speciation.age=spenode.older.age,
                                              predup.speciation.node.depth=spenode.older.depth,
                                              duplication.node=dupnode.young,
                                              postdup.speciation=spenode.younger,
                                              postdup.speciation.age=spenode.younger.age,
                                              postdup.speciation.node.depth=spenode.younger.depth,
                                              no.of.duplication= as.numeric(spenode.younger.index - dup.index)))
      }
    }
  }
  
  ## Removing the row with NA value
  paired.df<-paired.df[-1,]
  
  ##Removing redundant rows
  if((trait=="pic_Exp") & (nrow(paired.df) >= 1))
  {
    paired.df<- paired.df[!duplicated(paired.df),]
    
    ##Returning the tibble
    return(tibble(
      tree.num = count,
      gene = digest(tree),
      youngdup.age = paired.df$youngdup.age,
      predup.spenode = paired.df$predup.speciation,
      predup.spe.age = paired.df$predup.speciation.age,
      predup.spe.label=gene.data$label[paired.df$predup.speciation],
      predup.spe.pic = abs(gene.data$pic_Exp[paired.df$predup.speciation]),
      predup.spe.node.depth = paired.df$predup.speciation.node.depth,
      predup.spe.var = gene.data$variance[paired.df$predup.speciation],
      duplication.pic = abs(gene.data$pic_Exp[paired.df$duplication.node]),
      dupnode=paired.df$duplication.node,
      duplication.label=gene.data$label[paired.df$duplication.node],
      postdup.spenode = paired.df$postdup.speciation,
      postdup.spe.age = paired.df$postdup.speciation.age,
      postdup.spe.label=gene.data$label[paired.df$postdup.speciation],
      postdup.spe.pic = abs(gene.data$pic_Exp[paired.df$postdup.speciation]),
      postdup.spe.node.depth = paired.df$postdup.speciation.node.depth,
      postdup.spe.var = gene.data$variance[paired.df$postdup.speciation],
      duplication.num = paired.df$no.of.duplication))
  }
  if((trait=="pic_Exp_rand") & (nrow(paired.df) >= 1))
  {
    paired.df<- paired.df[!duplicated(paired.df),]
    
    ##Returning the tibble
    return(tibble(
      tree.num = count,
      gene = digest(tree),
      youngdup.age = paired.df$youngdup.age,
      predup.spenode = paired.df$predup.speciation,
      predup.spe.age = paired.df$predup.speciation.age,
      predup.spe.label=gene.data$label[paired.df$predup.speciation],
      predup.spe.pic = abs(gene.data$pic[paired.df$predup.speciation]),
      predup.spe.node.depth = paired.df$predup.speciation.node.depth,
      predup.spe.var = gene.data$variance[paired.df$predup.speciation],
      duplication.pic = abs(gene.data$pic[paired.df$duplication.node]),
      dupnode=paired.df$duplication.node,
      duplication.label=gene.data$label[paired.df$duplication.node],
      postdup.spenode = paired.df$postdup.speciation,
      postdup.spe.age = paired.df$postdup.speciation.age,
      postdup.spe.label=gene.data$label[paired.df$postdup.speciation],
      postdup.spe.pic = abs(gene.data$pic[paired.df$postdup.speciation]),
      postdup.spe.node.depth = paired.df$postdup.speciation.node.depth,
      postdup.spe.var = gene.data$variance[paired.df$postdup.speciation],
      duplication.num = paired.df$no.of.duplication))
  }
  if((trait=="pic_Tau") & (nrow(paired.df) >= 1))
  {
    paired.df<- paired.df[!duplicated(paired.df),]
    
    ##Returning the tibble
    return(tibble(
      tree.num = count,
      gene = digest(tree),
      youngdup.age = paired.df$youngdup.age,
      predup.spenode = paired.df$predup.speciation,
      predup.spe.age = paired.df$predup.speciation.age,
      predup.spe.label=gene.data$label[paired.df$predup.speciation],
      predup.spe.pic = abs(gene.data$pic_Tau[paired.df$predup.speciation]),
      predup.spe.node.depth = paired.df$predup.speciation.node.depth,
      predup.spe.var = gene.data$variance[paired.df$predup.speciation],
      duplication.pic = abs(gene.data$pic_Tau[paired.df$duplication.node]),
      dupnode=paired.df$duplication.node,
      duplication.label=gene.data$label[paired.df$duplication.node],
      postdup.spenode = paired.df$postdup.speciation,
      postdup.spe.age = paired.df$postdup.speciation.age,
      postdup.spe.label=gene.data$label[paired.df$postdup.speciation],
      postdup.spe.pic = abs(gene.data$pic_Tau[paired.df$postdup.speciation]),
      postdup.spe.node.depth = paired.df$postdup.speciation.node.depth,
      postdup.spe.var = gene.data$variance[paired.df$postdup.speciation],
      duplication.num = paired.df$no.of.duplication))
  }
}


## This function helps to make a tibble to find the phylogenetc independent contrast (PIC) of a pre-duplication speciation node and a post-duplication speciation node 
## Considering the most recent duplication events

pic.spe.dup.pair.young<- function(tree,trait)
{
  count<<-count+1
  print(count)
  
  ##Collecting gene tree and their data and ntips (number of tips)
  gene.tree<-tree@phylo
  gene.data<-tree@data
  ntips<-as.numeric(length(gene.tree$tip.label))
  
  ## Identifyting internal nodes and root node 
  Internal.node<-gene.data$node[!is.na(gene.data$events)]
  root.node<-ntips + 1
  
  ## Declaring dataframe
  paired.df<-data.frame(youngdup.age=NA,
                        predup.speciation=NA,
                        predup.speciation.age=NA,
                        predup.speciation.node.depth=NA,
                        duplication.node=NA,
                        postdup.speciation=NA,
                        postdup.speciation.age=NA,
                        postdup.speciation.node.depth=NA,
                        no.of.duplication=NA)
  
  
  ## Next we plan to get a list of node path to get all possible predup and postdup speciation nodes 
  ## Then we will individually check each path of the list to find if there is a duplication events
  list<-phytools::node.paths(gene.tree,root.node)
  for(i in 1:length(list))
  {
    ## Declaring variables
    element<-NULL 
    dupnode.old<-NULL
    dupnode.young<-NULL
    dup.index.young<-NULL
    spenode.older<-NULL
    spenode.older.age<-NULL
    spenode.older.depth<-NULL
    spenode.younger<-NULL
    spenode.younger.age<-NULL
    spenode.younger.depth<-NULL
    spenode.older.index<-NULL
    dupnode.young.age<-NULL
    
    ## reading the elements of the list as a dataframe
    element<-data.frame(node=root.node)
    element<-rbind(element,data.frame(node=as.numeric(list[[i]])))
    element$events<-gene.data$events[element$node]
    element$age<-gene.data$node.age[element$node]
    element$depth<-gene.data$depth_node[element$node]
    element$index<-as.numeric(row.names(element))
    
    ## Here we plan to identify youngest duplication node (among all duplication nodes) to study post duplication bias
    dupnode.all<-as.vector(element$node[which(element$events %in% "duplication")])
    total.dupnodes<-length(dupnode.all)
    
    ## Collecting the oldest duplication node 
    dupnode.old<-dupnode.all[1]
    
    ## Collecting the youngest duplication node 
    dupnode.young<-dupnode.all[total.dupnodes]
    dup.index.young<- element$index[which(element$node==dupnode.young)]
    dupnode.young.age<-element$age[dup.index.young]
    spenode.younger<-element$node[dup.index.young+1]
    spenode.younger.age<-element$age[dup.index.young+1]
    spenode.younger.depth<-element$depth[dup.index.young+1]
    
    ## Identifying the next speciation node and number of duplication in between
    if((length(dup.index.young)> 0) & (length(spenode.younger) > 0))
    {
      for(j in (dup.index.young-1):1)
      {
        new.event<-element$events[j] 
        if(new.event %in% "speciation")
        {
          spenode.older<-element$node[j]
          spenode.older.age<-element$age[j]
          spenode.older.depth<-element$depth[j]
          spenode.older.index<-element$index[j]
          break
        }
      }
    }
    
    ## Collection and pairing up preduplication speciation node & postduplication speciation node
    if((nrow(element)>=4) & (length(spenode.older) > 0) & (length(spenode.younger > 0)))
    {
      if((!is.na(dupnode.old)) & (!is.na(dupnode.young)) & (spenode.younger > ntips))
      {
        paired.df<-rbind(paired.df,data.frame(youngdup.age=dupnode.young.age,
                                              predup.speciation=spenode.older,
                                              predup.speciation.age=spenode.older.age,
                                              predup.speciation.node.depth=spenode.older.depth,
                                              duplication.node=dupnode.young,
                                              postdup.speciation=spenode.younger,
                                              postdup.speciation.age=spenode.younger.age,
                                              postdup.speciation.node.depth=spenode.younger.depth,
                                              no.of.duplication= as.numeric(dup.index.young - spenode.older.index)))
      }
    }
  }
  
  ## Removing the row with NA value
  paired.df<-paired.df[-1,]
  
  ##Removing redundant rows
  if((trait=="pic_Exp") & (nrow(paired.df) >= 1))
  {
    paired.df<- paired.df[!duplicated(paired.df),]
    
    ##Returning the tibble
    return(tibble(
      tree.num=count,
      gene = digest(tree),
      recentdup.age = paired.df$youngdup.age,
      predup.spenode = paired.df$predup.speciation,
      predup.spe.age = paired.df$predup.speciation.age,
      predup.spe.label=gene.data$label[paired.df$predup.speciation],
      predup.spe.pic = abs(gene.data$pic_Exp[paired.df$predup.speciation]),
      predup.spe.node.depth = paired.df$predup.speciation.node.depth,
      predup.spe.var = gene.data$variance[paired.df$predup.speciation],
      duplication.pic = abs(gene.data$pic_Exp[paired.df$duplication.node]),
      dupnode=paired.df$duplication.node,
      duplication.label=gene.data$label[paired.df$duplication.node],
      postdup.spenode = paired.df$postdup.speciation,
      postdup.spe.age = paired.df$postdup.speciation.age,
      postdup.spe.label=gene.data$label[paired.df$postdup.speciation],
      postdup.spe.pic = abs(gene.data$pic_Exp[paired.df$postdup.speciation]),
      postdup.spe.node.depth = paired.df$postdup.speciation.node.depth,
      postdup.spe.var = gene.data$variance[paired.df$postdup.speciation],
      duplication.num = paired.df$no.of.duplication))
  }
  
  ##Removing redundant rows
  if((trait=="pic_Tau") & (nrow(paired.df) >= 1))
  {
    paired.df<- paired.df[!duplicated(paired.df),]
    
    ##Returning the tibble
    return(tibble(
      tree.num=count,
      gene = digest(tree),
      recentdup.age = paired.df$youngdup.age,
      predup.spenode = paired.df$predup.speciation,
      predup.spe.age = paired.df$predup.speciation.age,
      predup.spe.label=gene.data$label[paired.df$predup.speciation],
      predup.spe.pic = abs(gene.data$pic_Tau[paired.df$predup.speciation]),
      predup.spe.node.depth = paired.df$predup.speciation.node.depth,
      predup.spe.var = gene.data$variance[paired.df$predup.speciation],
      duplication.pic = abs(gene.data$pic_Tau[paired.df$duplication.node]),
      dupnode=paired.df$duplication.node,
      duplication.label=gene.data$label[paired.df$duplication.node],
      postdup.spenode = paired.df$postdup.speciation,
      postdup.spe.age = paired.df$postdup.speciation.age,
      postdup.spe.label=gene.data$label[paired.df$postdup.speciation],
      postdup.spe.pic = abs(gene.data$pic_Tau[paired.df$postdup.speciation]),
      postdup.spe.node.depth = paired.df$postdup.speciation.node.depth,
      postdup.spe.var = gene.data$variance[paired.df$postdup.speciation],
      duplication.num = paired.df$no.of.duplication))
  }
}

## This function is written to paint tree (plotSimmap)
paint.tree<-function(tree)
{
  ## Considering real tree data
  gene.tree<-tree@phylo
  gene.data<-tree@data
  
  ## Identifiying duplication nodes and edges to paint them to simulate them at specified trait evolutionary rate
  dup.nodes <- gene.data$node[which(gene.data$events=="duplication")]
  dup.edges <- unique(gene.tree$edge[which(gene.tree$edge[,1] %in% dup.nodes), 2])
  tree.painted <- paintBranches (gene.tree, edge=dup.edges, "D", anc.state="S")
  return(tree.painted)
}

####This function helps to detect basal polytomy tree
### Although CAIC can handle polytomy, it cannot handle polytomy at root
basal.polytomy<-function(tree)
{
  node.root<-root.node(tree)
  daughter.root<-as.numeric(length(tree@phylo$edge[which(tree@phylo$edge[,1]==node.root),2]))
  if(daughter.root == 2){return(tree)}
  if(daughter.root >= 3){return(NA)}
}

## Computing average gene expression level and tissue specificity for our data
## For this we used 4 functions created by Nadezda

#Mean value per gene is calculated
fmean <- function(x)
{
  if(!all(is.na(x))) {
    res <- mean(x, na.rm=TRUE)
  } else {
    res <- NA
  }
  return(res)
}

##Function require a vector with expression of one gene in different tissues.
#Max value is calculated.
fmax<- function(x)
{
  if(!all(is.na(x))) {
    res <- max(x, na.rm=TRUE)
  } else {
    res <- NA
  }
  return(res)
}

##Function require a vector with expression of one gene in different tissues.
#Min value is calculated.
fmin <- function(x)
{
  if(!all(is.na(x))) {
    res <- min(x, na.rm=TRUE)
  } else {
    res <- NA
  }
  return(res)
}

#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
#Minimum 2 tissues
fTau <- function(x)
{
  if(all(!is.na(x))) {
    if(min(x, na.rm=TRUE) >= 0) {
      if(max(x)!=0) {
        x <- (1-(x/max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res/(length(x)-1)
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    } 
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  } 
  return(res)
}


## This function is written to match tree index
Info_collection_tree<-function(tree)
{
  ##Collecting gene tree data
  gene_tree<-tree@phylo
  gene_data<-tree@data
  ntips <-length(gene_tree$tip.label) ## Number of tips
  Internal_event_number<-as.numeric(length(gene_data$node[which(!is.na(gene_data$Tau))]))
  return(tibble(tip_num=ntips))
  
}

## This function is written to identify 1:1 orthologs based on species trees
ortholog.finder<-function(tree)
{
  ortholog_ID<-NULL
  
  ## Collecting tip labels
  tip_label<-tree@data$label[which(is.na(tree@data$events))]
  #gar_id<-tip_label[grepl("ENSLOCG",tip_label)]
  human_id<-tip_label[grepl("ENSG000",tip_label)]
  mouse_id<-tip_label[grepl("ENSMUSG",tip_label)]
  zfish_id<-tip_label[grepl("ENSDARG",tip_label)]
  medaka_id<-tip_label[grepl("ENSORLG00000",tip_label)]
  cavefish_id<-tip_label[grepl("ENSAMXG",tip_label)]
  pike_id<-tip_label[grepl("ENSELUG",tip_label)]
  tilapia_id<-tip_label[grepl("ENSONIG",tip_label)]
  rat_id<-tip_label[grepl("ENSRNOG",tip_label)]
  cow_id<-tip_label[grepl("ENSBTAG",tip_label)]
  macaque_id<-tip_label[grepl("ENSMMUG",tip_label)]
  #cod_id<-tip_label[grepl("ENSGMOG",tip_label)]
  chicken_id<-tip_label[grepl("ENSGALG",tip_label)]
  dog_id<-tip_label[grepl("ENSCAFG",tip_label)]
  rabbit_id<-tip_label[grepl("ENSOCUG",tip_label)]
  ferret_id<-tip_label[grepl("ENSMPUG",tip_label)]
  opossum_id<-tip_label[grepl("ENSMODG",tip_label)]
  
  ## Counting number
  freq<<-freq+1
  
  ## Getting ortholog groups ID
  ortholog_ID<-paste0("OG",freq,sep="")
  
  return(tibble(group=ortholog_ID,
                #spotted_gar.ID=gar_id,
                human.ID=human_id,
                mouse.ID=mouse_id,
                zebrafish.ID=zfish_id,
                medaka.ID=medaka_id,
                cavefish.ID=cavefish_id,
                northern_pike.ID=pike_id,
                tilapia.ID=tilapia_id,
                rat.ID=rat_id,
                cow.ID=cow_id,
                macaque.ID=macaque_id,
                chicken.ID=chicken_id,
                dog.ID=dog_id,
                rabbit.ID=rabbit_id,
                ferret.ID=ferret_id,
                opossum.ID=opossum_id
                #atlantic_cod.ID=cod_id
  ))
}

## This function is written to identify orthogroups for fish whole genome duplicates
orthogroup.finder<-function(tree)
{
  orthogroup_ID<-NULL
  
  ## Collecting tip labels
  tip_label<-tree@phylo$tip.label
  zfish_id<-tip_label[grepl("ENSDARG",tip_label)]
  medaka_id<-tip_label[grepl("ENSORLG",tip_label)]
  cavefish_id<-tip_label[grepl("ENSAMXG",tip_label)]
  pike_id<-tip_label[grepl("ENSELUG",tip_label)]
  tilapia_id<-tip_label[grepl("ENSONIG",tip_label)]
  
  ## Counting number
  freq<<-freq+1
  
  ## Getting ortholog groups ID
  orthogroup_ID<-paste0("OG",freq,sep="")
  
  ##Returning the tibble
  return(tibble(group=orthogroup_ID,
                zebrafish.1.ID=zfish_id[1],
                zebrafish.2.ID=zfish_id[2],
                medaka.1.ID=medaka_id[1],
                medaka.2.ID=medaka_id[2],
                cavefish.1.ID=cavefish_id[1],
                cavefish.2.ID=cavefish_id[2],
                northern_pike.1.ID=pike_id[1],
                northern_pike.2.ID=pike_id[2],
                tilapia.1.ID=tilapia_id[1],
                tilapia.2.ID=tilapia_id[2]))
}

## This function is written to collect expression levels of tissues of orthologs or orthogroups before performing PCA analysis
## We need to provide the dataframe of ortholog or orthogroup as input
## so we need to mention the type as ortholog or orthogroup 
expression.all.sample<-function(dataframe,type)
{
  tissues.of.interest<-c("brain","heart","kidney","liver","muscle","testis")
  hash.table=hash(keys=c("human","mouse","zebrafish","medaka","cavefish","northern_pike","tilapia","rat","cow","macaque","chicken","dog","rabbit","ferret","opossum"),
                  values=c("human.ID","mouse.ID","zebrafish.ID","medaka.ID","cavefish.ID","northern_pike.ID","tilapia.ID","rat.ID","cow.ID","macaque.ID","chicken.ID","dog.ID","rabbit.ID","ferret.ID","opossum.ID"))
  
  species<- sapply(colnames(dataframe)[-1],function(x){x <- unlist(strsplit(toString(x), split='.', fixed=TRUE))[1]})
  pca.data<-NULL
  start<-NULL
  end<-NULL
  for(a in 1:length(tissues.of.interest))
  {
    tissue<-NULL
    spec<-NULL
    new.data<-NULL
    final<-NULL
    tissue<-paste0("TPM.",tissues.of.interest[a],sep="")
    for(j in 1:length(species))
    {
      data<-NULL
      filename<-NULL
      spec<-species[j]
      
      ## For orthologs reading files and generating dataframe including all tissues
      if(type=="ortholog")
      {
        filename<-paste0(folder3,spec,"_tau.txt",sep="")
        data<-read.table(filename,header = T,sep = "\t")
        data1<-data.frame(ID=data[,1], data[tissue])
        colnames(data1)<-c(hash.table[[spec]],(paste0(spec,".",tissues.of.interest[a],sep="")))
        if(j==1){new.data<-dataframe}
        new.data<-merge(new.data,data1,by= colnames(data1)[1])
      }
      
      if(type=="orthogroup")
      {
        filename<-paste0(folder1,spec,"_normalized_exp.txt",sep="")
        data<-read.table(filename,header = T,sep = "\t")
        data1<-data.frame(ID=data[,1], data[tissue])
        colnames(data1)<-c(hash.table[[spec]],(paste0(spec,".",tissues.of.interest[a],sep="")))
        #if(j==1){colnames(data1)<-c(hash.table[[spec]],(paste0(spec,".",tissues.of.interest[a],sep="")))}
        if(j>=1 & j%%2 == 1){colnames(data1)<-c(paste0(spec,".1.ID",sep=""),(paste0(spec,".1.",tissues.of.interest[a],sep="")))}
        if(j>=1 & j%%2 == 0){colnames(data1)<-c(paste0(spec,".2.ID",sep=""),(paste0(spec,".2.",tissues.of.interest[a],sep="")))}
        if(j==1){new.data<-dataframe}
        new.data<-merge(new.data,data1,by= colnames(data1)[1])
      }
    }
    start<-as.numeric(length(colnames(dataframe)))
    end<-as.numeric(length(colnames(new.data)))
    
    ## collecting data for all samples
    final<-as.matrix(new.data[c((start+1):end)])
    row.names(final)<-new.data[,start]
    
    ## Collecting data and scaling factors for samples
    all.data <- data.frame(group=new.data[,start])
    all.data <-cbind(all.data,final)
    
    if(a==1){pca.data<-all.data}
    if(a>1){pca.data<-cbind(pca.data,all.data[,c(-1)])}
  }
  return(pca.data)
}

## Following the Brawand et al. (2011) method, this function is written to precisely derive scaling factor for samples and to normalize dataset
## First we aim to identify conserved genes with expression values in interqaurtile range
## Then we need to identify top conserved genes that have most conserved ranks amon samples
normalize.explevel.new<-function(dataset)
{
  columns<-as.numeric(ncol(dataset))
  
  ##Genes that are expressed in all samples
  #dataset.nonzero<-dataset[which(apply(dataset[,2:columns],1,min)>0 ),]
  
  ##Transforming the expression levels into ranks
  dataset.ranks<-apply(dataset[,2:columns],2,rank)
  
  ## Computing the variance of the ranks for each gene
  dataset.ranks.var<-apply(dataset.ranks,1,var,na.rm=T)
  
  ## Ranking the genes according to their variance
  dataset.consrank<-rank(dataset.ranks.var)
  
  ## Computing the median rank over all samples, for each gene
  median.rank<-apply(dataset.ranks,1,median,na.rm=T)
  
  ## Getting the genes which have a median rank in the 25%-75% range
  interquartile<-median.rank>(0.25*length(median.rank)) & median.rank<(0.75*length(median.rank))
  
  ##Getting the top 500 orthologous genes with the most conserved ranks, among those that fall in the interquartile range
  conservedgenes<-names(sort(dataset.consrank[interquartile]))[1:500]
  common.conserve.genes<-dataset[which(dataset$group %in% conservedgenes),]
  
  #print(nrow(common.conserve._genes))
  
  ## Calculating the normalization coefficient for each sample, i.e. the median of the TPM of the conserved orthologs in that sample
  sample.median<-apply(common.conserve.genes[,2:columns],2,median)
  
  ## Bringing all of the medians at the average of the medians 
  #normcoeff<-sample_median/mean(sample_median) 
  #normcoeff<-sample_median/mean(apply(dataset.nonzero[,2:columns],2,median)) ## Normalizing factor
  
  ## Matrix to be normalized
  final<-as.matrix(dataset[,2:columns])
  row.names(final)<-dataset[,1]
  
  ## Deriving scaling factor for each sample/tissue
  scale<-NULL
  tpm.colnames <-sapply(colnames(dataset[,2:columns]),
                        function(x) 
                        {x <- unlist(strsplit(toString(x), split='.', fixed=TRUE))[2]})
  tissuesNames<-as.vector(unique(tpm.colnames))
  end<-0
  for(i in 1:length(tissuesNames))
  {
    factor<-NULL
    start<-end+1
    factor<-sample.median[c(start:(start+14))]
    
    ## Collecting scaling factor for all tissue 
    scale<-append(scale,sapply(seq(1, length(factor)),function(x) factor[x] / factor[1]))
    end<-start+14
  }
  
  ## Using the scaling factor normalizing the ortholog set
  data.norm <- sapply(seq(1, length(scale)), function(x) final[,x] / scale[x])
  
  ## Collecting normalized ortholog set
  normalized.pca.ortholog<-data.frame(dataset[,1])
  normalized.pca.ortholog<-cbind(normalized.pca.ortholog,data.norm)
  colnames(normalized.pca.ortholog)<-colnames(dataset)
  
  ## Saving the species specific scaling factor (which can be used to normalize expression of the tissue of the species) & returning the normalized data for ortholog set
  scale<-data.frame(scale)
  colnames(scale)<-c("scaling.factor")
  outfile<- paste0(folder1,"/OG.scaling_factor.txt",sep="")
  write.table(scale,outfile, row.names = T, sep="\t",quote = F)
  return(normalized.pca.ortholog)
}


## This function is reads the scored scaling factor for the samples/tissues in each species from above function, and generating normalize dataset for each species
normalized.species.data<-function(scaling.factor)
{
  ##Initializing variable
  data.final<-NULL
  
  ## Reading scaling factors for species 
  input<-read.table(paste0(folder1,"/OG.scaling_factor.txt",sep=""),header = T,sep="\t")
  input$row.names<-row.names(input)
  row.names(input)<-NULL
  species<- unique(sapply(input$row.names,function(x){x <- unlist(strsplit(toString(x), split='.', fixed=TRUE))[1]}))
  
  ## Reading species data table to normalize
  for(i in 1:length(species))
  {
    spec<-NULL
    sf<-NULL
    data1<-NULL
    data.norm<-NULL
    spec<-species[i]
    
    ## Getting tissue scaling factor for the species
    sf<-input$scaling.factor[grepl(spec,input$row.names)] 
    names(sf)<-input$row.names[grepl(spec,input$row.names)]
    
    ## Reading raw expression data
    data<-read.table(paste0(folder3,spec,"_tau.txt",sep=""),header = T,sep="\t")
    
    ## Removing last column with tau data
    data<-data[c(1:(ncol(data)-1))]
    
    ## Removing column data for liver
    #drops <- c("TPM.liver")
    #data<-data[ , !(names(data) %in% drops)]
    
    ## Arranging tissues in proper order before normalization
    tissues<- sapply(names(sf),function(x){x <- unlist(strsplit(toString(x), split='.', fixed=TRUE))[2]})
    tissuef<-paste0("TPM.",tissues,sep="")
    data1<-data[,1]
    data1<-cbind(data1,data[tissuef])
    
    ## matrix to be normalized
    final<-as.matrix(data1[,2:(ncol(data1))])
    row.names(final)<-data1[,1]
    
    ## Using the scaling factor normalizing species expression dataset
    data.norm <- sapply(seq(1, length(sf)), function(x) final[,x] / sf[x])
    normalized_data<-data.frame(Ens.ID=data[,1])
    normalized_data<-cbind(normalized_data,data.norm)
    colnames(normalized_data)<-colnames(data1)
    
    ## Saving the normalized dataset for each species
    outfile<- paste0(folder1,spec,"_normalized_exp.txt",sep="")
    write.table(normalized_data,outfile, row.names = F, sep="\t",quote = F)
    data.final<-rbind(data.final,normalized_data)
  }
  row.names(data.final)<-NULL
  return(data.final)
}


## this function is written to identify highly and lowly expressed WGDs
duplicate.classfication1<-function(orthogroup,dataframe)
{
  ## We required IDs of orthogroup before classification, so we merge it
  pca.WGD1<-NULL
  pca.WGD1<-merge(orthogroup,dataframe,by="group")
  
  ## Creating new dataframe
  pca_WGD_highly_lowly_expressed_duplicates<-NULL
  pca_WGD_highly_lowly_expressed_duplicates<-pca.WGD1[c(1)]
  
  ## selecting species with duplicates
  spec <- unique(sapply(colnames(dataframe),function(x){x <- unlist(strsplit(toString(x), split='.', fixed=TRUE))[1]}))
  spec <- spec[-1]
  
  ## Now calculating average normalized gene expression level per gene to classify genes
  for(i in 1:length(spec))
  {
    ## Initializing variables
    species<-NULL
    species_data1<-NULL
    species_data2<-NULL
    var1<-NULL
    var2<-NULL
    
    ## Collecting data and after classification collecing the dataframe
    species<-spec[i]
    var1<-paste0(species,".high.exp.dup",sep="")
    var2<-paste0(species,".low.exp.dup",sep="")
    species_data1<-data.frame(pca.WGD1[,regexpr(paste0(species,".1",sep=""),colnames(pca.WGD1))>0])
    species_data1$avg.exp1<-apply(species_data1[,2:ncol(species_data1)],1,fmean)
    species_data1<-species_data1[c(1,ncol(species_data1))]
    species_data2<-data.frame(pca.WGD1[,regexpr(paste0(species,".2",sep=""),colnames(pca.WGD1))>0])
    species_data2$avg.exp2<-apply(species_data2[,2:ncol(species_data2)],1,fmean)
    species_data2<-species_data2[c(1,ncol(species_data2))]
    
    ## Classifying highly expressed and lowly expressed duplicates
    species_data1<-cbind(species_data1,species_data2)
    species_data1$High.Exp.copy<-ifelse(species_data1$avg.exp1>=species_data1$avg.exp2, species_data1[,1],species_data1[,3])
    species_data1$Low.Exp.copy<-ifelse(species_data1$avg.exp1<species_data1$avg.exp2, species_data1[,1],species_data1[,3])
    
    ## Adding information to the orthogroup dataframe
    pca_WGD_highly_lowly_expressed_duplicates[[var1]]<-species_data1$High.Exp.copy
    pca_WGD_highly_lowly_expressed_duplicates[[var2]]<-species_data1$Low.Exp.copy
  }
  return(pca_WGD_highly_lowly_expressed_duplicates)
}


## This function is written to generate PCA plot
PCA.plot <- function(dataframe,type,group.by,method)
{
  ## removing rows with no variance
  total.column<-ncol(dataframe)
  pca1 <- dataframe[apply(dataframe[,c(2:total.column)], 1, var) != 0, ]
  qmatrix<-as.matrix(pca1[,c(2:total.column)])
  tpm.colnames <-NULL
  tissuesNames<-NULL
  
  ## datatype
  if(type=="ortholog")
  {
    number.tissues<-6 # change here if needed
    number.species<-15
    ## identifying tissue names
    tpm.colnames <-sapply(colnames(dataframe[,2:total.column]),
                          function(x) 
                          {x <- unlist(strsplit(toString(x), split='.', fixed=TRUE))[2]})
    tissuesNames<-as.vector(tpm.colnames)
    sex<-c("F","M","F","F","F","F","F","M","M","M","M","M","M","M","M",
          "M","M","Mix","F","F","F","M","M","M","M","M","M","M","M","M",
          "F","M","Mix","F","F","F","F","M","M","M","M","M","M","M","M",
          "M","M","F","F","F","F","F","M","M","M","M","M","M","M","M",
          "M","M","F","F","F","F","F","M","M","M","M","M","M","M","M",
          "M","M","M","M","M","M","M","M","M","M","M","M","M","M","M")
    
    ## "group1" is required for grouping by "species" & "group2" is required for grouping by "tissues"
    group1<- rep(c("human","mouse","zebrafish","medaka","cavefish","pike","tilapia","rat","cow","macaque","chicken","dog","rabbit","ferret","opossum"),times=number.tissues)
    group2<-rep(tissuesNames,times=1)
    group3<-rep(sex,times=1)
    #pch<-rep(c(16,16,15,15,15,15,15,16,16,16,16,16,16,16,16), times=number.tissues)
    
    tissue.color<-c(rep("#0076b6",times=number.species),rep("#ff0000",times=number.species),
                    rep("#b5651d",times=number.species),rep("#008000",times=number.species),
                    rep("#CA2C92",times=number.species),rep("#ff5f1f",times=number.species))
                    #rep("#0091b3",times=number.species),
                    #rep("#dd6e6e",times=number.species),rep("#f7d031",times=number.species),
                    #rep("#fdae6b",times=number.species),rep("#fa9fb5",times=number.species) # change here if needed
    names(tissue.color)<-tissuesNames
    species<-c("human","mouse","zebrafish","medaka","cavefish","pike","tilapia","rat","cow","macaque","chicken","dog","rabbit","ferret","opossum") 
    species.color<-rep(c("#e31a1c","#a8ddb5","#e78ac3","#666666","#ff7f00","#386cb0","#c12828","#fc9272","#1c9099","#636363","#d8b365","#e9a3c9","#a1d76a","#fc8d59","#ffffb3"),times=number.tissues)
    names(species.color)<-rep(species,times=number.tissues)
    species.sign<-rep(c(4,5,6,15,16,17,18,7,8,9,10,11,12,13,14),times=number.tissues)
    sex.color<-c("#dd1c77","#3182bd","#dd1c77","#dd1c77","#dd1c77","#dd1c77","#dd1c77","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd",
                 "#3182bd","#3182bd","#31a354","#dd1c77","#dd1c77","#dd1c77","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd",
                 "#dd1c77","#3182bd","#31a354","#dd1c77","#dd1c77","#dd1c77","#dd1c77","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd",
                 "#3182bd","#3182bd","#dd1c77","#dd1c77","#dd1c77","#dd1c77","#dd1c77","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd",
                 "#3182bd","#3182bd","#dd1c77","#dd1c77","#dd1c77","#dd1c77","#dd1c77","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd",
                 "#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd","#3182bd")
    #names(sex.color)<-c("F","M","Mix")
    names(sex.color)<-sex
  }
  
  if(type=="common.WGD")
  {
    number.tissues<-6 ## Change here if needed
    number.species<-5 ## Change here if needed
    species.with.duplication<-5
    
    ## Identifying tissuenames
    tissuesNames<-c(rep(c("brain1","brain2"),times=number.species),
                    rep(c("heart1","heart2"),times=number.species),
                    rep(c("kidney1","kidney2"),times=number.species),
                    rep(c("liver1","liver2"),times=number.species),
                    rep(c("muscle1","muscle2"),times=number.species),
                    rep(c("testis1","testis2"),times=number.species))
    
    ## "group1" is required for grouping by "species" & "group2" is required for grouping by "tissues"
    group1<- rep(c("zebrafish","zebrafish","medaka","medaka","cavefish","cavefish","pike","pike","tilapia","tilapia"),times=number.tissues)
    group2<-rep(tissuesNames,times=1)
    
    tissue.color<-c(rep(c("#0076b6","#00f9ff"),times=species.with.duplication),
                    rep(c("#ff0000","#ff817e"),times=species.with.duplication),
                    rep(c("#b5651d","#c19a6b"),times=species.with.duplication),
                    rep(c("#008000","#00FF00"),times=species.with.duplication),
                    rep(c("#CA2C92","#ff00ff"),times=species.with.duplication),
                    rep(c("#ff5f1f","#ffd7b5"),times=species.with.duplication)) ##change here if needed
    
    names(tissue.color)<-tissuesNames
    species<-c("zebrafish","zebrafish","medaka","medaka","cavefish","cavefish","pike","pike","tilapia","tilapia")
    species.color<-rep(c("#e31a1c","#e31a1c","#e78ac3","#e78ac3","#666666","#666666","#ff7f00","#ff7f00","#386cb0","#386cb0"),times=number.tissues)
    names(species.color)<-rep(species,times=number.tissues)
    species.sign<-rep(c(6,6,15,15,16,16,17,17,18,18),times=number.tissues)
  }
  ## Performing Principle Component Analysis
  pca.data<-t(qmatrix) 
  pca <- prcomp(pca.data,scale = T, center = T)
  
  xlab=paste0("PC1: ", round(summary(pca)$importance[2,1],1)*100, "% variance explained")
  ylab=paste0("PC2: ", round(summary(pca)$importance[2,2],2)*100, "% variance explained")
  
  if(group.by=="tissue" & method=="normal.2D")
  {
    if(type=="ortholog")
    {
      # Add extra space to right of plot area; change clipping to figure
      par(mar=c(4, 4, 4, 8.1), xpd=TRUE) 
      plot(pca$x[,1], pca$x[,2], xlab=xlab, ylab=ylab,pch=species.sign,col=rep(tissue.color,times=1), font.lab=2,cex.lab=0.8, cex=1)
      legend("topleft", inset = c(0,-0.15), legend=as.factor(unique(unlist(tissuesNames))), ncol = 1, cex=0.7,pch=19, col=unique(tissue.color),text.font=2, bty="n", horiz=T)
      legend("topright", inset = c(-0.5,0), legend=as.factor(unique(species)), cex=0.7,pch=c(4,5,6,15,16,17,18,7,8,9,10,11,12,13,14), col ="black",text.font=2, bty="n", ncol=2)
    }
    if(type=="common.WGD")
    {
      # Add extra space to right of plot area; change clipping to figure
      par(mar=c(4, 4, 4, 8.1), xpd=TRUE) 
      plot(pca$x[,1], pca$x[,2], xlab=xlab, ylab=ylab,pch=species.sign,col=rep(tissue.color,times=1), font.lab=2,cex.lab=0.8, cex=1)
      legend("topright", inset = c(-0.5,0.2), legend=c("brain","heart","kidney","liver","muscle","testis","brain","heart","kidney","liver","muscle","testis"),cex=0.7,pch=19,col=c("#0076b6","#ff0000","#b5651d","#008000","#CA2C92","#ff5f1f","#00f9ff","#ff817e","#c19a6b","#00FF00","#ff00ff","#ff9248"), text.font=2, bty="n",ncol=2,title="High exp  Low exp")
      legend("topleft", inset = c(0,-0.15), legend=as.factor(unique(species)), cex=0.7,pch=c(6,15,16,17,18), col ="black",text.font=2, bty="n", horiz=T)
    }
  }
  if(group.by=="sex" & method=="normal.2D")
  {
    
    plot(pca$x[,3], pca$x[,4], xlab=xlab, ylab=ylab,pch=species.sign,col=rep(sex.color,times=1),cex=1.2,mar=c(5,2,6,2))
    par(xpd=T)
    legend("top", inset = c(0.03,-0.10), legend=as.factor(unique(unlist(sex))), ncol = 1, cex=0.6,pch=19, col=unique(sex.color), bty="y",horiz=T)
    legend("top", inset = c(0.03,-0.20), legend=as.factor(unique(species)), ncol = 1, cex=0.6,pch=c(4,5,6,15,16,17,18,7,8,9,10,11,12,13,14), col ="black", bty="y",horiz=T)
  }
  
  if(group.by=="species" & method=="normal.2D")
  {
    plot(pca$x[,3], pca$x[,4], xlab=xlab, ylab=ylab,pch=species.sign,col=rep(species.color,times=1),cex=1.2,mar=c(5,2,6,2))
    par(xpd=T)
    legend("top", inset = c(0.02,-0.10), legend=as.factor(unique(unlist(tissuesNames))), ncol = 1, cex=0.8,pch=19, col=unique(tissue.color), bty="y",horiz=T)
    legend("top", inset = c(0.02,-0.20), legend=as.factor(unique(species)), ncol = 1, cex=0.8,pch=c(4,5,6,15,16,17,18,7,8,9,10,11,12,13,14), col ="black", bty="y")
  }
}
  
## Function to create boxplot result of our data
boxplot.new<-function(dataframe,pval,trait,med)
{
  dodge <- position_dodge(width = 0.41)
  
  if(trait=="Tau")
  {
    plot<-ggplot(dataframe,aes(x=Event, y=Tau.abs, fill=Event)) + 
      #geom_violin(position = dodge, alpha=0.3)+
      guides(colour = guide_legend(override.aes = list(shape = 16))) +
      geom_boxplot(width=0.5, position = dodge, outlier.colour=NA,notch = T) +
      geom_text(data = med, size=3.5, aes(label = pic.round),fontface = 2,hjust=0.5,vjust =-1)+
      xlab( NULL ) +
      ylab(expression(bold(paste("PICs of ",tau)))) +
      theme_classic()+
      theme(legend.title=element_blank(),legend.position="none") +
      theme(axis.text=element_text(size=10,face="bold")) +
      theme(legend.text=element_text(size=10,face="bold"))+
      annotate("rect", xmin = 1.1, xmax = 1.9, ymin = 0.03, ymax =0.03, alpha=0.8,colour = "blue")+
      annotate("text", x = 1.5, y = 0.033, label= as.character(paste0("italic(",pval,")")), parse=TRUE, size=4)+ #using plotmath syntax with parse = TRUE
      annotate("text", x = 1.5, y = 0.033, label= as.character(paste0("italic(",pval,")")), parse=TRUE, size=4.02)+ 
      annotate("text", x = 1.5, y = 0.033, label= as.character(paste0("italic(",pval,")")), parse=TRUE, size=4.04) 
  }
  
  if(trait=="Exp")
  {
    plot<-ggplot(dataframe,aes(x=Event, y=Exp.abs, fill=Event)) + 
      guides(colour = guide_legend(override.aes = list(shape = 16))) +
     #geom_boxplot(width=0.5,outlier.colour=NA, position = dodge, notch = T) +
      geom_boxplot(width=0.5,position = dodge, outlier.colour=NA,notch = T) +
      geom_text(data = med, size=3.5, aes(label = pic.round),fontface = 2,hjust=0.5,vjust =-1)+
      xlab( NULL ) +
      ylab(expression(bold(paste("PICs of exp levels")))) +
      theme_classic()+
      theme(legend.title=element_blank(),legend.position="none") +
      theme(axis.text=element_text(size=10,face="bold")) +
      theme(legend.text=element_text(size=10,face="bold"))+
      annotate("rect", xmin = 1.1, xmax = 1.9, ymin = 0.18, ymax =0.18, alpha=0.8,colour = "blue")+
      #annotate("text", x = 1.5, y = 0.2, label= pval, parse=TRUE)
      annotate("text", x = 1.5, y = 0.2, label= as.character(paste0("italic(",pval,")")), parse=TRUE, size=4)+ 
      annotate("text", x = 1.5, y = 0.2, label= as.character(paste0("italic(",pval,")")), parse=TRUE, size=4.02)+ 
      annotate("text", x = 1.5, y = 0.2, label= as.character(paste0("italic(",pval,")")), parse=TRUE, size=4.04) 
  }
  return(plot)
}

boxplot.new.3RWGD<-function(dataframe,pval1,pval2,trait,med)
{
  dodge <- position_dodge(width = 0.41)
  if(trait=="Tau.3RWGD")
  {
    plot<-ggplot(dataframe,aes(x=Event, y=Tau.abs, fill=Event)) + 
    guides(colour = guide_legend(override.aes = list(shape = 16))) +
    geom_boxplot(width=0.5, position = dodge, outlier.colour=NA,notch = T) +
    geom_text(data = med,size=3.2, aes(label = pic.round),fontface = 2,hjust=0.5,vjust =-1)+
    scale_fill_manual(values=c("#F8766D", "#00BFC4", "#008000"))+
    xlab( NULL ) +
    ylab(expression(bold(paste("PICs of ",tau)))) +
    theme_classic()+
    theme(legend.title=element_blank(),legend.position=c(9.9,9.9)) +
    theme(axis.text=element_text(size=10,face="bold")) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
    theme(legend.text=element_text(size=10,face="bold"))+
    annotate("rect", xmin = 1.1, xmax = 1.9, ymin = 0.021, ymax =0.021, alpha=0.8,colour = "blue")+
    #annotate("text", x = 1.5, y = 0.024, label= pval2, fontface = 4,parse=TRUE)+
    annotate("text", x = 1.5, y = 0.024, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4)+ 
    annotate("text", x = 1.5, y = 0.024, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4.02)+ 
    annotate("text", x = 1.5, y = 0.024, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4.04)+ 
    annotate("rect", xmin = 1.1, xmax = 2.9, ymin = 0.03, ymax =0.03, alpha=0.8,colour = "blue")+
    annotate("text", x = 2, y = 0.033, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4)+
    annotate("text", x = 2, y = 0.033, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4.02)+
    annotate("text", x = 2, y = 0.033, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4.04)
  }
  if(trait=="Exp.3RWGD")
  {
    plot<-ggplot(dataframe,aes(x=Event, y=Exp.abs, fill=Event)) + 
      guides(colour = guide_legend(override.aes = list(shape = 16))) +
      geom_boxplot(width=0.5,outlier.colour=NA, position = dodge, notch = T) +
      geom_text(data = med, size=3.2,aes(label = pic.round),fontface = 2,hjust=0.5,vjust =-1)+
      scale_fill_manual(values=c("#F8766D", "#00BFC4", "#008000"))+
      xlab( NULL ) +
      ylab(expression(bold(paste("PICs of expression levels")))) +
      theme_classic()+
      theme(legend.title=element_blank(),legend.position=c(9.9,9.9)) +
      theme(legend.position="none")+
      theme(axis.text=element_text(size=10,face="bold")) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
      theme(legend.text=element_text(size=10,face="bold")) +
      annotate("rect", xmin = 1.1, xmax = 1.9, ymin = 0.13, ymax =0.13, alpha=0.8,colour = "blue")+
      #annotate("text", x = 1.5, y = 0.15, label= pval2, fontface = 4,parse=TRUE)+
      annotate("text", x = 1.5, y = 0.15, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4)+ 
      annotate("text", x = 1.5, y = 0.15, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4.02)+ 
      annotate("text", x = 1.5, y = 0.15, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4.04)+ 
      annotate("rect", xmin = 1.1, xmax = 2.9, ymin = 0.18, ymax =0.18, alpha=0.8,colour = "blue")+
      #annotate("text", x = 2.0, y = 0.2, label= pval1, fontface = 4,parse=TRUE)
      annotate("text", x = 2.0, y = 0.2, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4) 
      annotate("text", x = 2.0, y = 0.2, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4.02) 
      annotate("text", x = 2.0, y = 0.2, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4.04) 
  }
  return(plot)
}

## Function to create violin plot for FishWGD trees with jitter 
vioplot.new.3RWGD2<-function(dataframe,pval1,pval2,trait,med)
{
  dodge <- position_dodge(width = 0.21)
  if(trait=="Tau.3RWGD")
  {
    plot<-ggplot(dataframe,aes(x=Event, y=Tau.abs, shape=Event,col=Event)) + 
      geom_violin(position = dodge,alpha=0.2)+
      geom_jitter(position=position_jitter(0.05),alpha=0.2,size=1)+
      scale_color_manual(values=c("#F8766D", "#00BFC4", "#008000"))+
      geom_text(data = med, aes(label = pic.round),fontface = 2,size = 3.2, 
                hjust=0.3,vjust =-2.5,col= "black")+
      geom_crossbar(data=med, aes(ymin = Tau.abs, ymax = Tau.abs),
                    size=0.7,col= "blue", width = 0.15)+
      xlab( NULL ) +
      ylab(expression(bold(paste("PICs of ",tau)))) +
      theme_classic()+
      theme(legend.title=element_blank(),legend.position=c(9.9,9.9)) +
      theme(axis.text=element_text(size=10,face="bold")) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
      theme(legend.text=element_text(size=10,face="bold"))
  }
  if(trait=="Exp.3RWGD")
  {
    plot<-ggplot(dataframe,aes(x=Event, y=Exp.abs, shape=Event,col=Event)) + 
      geom_violin(position = dodge,alpha=0.2)+
      geom_jitter(position=position_jitter(0.05),alpha=0.2,size=1)+
      scale_color_manual(values=c("#F8766D", "#00BFC4", "#008000"))+
      geom_text(data = med, aes(label = pic.round),fontface = 2,size = 3.2, 
                hjust=0.3,vjust =-2.5,col= "black")+
      geom_crossbar(data=med, aes(ymin = Exp.abs, ymax = Exp.abs),
                    size=0.7,col= "blue", width = 0.15)+
      xlab( NULL ) +
      ylab(expression(bold(paste("PICs of expression levels")))) +
      theme_classic()+
      theme(legend.title=element_blank(),legend.position=c(9.9,9.9)) +
      theme(legend.position="none")+
      theme(axis.text=element_text(size=10,face="bold")) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
      theme(legend.text=element_text(size=10,face="bold")) 
  }
  return(plot)
}

## Function to create violin plot for SSD trees  with jitter 
vioplot.new2<-function(dataframe,pval,trait,med)
{
  dodge <- position_dodge(width = 0.21)
  
  if(trait=="Tau")
  {
    plot<-ggplot(dataframe,aes(x=Event, y=Tau.abs, shape=Event,col=Event)) + 
      geom_violin(position = dodge,alpha=0.2)+
      #geom_area(position = dodge, alpha=1)+
      geom_jitter(position=position_jitter(0.05),alpha=0.2,size=1)+
      scale_color_manual(values=c("#F8766D", "#00BFC4"))+
      geom_text(data = med, aes(label = pic.round),fontface = 2,size = 3.2, 
                hjust=0.3,vjust =-2.5,col= "black")+
      geom_crossbar(data=med, aes(ymin = Tau.abs, ymax = Tau.abs),
                    size=0.7,col= "blue", width = 0.15)+
      xlab( NULL ) +
      ylab(expression(bold(paste("PICs of ",tau)))) +
      theme_classic()+
      theme(legend.title=element_blank(),legend.position="none") +
      theme(axis.text=element_text(size=10,face="bold")) +
      theme(legend.text=element_text(size=10,face="bold"))
  }
  
  if(trait=="Exp")
  {
    plot<-ggplot(dataframe,aes(x=Event, y=Exp.abs, shape=Event,col=Event)) + 
      geom_violin(position = dodge,alpha=0.2)+
      #geom_area(position = dodge, alpha=1)+
      geom_jitter(position=position_jitter(0.05),alpha=0.2,size=1)+
      scale_color_manual(values=c("#F8766D", "#00BFC4"))+
      geom_text(data = med, aes(label = pic.round),fontface = 2,size = 3.2, 
                hjust=0.3,vjust =-2.5,col= "black")+
      geom_crossbar(data=med, aes(ymin = Exp.abs, ymax = Exp.abs),
                    size=0.7,col= "blue", width = 0.15)+
      xlab( NULL ) +
      ylab(expression(bold(paste("PICs of exp levels")))) +
      theme_classic()+
      theme(legend.title=element_blank(),legend.position="none") +
      theme(axis.text=element_text(size=10,face="bold")) +
      theme(legend.text=element_text(size=10,face="bold"))
  }
  return(plot)
}


## Boxplot function for our improved PIC  method 
boxplot.our.method<-function(dataframe,pval1,pval2,trait)
{
  dodge <- position_dodge(width = 0.41)
  
  if(trait=="Tau")
  {
    plot<-ggplot(dataframe,aes(x=Event, y=pic,fill=Event))+   
      guides(colour = guide_legend(override.aes = list(shape = 16))) +
      geom_boxplot(width=0.3,outlier.colour=NA, position = dodge, notch = T,stat="boxplot") +
      geom_violin(position = dodge, alpha=0.3)+
      xlab( NULL ) +
      ylab(expression(bold(paste("PICs of ", tau)))) +
      coord_cartesian(ylim=c(0, 0.035)) +
      theme_classic()+
      theme(legend.title=element_blank(),legend.position=c(9.9,9.9)) +
      theme(legend.position="none")+
      theme(axis.text=element_text(size=9,face="bold")) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
      theme(legend.text=element_text(size=10,face="bold")) +
      scale_fill_manual(values=c("#F8766D", "#00BFC4","#F8766D"))+
      annotate("rect", xmin = 1.1, xmax = 1.8, ymin = 0.022, ymax =0.022, alpha=0.8,colour = "blue")+
      annotate("text", x = 1.5, y = 0.025, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4)+ 
      annotate("text", x = 1.5, y = 0.025, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4.02)+
      annotate("text", x = 1.5, y = 0.025, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4.04)+ 
      annotate("rect", xmin = 1.2, xmax = 2.8, ymin = 0.031, ymax =0.031, alpha=0.8,colour = "blue")+
      annotate("text", x = 2, y = 0.034, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4)+ 
      annotate("text", x = 2, y = 0.034, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4.02)+ 
      annotate("text", x = 2, y = 0.034, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4.04) 
  }
  if(trait=="Exp")
  {
    plot<-ggplot(dataframe,aes(x=Event, y=pic, fill=Event)) + 
      guides(colour = guide_legend( override.aes = list(shape = 16))) +
      geom_boxplot(width=0.3,outlier.colour=NA, position = dodge, notch = T) +
      geom_violin(position = dodge, alpha=0.3)+
      xlab( NULL ) +
      ylab(expression(bold(paste("PICs of expression levels")))) +
      coord_cartesian(ylim=c(0, 0.25)) +
      theme_classic()+
      theme(legend.title=element_blank(),legend.position=c(9.9,9.9)) +
      theme(legend.position="none")+
      theme(axis.text=element_text(size=9,face="bold")) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
      theme(legend.text=element_text(size=10,face="bold")) +
      scale_fill_manual(values=c("#F8766D", "#00BFC4","#F8766D"))+
      annotate("rect", xmin = 1.1, xmax = 1.8, ymin = 0.16, ymax =0.16, alpha=0.8,colour = "blue")+
      annotate("text", x = 1.5, y = 0.18, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4)+ 
      annotate("text", x = 1.5, y = 0.18, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4.02)+ 
      annotate("text", x = 1.5, y = 0.18, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4.04)+ 
      annotate("rect", xmin = 1.2, xmax = 2.8, ymin = 0.22, ymax =0.22, alpha=0.8,colour = "blue")+
      annotate("text", x = 2.1, y = 0.24, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4)+ 
      annotate("text", x = 2.1, y = 0.24, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4.02)+ 
      annotate("text", x = 2.1, y = 0.24, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4.04) 
  }
  return(plot)
}

## Boxplot function for our improved PIC method for 3RWGD 
boxplot.our.method.3RWGD<-function(dataframe,pval1,pval2,trait)
{
  dodge <- position_dodge(width = 0.41)
  
  if(trait=="Tau")
  {
    plot<-ggplot(dataframe,aes(x=Event, y=pic,fill=Event))+   
      guides(colour = guide_legend(override.aes = list(shape = 16))) +
      geom_boxplot(width=0.3,outlier.colour=NA, position = dodge, notch = T,stat="boxplot") +
      geom_violin(position = dodge, alpha=0.3)+
      xlab( NULL ) +
      ylab(expression(bold(paste("PICs of ", tau)))) +
      coord_cartesian(ylim=c(0, 0.035)) +
      theme_classic()+
      theme(legend.title=element_blank(),legend.position=c(9.9,9.9)) +
      theme(legend.position="none")+
      theme(axis.text=element_text(size=9,face="bold")) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
      theme(legend.text=element_text(size=10,face="bold")) +
      scale_fill_manual(values=c("#F8766D", "#008000","#F8766D"))+
      annotate("rect", xmin = 1.1, xmax = 1.8, ymin = 0.022, ymax =0.022, alpha=0.8,colour = "blue")+
      annotate("text", x = 1.5, y = 0.025, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4)+ 
      annotate("text", x = 1.5, y = 0.025, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4.02)+ 
      annotate("text", x = 1.5, y = 0.025, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4.04)+ 
      annotate("rect", xmin = 1.2, xmax = 2.8, ymin = 0.031, ymax =0.031, alpha=0.8,colour = "blue")+
      annotate("text", x = 2, y = 0.034, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4)+ 
      annotate("text", x = 2, y = 0.034, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4.02)+ 
      annotate("text", x = 2, y = 0.034, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4.04) 
  }
  if(trait=="Exp")
  {
    plot<-ggplot(dataframe,aes(x=Event, y=pic, fill=Event)) + 
      guides(colour = guide_legend( override.aes = list(shape = 16))) +
      geom_boxplot(width=0.3,outlier.colour=NA, position = dodge, notch = T) +
      geom_violin(position = dodge, alpha=0.3)+
      xlab( NULL ) +
      ylab(expression(bold(paste("PICs of expression levels")))) +
      coord_cartesian(ylim=c(0, 0.25)) +
      theme_classic()+
      theme(legend.title=element_blank(),legend.position=c(9.9,9.9)) +
      theme(legend.position="none")+
      theme(axis.text=element_text(size=9,face="bold")) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
      theme(legend.text=element_text(size=10,face="bold")) +
      scale_fill_manual(values=c("#F8766D", "#008000","#F8766D"))+
      annotate("rect", xmin = 1.1, xmax = 1.8, ymin = 0.16, ymax =0.16, alpha=0.8,colour = "blue")+
      annotate("text", x = 1.5, y = 0.18, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4)+ 
      annotate("text", x = 1.5, y = 0.18, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4.02)+ 
      annotate("text", x = 1.5, y = 0.18, label= as.character(paste0("italic(",pval1,")")), parse=TRUE, size=4.04)+ 
      annotate("rect", xmin = 1.2, xmax = 2.8, ymin = 0.22, ymax =0.22, alpha=0.8,colour = "blue")+
      annotate("text", x = 2.1, y = 0.24, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4)+ 
      annotate("text", x = 2.1, y = 0.24, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4.02)+ 
      annotate("text", x = 2.1, y = 0.24, label= as.character(paste0("italic(",pval2,")")), parse=TRUE, size=4.04) 
  }
  return(plot)
}

## This function is written to map gene name to orthogroup index
group2name<-function(tree)
{
  Count <<- Count + 1
  print(Count)
  OGroup<-paste0("OG",Count,sep="")
  gene_tree<-tree@phylo
  gene_data<-tree@data
  tiplabels<-gene_data$label[which(is.na(gene_data$events))]
  gar_id<-tiplabels[grepl("ENSLOCG",tiplabels)]
  
  ##reading gene name
  gar_name<-read.table("gar_name.txt",header = T,sep="\t")
  
  df<-NULL
  df<-data.frame(Gene.stable.ID=gar_id)
  df<-merge(df,gar_name,by=c("Gene.stable.ID"))
  
  return(tibble(OG=OGroup,Gene=df$Gene.name))
}

## This function to check contrasts independence of internal parameters
diagnostic.plot.test<- function(tree,nodelimit)
{
  gene.tree<-tree@phylo
  gene.tree$node.label<-NULL
  count<<-count+1
  print(count)
  
  ## Collecting data
  data<-tree@data
  data.new<-data[which(!is.na(data$Tau)),] 
  data.new1<-data.frame(label=data.new$label,Tau=data.new$Tau, Mean.Exp=data.new$avg.exp)
  rownames(data.new1)<-data.new1$label
  treedata<-comparative.data(gene.tree,data.new1,label,vcv=T)
  test<-tryCatch(crunch(Mean.Exp~Tau,data=treedata,equal.branch.length=F), error = function(e) {"Error"})
  if((test=="Error") || (is.na(test$mod[1]$coefficients[1])))
  { print ("Error obtained in crunch")
    return(NA)
  }
  if(test!="Error")
  {
    ## We need to check for diagnostic plots for expected variance and node age from this Caper package
    contrast<-caic.table(test)
    contrast$node<-as.numeric(rownames(contrast))
    if(nodelimit==0){contrast$height<-data$nodeheight[which(data$node.age>0)]}
    if(nodelimit >0){contrast$height<-data$nodeheight[which(!is.na(data$events))]}
    
    ## Considering contrasts data without any age limit
    if(nodelimit==0){contrast.new<-contrast}
    if(nodelimit>0){contrast.new<-contrast[contrast$nodeAge<=nodelimit,]}
    
    ## length of the dataframe should be higher than 3
    if(nrow(contrast.new) >= 3)
    {
      ## Diagnstic tests for expected variance, node age, node depth and node height as recommended in earlier studies
      test.expvar<-summary(lm((abs(contrast.new$Tau))~sqrt(contrast.new$contrVar)))
      test.age<-summary(lm((abs(contrast.new$Tau))~log(contrast.new$nodeAge)))
      test.depth<-summary(lm((abs(contrast.new$Tau))~contrast.new$nodeDepth))
      
      if((!(is.na(test.expvar$coefficients[2,4]))) || (!(is.na(test.age$coefficients[2,4]))) || (!(is.na(test.depth$coefficients[2,4]))))
      {
        ## To remove the trees for which contrast value is not properly standized, 
        ## we need to collect the P values of the diagnostic plot
        p.SDT<-test.expvar$coefficients[2,4]
        p.AgeT<-test.age$coefficients[2,4]
        p.depthT<-test.depth$coefficients[2,4]
       # p.heightT<-test.height$coefficients[2,4]
        
        ##If contrasts are not properly standized, we return NA
        if((p.SDT < 0.05) || (p.AgeT < 0.05) || (p.depthT < 0.05)){return (NA) }
        
        ## Collecting the initial pic data in a different variable name
        if((p.SDT >= 0.05) && (p.AgeT >= 0.05) && (p.depthT >= 0.05))
        {
          tree@data$pic_Tau<-NULL
          tree@data$variance<-NULL
          tree@data$pic_Exp<-NULL
          
          ## Now we return the trees dataframe
          tree@data$pic_Tau <- c(rep(NA, length(gene.tree$tip.label)), (contrast$Tau))
          tree@data$variance<- c(rep(NA, length(gene.tree$tip.label)), contrast$contrVar)
          tree@data$pic_Exp<-  c(rep(NA, length(gene.tree$tip.label)), (contrast$Mean.Exp))
          return(tree)
        }}}}
}

diagnostic.plot.test.exp<- function(tree,nodelimit)
{
  gene.tree<-tree@phylo
  gene.tree$node.label<-NULL
  count<<-count+1
  print(count)
  
  ## Collecting data
  data<-tree@data
  data.new<-data[which(!is.na(data$Tau)),] 
  data.new1<-data.frame(label=data.new$label,Tau=data.new$Tau, Mean.Exp=data.new$avg.exp)
  rownames(data.new1)<-data.new1$label
  treedata<-comparative.data(gene.tree,data.new1,label)
  test<-tryCatch(crunch(Mean.Exp~Tau,data=treedata,equal.branch.length=F), error = function(e) {"Error"})
  if((test=="Error") || (is.na(test$mod[1]$coefficients[1])))
  { print ("Error obtained in crunch")
    return(NA)
  }
  if(test!="Error")
  {
    ## We need to check for diagnostic plots for expected variance and node age from this Caper package
    contrast<-caic.table(test)
    contrast$node<-as.numeric(rownames(contrast))
    if(nodelimit==0){contrast$height<-data$nodeheight[which(data$node.age>0)]}
    if(nodelimit >0){contrast$height<-data$nodeheight[which(!is.na(data$events))]}
    
    
    ## Considering contrasts data without any age limit
    if(nodelimit==0){contrast.new<-contrast}
    if(nodelimit>0){contrast.new<-contrast[contrast$nodeAge<=nodelimit,]}
    
    ## length of the dataframe should be higher than 3
    if(nrow(contrast.new) >= 3)
    {
      ## Diagnstic tests for expected variance, node age, node depth and node height as recommended in earlier studies
      test.expvar<-summary(lm((abs(contrast.new$Mean.Exp))~sqrt(contrast.new$contrVar)))
      test.age<-summary(lm((abs(contrast.new$Mean.Exp))~log(contrast.new$nodeAge)))
      test.depth<-summary(lm((abs(contrast.new$Mean.Exp))~contrast.new$nodeDepth))
      
      if((!(is.na(test.expvar$coefficients[2,4]))) || (!(is.na(test.age$coefficients[2,4]))) || (!(is.na(test.depth$coefficients[2,4]))))
      {
        ## To remove the trees for which contrast value is not properly standized, 
        ## we need to collect the P values of the diagnostic plot
        p.SDE<-test.expvar$coefficients[2,4]
        p.AgeE<-test.age$coefficients[2,4]
        p.depthE<-test.depth$coefficients[2,4]
        
        ##If contrasts are not properly standized, we return NA
        if((p.SDE < 0.05) || (p.AgeE < 0.05) || (p.depthE < 0.05)){return (NA) }
        
        ## Collecting the initial pic data in a different variable name
        if((p.SDE >= 0.05) && (p.AgeE >= 0.05) && (p.depthE >= 0.05))
        {
          tree@data$pic_Exp<-NULL
          tree@data$variance<-NULL
          tree@data$pic_Tau<-NULL
          
          ## Now we return the trees dataframe
          tree@data$pic_Exp <- c(rep(NA, length(gene.tree$tip.label)), (contrast$Mean.Exp))
          tree@data$variance<- c(rep(NA, length(gene.tree$tip.label)), contrast$contrVar)
          tree@data$pic_Tau<- c(rep(NA, length(gene.tree$tip.label)), (contrast$Tau))
          return(tree)
        }}}}
}

diagnostic.plot.test.exp.all.tissue<- function(tree,tissue)
{
  gene.tree<-tree@phylo
  gene.tree$node.label<-NULL
  count<<-count+1
  print(count)
  
  ## Collecting data
  data<-tree@data
  data.new<-data[which(!is.na(data$TPM.brain)),] 
  tissue.new<-unlist(strsplit(tissue,"TPM."))[2]
  data.new1<-data.frame(label=data.new$label,Exp.Tissue=data.new[[tissue]],Mean.Exp=data.new$avg.exp)
  rownames(data.new1)<-data.new1$label
  treedata<-comparative.data(gene.tree,data.new1,label)
  test<-NULL
  test<-tryCatch(crunch(Exp.Tissue~Mean.Exp,data=treedata,equal.branch.length=F), error = function(e) {"Error"})
  if((test=="Error") || (is.na(test$mod[1]$coefficients[1])))
  { print ("Error obtained in crunch")
    return(NA)
  }
  if(test!="Error")
  {
    ## We need to check for diagnostic plots for expected variance and node age from the Caper package
    contrast<-caic.table(test)
    contrast$node<-as.numeric(rownames(contrast))
    contrast$height<-data$nodeheight[which(data$node.age>0)]
    
    ## Considering contrasts data 
    contrast.new<-contrast
    
    
    ## length of the dataframe should be higher than 3
    if(nrow(contrast.new) >= 3)
    {
      ## Diagnstic tests for expected variance, node age, node depth and node height as recommended in earlier studies
      test.expvar<-summary(lm((abs(contrast.new$Exp.Tissue))~sqrt(contrast.new$contrVar)))
      test.age<-summary(lm((abs(contrast.new$Exp.Tissue))~log(contrast.new$nodeAge)))
      test.depth<-summary(lm((abs(contrast.new$Exp.Tissue))~contrast.new$nodeDepth))
      
      
      if((!(is.na(test.expvar$coefficients[2,4]))) || (!(is.na(test.age$coefficients[2,4]))) || (!(is.na(test.depth$coefficients[2,4]))))
      {
        ## To remove the trees for which contrast value is not properly standized, 
        ## we need to collect the P values of the diagnostic plot
        p.SDE<-test.expvar$coefficients[2,4]
        p.AgeE<-test.age$coefficients[2,4]
        p.depthE<-test.depth$coefficients[2,4]
        
        ##If contrasts are not properly standized, we return NA
        if((p.SDE < 0.05) || (p.AgeE < 0.05) || (p.depthE < 0.05)){return (NA) }
        
        ## Collecting the initial pic data in a different variable name
        if((p.SDE >= 0.05) && (p.AgeE >= 0.05) && (p.depthE >= 0.05))
        {
          var<-paste0("pic.",tissue.new,sep="")
          tree@data[[var]]<-NULL
          tree@data$variance<-NULL
          tree@data$pic.new<-NULL
          
          ## Now we return the trees dataframe
          tree@data[[var]] <- c(rep(NA, length(gene.tree$tip.label)), (contrast.new$Exp.Tissue))
          tree@data$variance<- c(rep(NA, length(gene.tree$tip.label)), contrast.new$contrVar)
          tree@data$pic<- c(rep(NA, length(gene.tree$tip.label)), abs(contrast.new$Exp.Tissue))
          return(tree)
        }}}}
}

## Dropping tips without average exp data for species
drop.leaf.no.exp <- function(tree) 
{
  ##Reading the tips labels and data table
  tree_tip <- tree@phylo$tip.label
  gene_tree_data <- tree@data
  
  ##Identify the tips without expression data
  tips_to_remove <- gene_tree_data$label[which((is.na(gene_tree_data$avg.exp)) & (is.na(gene_tree_data$events)))]
  
  ##Identify the tips with expression data
  tips_with_exp <- gene_tree_data$label[which((!is.na(gene_tree_data$avg.exp)) & (is.na(gene_tree_data$events)))]
  
  
  ##Returning the pruned tree if atleast 4 tips have Tau data
  if(length(tips_with_exp) >= 4)
  {
    pruned.tree<-treeio::drop.tip(tree, tips_to_remove)
    return(pruned.tree)
  }
  else{return(NA)}
}


##This function is written to return summary function
summary.tree<-function(tree)
{
  summary.all.tree<-bind_rows(lapply(tree,
                                     function(tree){
                                       gene.data<-tree@data 
                                       gene.data$gene<-digest(tree) ## Creating hash index unique to each gene tree
                                       drops <- c("B","SIS")  # B and SIS have inconsistent types in binding rows, so we remove those columns
                                       tree@data<-gene.data[ , !(names(gene.data) %in% drops)]
                                       return(tree@data)}))
  return(summary.all.tree)
}

## This function is written to extract PIC values of speciation nodes
## frame is the dataframe and  x should be mentioned the column name (ex: 'Tau.abs' or 'Exp.abs')
speciation.contrast<-function(frame, x)
{
  contrast.spe <- frame[[x]][which((frame$Event=="speciation") | (frame$Event=="Speciation"))]
  contrast.spe <- abs(contrast.spe)
  return(contrast.spe)
}

## This function is written to extract PIC values of duplication nodes
## frame is the dataframe and  x should be mentioned the column name (ex: 'pic_Tau' or 'pic_Exp')
duplication.contrast<-function(frame, x)
{
  contrast.dup <- frame[[x]][which((frame$Event=="duplication") | (frame$Event=="Duplication"))]
  contrast.dup <- abs(contrast.dup)
  return(contrast.dup)
}

## This function is written to extract PIC values for duplication nodes of FishWGDs
## frame is the dataframe and  x should be mentioned the column name (ex: 'pic_Tau' or 'pic_Exp')
duplication.contrast.3R<-function(frame, x)
{
  contrast.dup <- frame[[x]][which((frame$Event=="FishWGD" & frame$label=="Clupeocephala") | (frame$Event=="FishWGD" & frame$label=="Osteoglossocephalai"))]
  contrast.dup <- abs(contrast.dup)
  return(contrast.dup)
}

## This function is written to generate nexux format of trees with Tau data for applying levolution
nexux.tree.tau<-function(tree)
{
  Frequency <<- Frequency + 1
  gene.tree<-tree@phylo
  gene.data<-tree@data
  data.mod<-gene.data[which(is.na(gene.data$events)),]
  data.mod.tau<-data.mod[c("label","Tau")]
  
  
  folder.new<-paste0(folder1,"To_pablo/Tau/", sep="")
  out.tree<-paste0(folder.new,"genetree.",Frequency,sep="")
  out.tau<-paste0(folder.new,"genetree.",Frequency,".Tau",sep="")
  
  ## Writing the file
  write.tree(gene.tree,out.tree)
  write.table(data.mod.tau,out.tau,row.names = F,col.names = F, sep="\t",quote = F)
  
}

## This function is written to generate nexux format of FishWGD trees with Tau data for applying levolution
nexux.tree.3Rtau<-function(tree)
{
  Frequency <<- Frequency + 1
  gene.tree<-tree@phylo
  gene.data<-tree@data
  data.mod<-gene.data[which(is.na(gene.data$events)),]
  data.mod.tau<-data.mod[c("label","Tau")]
  
  
  folder.new<-paste0(folder1,"To_pablo/Tau_3R/", sep="")
  out.tree<-paste0(folder.new,"genetree.",Frequency,sep="")
  out.tau<-paste0(folder.new,"genetree.",Frequency,".Tau",sep="")
  
  ## Writing the file
  write.tree(gene.tree,out.tree)
  write.table(data.mod.tau,out.tau,row.names = F,col.names = F, sep="\t",quote = F)
  
}

## This function is written to generate nexux format of FishWGD trees with average expression data for applying levolution
nexux.tree.3Rexp<-function(tree)
{
  Frequency <<- Frequency + 1
  gene.tree<-tree@phylo
  gene.data<-tree@data
  data.mod<-gene.data[which(is.na(gene.data$events)),]
  data.mod.tau<-data.mod[c("label","avg.exp")]
  
  
  folder.new<-paste0(folder1,"To_pablo/Exp_3R/", sep="")
  out.tree<-paste0(folder.new,"genetree.",Frequency,sep="")
  out.tau<-paste0(folder.new,"genetree.",Frequency,".avgExp",sep="")
  
  ## Writing the file
  write.tree(gene.tree,out.tree)
  write.table(data.mod.tau,out.tau,row.names = F,col.names = F, sep="\t",quote = F)
  
}


## This function is written to generate nexux format of trees with average expression data for applying levolution
nexux.tree.exp<-function(tree)
{
  Frequency <<- Frequency + 1
  gene.tree<-tree@phylo
  gene.data<-tree@data
  data.mod<-gene.data[which(is.na(gene.data$events)),]
  data.mod.exp<-data.mod[c("label","avg.exp")]
  
  
  folder.new<-paste0(folder1,"To_pablo/Avg_exp_new/", sep="")
  out.tree<-paste0(folder.new,"genetree.",Frequency,sep="")
  out.exp<-paste0(folder.new,"genetree.",Frequency,".Avg.exp",sep="")
  
  ## Writing the file
  write.tree(gene.tree,out.tree)
  write.table(data.mod.exp,out.exp,row.names = F,col.names = F, sep="\t",quote = F)
  
}


## This function is written to generate nexux format of trees and different tissue expression data for applying levolution
nexux.tree.tissue<-function(tree,tissue)
{
  Frequency <<- Frequency + 1
  gene.tree<-tree@phylo
  gene.data<-tree@data
  data.mod<-gene.data[which(is.na(gene.data$events)),]
  tissuen<-tolower(tissue)
  colname.tissue<-paste0("TPM.",tissuen,sep="")
  tissuedata<-data.mod[c("label",colname.tissue)]
  
  
  folder.new<-paste0(folder1,"To_pablo/",tissue,"/", sep="")
  out.tree<-paste0(folder.new,"genetree.",Frequency,sep="")
  out.tissue.exp<-paste0(folder.new,"genetree.",Frequency,".",tissue,sep="")
  
  ## Writing the file
  write.tree(gene.tree,out.tree)
  write.table(tissuedata,out.tissue.exp,row.names = F,col.names = F, sep="\t",quote = F)
  
}


## Function to compute Wilcoxon test of data1 and data2 using two sided test
two.tailed.wilcox<-function (data1, data2)
{
  wilcox.oc.two.tailed <- wilcox.test(data1,data2,alternative="two.sided")$p.value
  
  #star <-gtools::stars.pval(wilcox.oc.two.tailed)
  #star0<- stars.pval(0)
  if (wilcox.oc.two.tailed == 0)
  {
    label.p <- paste0("P < 2.2%*%10^{-16}")
    return(label.p)
  }
  if ((wilcox.oc.two.tailed != 0))
  {
    if((wilcox.oc.two.tailed <= 0.001) & (wilcox.oc.two.tailed > 2.2e-16))
    {
      wilcox.oc.two.tailed1 <- format(wilcox.oc.two.tailed, digits= 3, scientific = TRUE)
      v1<-as.numeric(unlist(strsplit(toString(wilcox.oc.two.tailed1), split='e', fixed=TRUE))[1])
      v2<-as.numeric(unlist(strsplit(toString(wilcox.oc.two.tailed1), split='e', fixed=TRUE))[2])
      wilcox.oc.two.tailed1 <-paste0(v1,"%*%10^{",v2,"}")
      #wilcox.oc.two.tailed1 <- sub("e"," × 10^{", wilcox.oc.two.tailed1,"}")
      label.p <-paste0("P == ",wilcox.oc.two.tailed1)
      return(label.p)
    }
    if((wilcox.oc.two.tailed < 0.05) & (wilcox.oc.two.tailed > 0.001))
    {
      wilcox.oc.two.tailed1 <- format(wilcox.oc.two.tailed, digits= 3, scientific = FALSE)
      label.p <-paste0("P == ",wilcox.oc.two.tailed1)
      return(label.p)
    }
    if(wilcox.oc.two.tailed < 2.2e-16)
    {
      label.p <- paste0("P < 2.2%*%10^{-16}")
      return(label.p)
    }
    if(wilcox.oc.two.tailed >= 0.05)
    {
      wilcox.oc.two.tailed1 <- format(wilcox.oc.two.tailed, digits= 3, scientific = FALSE)
      label.p <-paste0("P == ",wilcox.oc.two.tailed1)
      return(label.p)
    }
  }
} 

## Function to compute paired Wilcoxon two-tailed test 
paired.wilcox<-function (data1, data2)
{
  # wilcox_paired <- wilcox.test(data1,data2,paired = TRUE,alternative = "less")$p.value
  wilcox.paired <- as.numeric(wilcox.test(data1,data2,paired = TRUE)$p.value)
  
  #star <-gtools::stars.pval(wilcox.paired)
  #star0<- stars.pval(0)
  
  if ((wilcox.paired == 0) | (wilcox.paired < 2.2e-16))
  {
    label.p<-NULL
    label.p <- paste0("P < 2.2%*%10^{-16}")
    return(label.p)
  }
  if((wilcox.paired <= 0.001) & (wilcox.paired > 2.2e-16))
    { 
      label.p<-NULL
      wilcox.paired1 <- format(wilcox.paired, digits= 3, scientific = TRUE)
      v1<-as.numeric(unlist(strsplit(toString(wilcox.paired1), split='e', fixed=TRUE))[1])
      v2<-as.numeric(unlist(strsplit(toString(wilcox.paired1), split='e', fixed=TRUE))[2])
      wilcox.paired1 <-paste0(v1,"%*%10^{",v2,"}")
      #wilcox.paired1 <- sub("e"," × 10", wilcox.paired1)
      label.p <-paste0("P == ",wilcox.paired1)
      return(label.p)
  }
  if ((wilcox.paired < 0.05) & (wilcox.paired > 0.001))
  {
    label.p<-NULL
    wilcox.paired1 <- format(wilcox.paired, digits= 3, scientific = FALSE)
    label.p <-paste0("P == ",wilcox.paired1)
    return(label.p)
  }
 
  if(wilcox.paired >= 0.05)
  {
      label.p<-NULL
      wilcox.paired1 <- format(wilcox.paired, digits= 3, scientific = FALSE)
      label.p =paste0("P == ",wilcox.paired1)
      return(label.p)
  }
}

## Function to compute two-proportions z-test
proportion.test<-function(x1,x2,n1,n2)
{
  pvalue.prop<-prop.test(x = c(x1, x2), n = c(n1, n2))$p.value
  if ((pvalue.prop == 0) | (pvalue.prop < 2.2e-16))
  {
    label.p<-NULL
    label.p <- paste0("P < 2.2%*%10^{-16}")
  }
 if((pvalue.prop <= 0.001) & (pvalue.prop > 2.2e-16))
    { 
      label.p<-NULL
      pvalue.prop1 <- format(pvalue.prop, digits= 3, scientific = TRUE)
      v1<-as.numeric(unlist(strsplit(toString(pvalue.prop1), split='e', fixed=TRUE))[1])
      v2<-as.numeric(unlist(strsplit(toString(pvalue.prop1), split='e', fixed=TRUE))[2])
      pvalue.prop11 <-paste0(v1,"%*%10^{",v2,"}")
      #pvalue.prop1 <- sub("e"," × 10", pvalue.prop1)
      label.p <-paste0("P == ",pvalue.prop1)
      return(label.p)
  }
  if ((pvalue.prop < 0.05) & (pvalue.prop > 0.001))
  { 
    label.p<-NULL
    pvalue.prop1 <- format(pvalue.prop, digits= 3, scientific = FALSE)
    label.p <-paste0("P == ",pvalue.prop1)
    return(label.p)
  }
 
  if(pvalue.prop >= 0.05)
  {
    pvalue.prop1 <- format(pvalue.prop, digits= 3, scientific = FALSE)
    #pvalue.prop1 <- sub("e", " × 10", pvalue.prop1)
    label.p =paste0("P == ",pvalue.prop1)
  }
  return(label.p)
}


## This function allows permutation of Tau data of the tips of a tree 
shuffling.tau<-function(tree)
{
  ## Reading @phylo and @data slots of a tree
  gene.tree<-tree@phylo
  gene.data<-tree@data
  
  if (class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  
  ## Collecting Tau data of the tip
  Tau.tree<- gene.data$Tau[which(is.na(gene.data$events))]
  
  ## Permuting the trait data and adding the shuffled Tau data before returning tree as a treeio::treedata object
  Tau.new<-sample(x=Tau.tree,size = length(Tau.tree), replace = FALSE)
  tree@data$Tau<- c(Tau.new,rep(NA, times=gene.tree$Nnode))
  return(tree)
}

## This function allows permutation of Exp data of the tips of a tree 
shuffling.exp<-function(tree)
{
  ## Reading @phylo and @data slots of a tree
  gene.tree<-tree@phylo
  gene.data<-tree@data
  
  if (class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  
  ## Collecting Exp data of the tip
  Exp.tree<- gene.data$avg.exp[which(is.na(gene.data$events))]
  
  ## Permuting the trait data and adding the shuffled Tau data before returning tree as a treeio::treedata object
  Exp.new<-sample(x=Exp.tree,size = length(Exp.tree), replace = FALSE)
  tree@data$Exp<- c(Exp.new,rep(NA, times=gene.tree$Nnode))
  return(tree)
}

##This function calculates PIC of nodes for each tree based on trait data of tips
contrast.calc<-function(tree,trait)
{
  Count<<-Count+1
  print(Count)
  ## Reading tree data slot
  gene.tree<-tree@phylo
  gene.data<-tree@data
  
  ## Collecting the Trait data of the tips of a tree
  Trait.tip<- gene.data[[trait]][which(is.na(gene.data$events))]
  
  ## Initializing variable
  abs.variable<-NULL
  abs.variable<-paste0("pic_",trait,sep="")
  tree@data[[abs.variable]] <- NULL
  tree@data$variance <- NULL
  
  ## Calculating the phylogenetic independent contrasts of tree
  ## Returning the results to "data" frame of the tree 
  ## To resolve polytomic tree we used "multi2di"
  pic.tree <- ape::pic(Trait.tip, multi2di(gene.tree), scaled = T, var.contrasts=TRUE)

  ## However this using "multi2di", adds PIC and variance for nodes(afer resolving) for which data are absent in @data slot
  ## New we only took pic and variance data for nodes for which duplication and speciation data available
  contra<-pic.tree[,1][which(row.names(pic.tree) %in% gene.data$label)]
  var<-pic.tree[,2][which(row.names(pic.tree) %in% gene.data$label)]
  tree@data[[abs.variable]] <- abs(c(rep(NA, length(gene.tree$tip.label)), contra)) ## absolute PIC values
  tree@data$pic <- abs(c(rep(NA, length(gene.tree$tip.label)), contra))
  tree@data$variance <- c(rep( NA, length(gene.tree$tip.label)),var)
  return(tree)
}

### This function is written to identify trees which have at leaast 3 tips with fish data and atleast 1 outgroup tetrapod species data
filter1<-function(tree)
{ 
  ##Collecting tip labels (Ensembl.ID) of species
  treep <- tree@phylo$tip.label
  
  ##Tips with species data to be retained
  fish.spec <- c("ENSDARG","ENSELUG","ENSAMXG","ENSORLG00000","ENSONIG")
  outgroup.spec <- c("ENSG000","ENSMUSG","ENSRNOG","ENSBTAG","ENSMMUG","ENSGALG","ENSCAFG","ENSOCUG","ENSMPUG","ENSMODG")
  search.string <- paste0(fish.spec, collapse = "|")
  search.outgroup.string <- paste0(outgroup.spec, collapse = "|")
  fish.species.data.tip <- treep[grepl(search.string,treep)] 
  outgroup.species.data.tip <- treep[grepl(search.outgroup.string,treep)] 
  
  ## We considered trees with atleast 3 tips of fish data
  if((length(fish.species.data.tip) >= 3) & (length(outgroup.species.data.tip) >= 1))
  {
    return(tree)
  }
  else{return(NA)}
}

## This function is written to process the output of trait jump model
process.jump.output<-function(dataframe)
{
  Frequency <<- Frequency + 1
  
  ## Initializing variables
  index.node<-NULL
  gene.tree<-NULL
  gene.data<-NULL
  node.event<-NULL
  root.node.event<-NULL
  clade<-NULL
  pic<-NULL
  jump<-NULL
  
  ## Extracting data
  index.tree<-as.numeric(dataframe$index.6923tree.tau[Frequency])
  index.node<-as.numeric(dataframe$From[Frequency])
  jump<-as.numeric(dataframe$JumpProbability[Frequency])
  gene.tree<-trees.all.interest[[index.tree]]@phylo
  #gene.tree<-Common.WGD.trees2[[index.tree]]@phylo
  root.node<-gene.tree$Nnode+2
  gene.data<-trees.all.interest[[index.tree]]@data
  #gene.data<-Common.WGD.trees2[[index.tree]]@data
  node.event<-as.character(gene.data$events[which(gene.data$node==index.node)])
  root.node.event<-as.character(gene.data$events[which(gene.data$node==root.node)])
  clade<-as.character(gene.data$label[which(gene.data$node==index.node)])
  pic<-as.numeric(gene.data$pic[which(gene.data$node==index.node)])
 
  ## Returning result
  return(tibble(index=index.tree,
                node.event=node.event,
                clade.name=clade,
                root.evet=root.node.event,
                pic=abs(pic),
                jump.probability=jump))
}

## This function is written to process the output of trait jump model
process.jump.output2<-function(DataFrame)
{
  Frequency <<- Frequency + 1
  
  ## Initializing variables
  index.node<-NULL
  gene.tree<-NULL
  gene.data<-NULL
  node.event<-NULL
  root.node.event<-NULL
  clade<-NULL
  pic<-NULL
  jump<-NULL
  
  ## Extracting data
  index.tree<-as.numeric(DataFrame$Gene[Frequency])
  index.node<-as.numeric(DataFrame$From[Frequency])
  jump<-as.numeric(DataFrame$JumpProbability[Frequency])
  gene.tree<-trees.all.interest[[index.tree]]@phylo
  root.node<-gene.tree$Nnode+2
  gene.data<-trees.all.interest[[index.tree]]@data
  node.event<-as.character(gene.data$events[which(gene.data$node==index.node)])
  root.node.event<-as.character(gene.data$events[which(gene.data$node==root.node)])
  clade<-as.character(gene.data$label[which(gene.data$node==index.node)])
  pic<-as.numeric(gene.data$pic[which(gene.data$node==index.node)])
  
  ## Returning result
  return(tibble(index=index.tree,
                node.event=node.event,
                clade.name=clade,
                root.evet=root.node.event,
                pic=abs(pic),
                jump.probability=jump))
}

## This function is written to obtain tree index info
tree.index.collect<-function(tree)
{
  count<<-count+1
  
  ##Reading gene tree data slot
  gene.data<-tree@data
  tree.index<-gene.data$index.tree
 
  ## Returning the tibble
  return(tibble(tree.num=count,
                Gene=unique(tree.index)))
}

## This function is written to calculate the proportion of duplication, speciation and NA events for each tree
tree.data.collection<-function(tree)
{
  ##Reading gene tree data slot
  gene.tree<-tree@phylo
  gene.data<-tree@data
  ntips <-length(gene.tree$tip.label) ## Number of tips
  root<-as.numeric(ntips+1) 
  root.event<-gene.data$events[root]
  root.age<-gene.data$node.age[root]
  Internal.event.number<-as.numeric(length(gene.data$node[which(!is.na(gene.data$pic))]))
  index.old<-gene.tree$index.tree
  
  ## initializing variables
  pdup<-0
  pspe<-0
  pNA<-0
  
  
  ##Since, many internal node events are assigned as "NA", sum of proportion of speciation and duplication may not be equal to 1
  dup<-as.numeric(length(gene.data$node[which((!is.na(gene.data$pic)) & (gene.data$events %in% "duplication"))])) 
  pdup<-round((dup/Internal.event.number),2) ## Proportion of duplication event
  spe<-as.numeric(length(gene.data$node[which((!is.na(gene.data$pic)) & (gene.data$events %in% "speciation"))])) 
  pspe<-round((spe/Internal.event.number),2) ## Proportion of speciation events
  NA.event<-as.numeric(length(gene.data$node[which((!is.na(gene.data$pic)) & (is.na(gene.data$events)))])) 
  pNA<-round((NA.event/Internal.event.number),2) ## Proportion of NA event
  spe.var<-round(median(gene.data$variance[which((!is.na(gene.data$pic)) & (gene.data$events %in% "speciation"))]),0) ## variance of speciation
  dup.var<-round(median(gene.data$variance[which((!is.na(gene.data$pic)) & (gene.data$events %in% "duplication"))]),0) ## variance of duplication
  if(pdup==0){dup.var<-NA} 
  
  ## Returning the tibble
  return(tibble(tip.num=ntips,
                internal.events=Internal.event.number,
                dup.num=dup,
                dup.prop=pdup,
                spe.prop=pspe,
                NA.prop=pNA,
                dup.var=dup.var,
                spe.var=spe.var,
                root.age=root.age,
                root.event=root.event,
                index.old=index.old))
  
}

## This function is written to calculate the proportion of duplication, speciation and NA events for FishWGD trees
tree.data.collection.3RWGD<-function(tree)
{
  ##Reading gene tree data slot
  gene.tree<-tree@phylo
  gene.data<-tree@data
  ntips <-length(gene.tree$tip.label) ## Number of tips
  root<-as.numeric(ntips+1) 
  root.event<-gene.data$events[root]
  root.age<-gene.data$node.age[root]
  Internal.event.number<-as.numeric(length(gene.data$node[which(!is.na(gene.data$pic))]))
  index.old<-gene.tree$index.tree
  
  ## initializing variables
  pdup<-0
  pWGD<-0
  pspe<-0
  pNA<-0
  
  
  ##Since, many internal node events are assigned as "NA", sum of proportion of speciation and duplication may not be equal to 1
  dup<-as.numeric(length(gene.data$node[which((!is.na(gene.data$pic)) & (gene.data$events %in% "duplication"))])) 
  pdup<-round((dup/Internal.event.number),2) ## Proportion of duplication event
  WGD<-as.numeric(length(gene.data$node[which((!is.na(gene.data$pic)) & (gene.data$events %in% "FishWGD"))])) 
  pWGD<-round((WGD/Internal.event.number),2) ## Proportion of whole genome duplication events
  spe<-as.numeric(length(gene.data$node[which((!is.na(gene.data$pic)) & (gene.data$events %in% "speciation"))])) 
  pspe<-round((spe/Internal.event.number),2) ## Proportion of speciation events
  NA.event<-as.numeric(length(gene.data$node[which((!is.na(gene.data$pic)) & (is.na(gene.data$events)))])) 
  pNA<-round((NA.event/Internal.event.number),2) ## Proportion of NA event
  spe.var<-round(median(gene.data$variance[which((!is.na(gene.data$pic)) & (gene.data$events %in% "speciation"))]),0) ## variance of speciation
  dup.var<-round(median(gene.data$variance[which((!is.na(gene.data$pic)) & (gene.data$events %in% "duplication"))]),0) ## variance of duplication
  WGD.var<-round(median(gene.data$variance[which((!is.na(gene.data$pic)) & (gene.data$events %in% "FishWGD"))]),0) ## variance of whole genome duplication
  if(pdup==0){dup.var<-NA}
  if(pWGD==0){WGD.var<-NA}
  
  ## Returning the tibble
  return(tibble(tip.num=ntips,
                internal.events=Internal.event.number,
                dup.num=dup,
                dup.prop=pdup,
                spe.prop=pspe,
                NA.prop=pNA,
                dup.var=dup.var,
                spe.var=spe.var,
                root.age=root.age,
                root.event=root.event,
                index.old=index.old,
                WGD.num=WGD,
                WGD.prop=pWGD))
  
}


## This function is written to calculate the proportion of duplication, speciation and NA events for teleosts specific nodes of FishWGD trees
tree.data.collection.3RWGD.fish<-function(tree)
{
  ##Reading gene tree data slot
  gene.tree<-tree@phylo
  gene.data<-tree@data
  ntips <-length(gene.tree$tip.label) ## Number of tips
  root<-as.numeric(ntips+1) 
  root.event<-gene.data$events[root]
  root.age<-gene.data$node.age[root]
  Internal.event.number<-as.numeric(length(gene.data$node[which(!is.na(gene.data$pic))]))
  index.old<-gene.tree$index.tree
  
  ## Initializing variables
  pdup.fish<-0
  pWGD.fish<-0
  pspe.fish<-0
  pNA.fish<-0
  Internal.eventnumber.fish<-0
  
  ## Fish clades we consider
  fish.clade<-c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio",
                "Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha",
                "Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae",
                "Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus",
                "Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae",
                "Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")
  
  Internal.eventnumber.fish<-as.numeric(length(gene.data$node[which((!is.na(gene.data$pic)) & (gene.data$label %in% fish.clade))]))
  
  if(Internal.eventnumber.fish >= 1)
  { 
    fish.dup<-as.numeric(length(gene.data$node[which((!is.na(gene.data$pic)) & (gene.data$events %in% "duplication") & (gene.data$label %in% fish.clade))])) 
    pdup.fish<-round((fish.dup/Internal.eventnumber.fish),2) ## Proportion of duplication event
    fish.WGD<-as.numeric(length(gene.data$node[which((!is.na(gene.data$pic)) & (gene.data$events %in% "FishWGD") & (gene.data$label %in% fish.clade))])) 
    pWGD.fish<-round((fish.WGD/Internal.eventnumber.fish),2) ## Proportion of whole genome duplication events
    fish.spe<-as.numeric(length(gene.data$node[which((!is.na(gene.data$pic)) & (gene.data$events %in% "speciation") & (gene.data$label %in% fish.clade))])) 
    pspe.fish<-round((fish.spe/Internal.eventnumber.fish),2) ## Proportion of speciation events
    NA.event.fish<-as.numeric(length(gene.data$node[which((!is.na(gene.data$pic)) & (is.na(gene.data$events)) & (gene.data$label %in% fish.clade))])) 
    pNA.fish<-round((NA.event.fish/Internal.eventnumber.fish),2) ## Proportion of NA event
    spe.var<-round(median(gene.data$variance[which((!is.na(gene.data$pic)) & (gene.data$events %in% "speciation") & (gene.data$label %in% fish.clade))]),0) ## variance of speciation
    dup.var<-round(median(gene.data$variance[which((!is.na(gene.data$pic)) & (gene.data$events %in% "duplication") & (gene.data$label %in% fish.clade))]),0) ## variance of duplication
    WGD.var<-round(median(gene.data$variance[which((!is.na(gene.data$pic)) & (gene.data$events %in% "FishWGD") & (gene.data$label %in% fish.clade))]),0) ## variance of whole genome duplication
    if(pdup.fish==0){dup.var<-NA}
    if(pWGD.fish==0){WGD.var<-NA}
    
    ## Returning the tibble
    return(tibble(tip.num=ntips,
                  internal.events.fish=Internal.eventnumber.fish,
                  dup.num=fish.dup,
                  dup.prop=pdup.fish,
                  spe.prop=pspe.fish,
                  NA.prop=pNA.fish,
                  dup.var=dup.var,
                  spe.var=spe.var,
                  root.age=root.age,
                  root.event=root.event,
                  index.old=index.old,
                  WGD.num=fish.WGD,
                  WGD.prop=pWGD.fish))
  }
  if(Internal.eventnumber.fish == 0) {return (tibble(tip.num=ntips,
                                                     internal.events.fish=Internal.eventnumber.fish,
                                                     dup.num=NA,
                                                     dup.prop=NA,
                                                     spe.prop=NA,
                                                     NA.prop=NA,
                                                     dup.var=NA,
                                                     spe.var=NA,
                                                     root.age=root.age,
                                                     root.event=root.event,
                                                     index.old=index.old,
                                                     WGD.num=NA,
                                                     WGD.prop=NA))
                                      }
}



## This function is written to calculate the proportion of duplication, speciation and NA events for fish clades for each tree 
tree.data.collection.fish<-function(tree)
{
  count<<-count + 1
  print(count)
  ##Reading gene tree data slot
  gene.tree<-tree@phylo
  gene.data<-tree@data
  ntips <-length(gene.tree$tip.label) ## Number of tips
  root<-as.numeric(ntips+1) 
  root.event<-gene.data$events[root]
  root.age<-gene.data$node.age[root]
  index.old<-gene.tree$index.tree
  
  ## initializing variables
  pdup.fish<-0
  pspe.fish<-0
  pNA.fish<-0
  Internal.eventnumber.fish<-0
  
  ## Fish clades we consider
  fish.clade<-c("Astyanax.mexicanus","Pseudocrenilabrinae","Danio.rerio",
                "Acanthomorphata","Oryzias.latipes","Characiphysae","Euteleosteomorpha",
                "Osteoglossocephalai","Characoidei","Oryzias","Atherinomorphae",
                "Clupeocephala","Otophysi","Oryzias.latipes.ASM223467v1","Oreochromis.niloticus",
                "Cichlidae","Percomorphaceae_A","Percomorphaceae_B" ,"Percomorphaceae",
                "Esox.lucius","Ovalentaria","NAME_3","NAME_9","Neopterygii")
  
  Internal.eventnumber.fish<-as.numeric(length(gene.data$node[which((!is.na(gene.data$pic)) & (gene.data$label %in% fish.clade))]))

  
  if(Internal.eventnumber.fish >= 1)
  {  
    
    fish.dup<-as.numeric(length(gene.data$node[which((!is.na(gene.data$pic)) & (gene.data$events %in% "duplication") & (gene.data$label %in% fish.clade))])) 
    pdup.fish<-round((fish.dup/Internal.eventnumber.fish),2) ## Proportion of duplication event
    fish.spe<-as.numeric(length(gene.data$node[which((!is.na(gene.data$pic)) & (gene.data$events %in% "speciation") & (gene.data$label %in% fish.clade))])) 
    pspe.fish<-round((fish.spe/Internal.eventnumber.fish),2) ## Proportion of speciation events
    NA.event.fish<-as.numeric(length(gene.data$node[which((!is.na(gene.data$pic)) & (is.na(gene.data$events)) & (gene.data$label %in% fish.clade))])) 
    pNA.fish<-round((NA.event.fish/Internal.eventnumber.fish),2) ## Proportion of NA event
    spe.var<-round(median(gene.data$variance[which((!is.na(gene.data$pic)) & (gene.data$events %in% "speciation") & (gene.data$label %in% fish.clade))]),0) ## variance of speciation
    dup.var<-round(median(gene.data$variance[which((!is.na(gene.data$pic)) & (gene.data$events %in% "duplication") & (gene.data$label %in% fish.clade))]),0) ## variance of duplication
    if(pdup.fish==0){dup.var<-NA} 
  
    ## Returning the tibble
    return(tibble(tip.num=ntips,
                  internal.events.fish=Internal.eventnumber.fish,
                dup.num=fish.dup,
                dup.prop=pdup.fish,
                spe.prop=pspe.fish,
                NA.prop=pNA.fish,
                dup.var=dup.var,
                spe.var=spe.var,
                root.age=root.age,
                root.event=root.event,
                index.old=index.old))
  }
  if(Internal.eventnumber.fish == 0) {return (tibble(tip.num=ntips,
                                                     internal.events.fish=Internal.eventnumber.fish,
                                                      dup.num=NA,
                                                      dup.prop=NA,
                                                      spe.prop=NA,
                                                      NA.prop=NA,
                                                      dup.var=NA,
                                                      spe.var=NA,
                                                      root.age=root.age,
                                                      root.event=root.event,
                                                     index.old=index.old))
                                      }
}

## This function is written to obtain original index of trees after removal of null tree 
null.remove<- function(tree,indexing)
{ 

  count<<-count + 1
  print(count)
  if(((is.null(tree)) & (length(tree)==0)))
   {
    indexing<-count
    return(indexing)
   }
}

#### This function is written to paint tree (plotSimmap) to separate WGD nodes from other duplication nodes of the tree
paint.tree.mod<-function(tree)
{
  #count <<-count+1
  #print(count)
  ## Considering real tree data
  gene.tree<-tree@phylo
  gene.data<-tree@data
  ntips <- length(gene.tree$tip.label) 
  Internal.nodedata <- nrow(gene.data)-ntips
  Duplication.nodelength <- length(gene.data$events[which(gene.data$events == "duplication")])
  WGD.nodelength <- length(gene.data$events[which((gene.data$events == "FishWGD"))])
  Speciation.nodelength <- length(gene.data$events[which(gene.data$events == "speciation")])
  
  ## Now checking for trees with all duplication events
  ## If does not match the criteria returns NA
  if((Internal.nodedata - WGD.nodelength == 0) | (Internal.nodedata - Duplication.nodelength == 0) | (Internal.nodedata - Speciation.nodelength == 0))
  {
    return(NA)
  }
  else
  {
    
    ## Identifiying duplication nodes and edges to paint them to simulate them at specified trait evolutionary rate
    dup.nodes <- as.numeric(gene.data$node[which(gene.data$events=="duplication")])
    WGdup.nodes <- as.numeric(gene.data$node[which((gene.data$events == "FishWGD"))])
    WGdup.edges <- unique(gene.tree$edge[which(gene.tree$edge[,1] %in% WGdup.nodes), 2])
    tree.painted <- paintBranches (gene.tree, edge=WGdup.edges, "FishWGD", anc.state="Speciation")
    
    # Checking if we have other duplication nodes in a tree except WGD nodes
    other.dupnodes <- as.numeric(dup.nodes[!(dup.nodes %in% WGdup.nodes)])
    if(length(other.dupnodes) > 0) 
    {
      other.dupedges <- unique(gene.tree$edge[which(gene.tree$edge[,1] %in% other.dupnodes), 2])
      tree.painted<-paintBranches (tree.painted, edge=other.dupedges, "Duplication", anc.state="Speciation")
    }
    if("phylo" %in% class(tree.painted))
    {
      #tree@phylo$painted.tree<- tree.painted
      return(tree.painted)
    }
    else{return(NA)}
  }
}

## This function is written to identify synteny based WGD identification
## We obtained 3RWGD from ohnolog database by selecting only "fishWGD" duplication time
validate.WGD.2species<-function(tree,list)
{
  ## Considering real tree data
  gene.tree<-tree@phylo
  data<-tree@data
  
  ## reading data
  Count <<-Count+1
  print(Count)
  Gene<- paste0("OG",Count,sep="")
  tip<-data$label[which(is.na(data$events))]
  root.node<-length(tip)+1
  root.age<-data$node.age[which(data$node==root.node)]
  zebrafish.duplicates<-tip[grepl("ENSDARG",tip)]
  zfish.dup<-length(zebrafish.duplicates)
  medaka.duplicates<-tip[grepl("ENSORLG",tip)]
  medaka.dup<-length(medaka.duplicates)

  if(zfish.dup==0 & medaka.dup==0){return(NA)}
  if(zfish.dup>0){
    if((zfish.dup==2) & (zebrafish.duplicates[1] %in% c(list$Ohno2,list$Ohno1)) & (zebrafish.duplicates[2] %in% c(list$Ohno2,list$Ohno1)))
    {
      return(tree)
    }
    if((zfish.dup==1) & (zebrafish.duplicates[1] %in% list$Ohno1) | (zebrafish.duplicates[1] %in% list$Ohno2))
    {
      return(tree)
    }
    if((zfish.dup>2) & (zebrafish.duplicates %in% c(list$Ohno2,list$Ohno1)))
    {
      return(tree)
    }
  }
  if(medaka.dup>0){
    if((medaka.dup==2) & (medaka.duplicates[1] %in% c(list$Ohno2,list$Ohno1)) & (medaka.duplicates[2] %in% c(list$Ohno2,list$Ohno1)))
    {
      return(tree)
    }
    if((medaka.dup==1) & (medaka.duplicates[1] %in% list$Ohno1) | (medaka.duplicates[1] %in% list$Ohno2))
    {
      return(tree)
    }
    if((medaka.dup>2) & (medaka.duplicates %in% c(list$Ohno2,list$Ohno1)))
    {
      return(tree)
    }
  }
}


## This function is written to identify fishWGD node using FishWGD trees for zebrafish and medaka
mrca.WGD<-function(tree,list)
{
  count<<-count+1
  print(count)
  info<-tree@data
  
  ## Collecting tip labels
  tiplabel<-tree@phylo$tip.label
  zfish.id<-tiplabel[grepl("ENSDARG",tiplabel)]
  medaka.id<-tiplabel[grepl("ENSORLG",tiplabel)]
  
  ##Matching ohnolog pairs for zebrafish and medaka from Ohnolog database
  zpair1<-zfish.id[which(zfish.id %in% list$Ohno1)]
  zpair2<-list$Ohno2[which(list$Ohno1 %in% zpair1)]
  mpair1<-medaka.id[which(medaka.id %in% list$Ohno1)]
  mpair2<-list$Ohno2[which(list$Ohno1 %in% mpair1)]
  
  tip.combinations<-NULL
  tip.combinations1<-NULL
  if(((length(zpair1)==1) & (length(zpair2)==1)) | ((length(mpair1==1) & (length(mpair2)==1))))
  {
    if((length(zpair1)==1) & (length(zpair2)==1))
    {
      req.tiplabels<-NULL
      for(i in 1:length(tiplabel))
      {
        if(tiplabel[i] %in% zpair1){req.tiplabels<-append(req.tiplabels,i) }
        if(tiplabel[i] %in% zpair2){req.tiplabels<-append(req.tiplabels,i) }
      }
    
      if(length(req.tiplabels)>=2)
      {
        ## Pairwise combinations of tips
        tip.combinations<-as_tibble(t(combn(req.tiplabels,2)))
      }
    }
    if((length(mpair1)==1) & (length(mpair2)==1))
    {
      req.tiplabels1<-NULL
      for(i in 1:length(tiplabel))
      {
        if(tiplabel[i] %in% mpair1){req.tiplabels1<-append(req.tiplabels1,i) }
        if(tiplabel[i] %in% mpair2){req.tiplabels1<-append(req.tiplabels1,i) }
      }
    
      if(length(req.tiplabels1)>=2)
      {
        ## Pairwise combinations of tips
        tip.combinations1<-as_tibble(t(combn(req.tiplabels1,2)))
      }
    }  
    tip.combinations<-bind_rows(tip.combinations,tip.combinations1)
    if(nrow(tip.combinations)>=1)
    {
      names(tip.combinations) = c("tip1", "tip2")
      tip.combinations$name1<-info$label[tip.combinations$tip1]
      tip.combinations$name2<-info$label[tip.combinations$tip2]
    
      # For each row, get the most recent common ancestor of the pairwise tip combination
      tip.combinations$node <- as.integer(apply(tip.combinations, 1, 
                                                function(x) ape::getMRCA( 
                                                  tree@phylo, 
                                                  c(as.integer(x['tip1']), as.integer(x['tip2'])))))
      
      tip.combinations$mrca<-info$label[tip.combinations$node]
      tip.combinations$event<-info$events[tip.combinations$node]
      tree@data$events<-as.character(tree@data$events)
      FishWGD.nodes <- as.numeric(tip.combinations$node)
      tree@data$events[FishWGD.nodes]<-"FishWGD"
      tree@data$events<- factor(tree@data$events, levels=c("speciation", "duplication","FishWGD"))
      return(tree)
    }
    else{return(NA)}
  }
}


## This function is written to identify teleost specific FishWGD trees 
mrca.WGD.types<-function(tree,list,type)
{
  count<<-count+1
  print(count)
  info<-tree@data
  
  ## Collecting tip labels
  tiplabel<-tree@phylo$tip.label
  zfish.id<-tiplabel[grepl("ENSDARG",tiplabel)]
  medaka.id<-tiplabel[grepl("ENSORLG",tiplabel)]
  
  ##Matching ohnolog pairs for zebrafish and medaka from Ohnolog database
  zpair1<-zfish.id[which(zfish.id %in% list$Ohno1)]
  zpair2<-list$Ohno2[which(list$Ohno1 %in% zpair1)]
  mpair1<-medaka.id[which(medaka.id %in% list$Ohno1)]
  mpair2<-list$Ohno2[which(list$Ohno1 %in% mpair1)]
  
  tip.combinations<-NULL
  tip.combinations1<-NULL
  if(((length(zpair1)==1) & (length(zpair2)==1)) | ((length(mpair1==1) & (length(mpair2)==1))))
  {
    if((length(zpair1)==1) & (length(zpair2)==1))
    {
      req.tiplabels<-NULL
      for(i in 1:length(tiplabel))
      {
        if(tiplabel[i] %in% zpair1){req.tiplabels<-append(req.tiplabels,i) }
        if(tiplabel[i] %in% zpair2){req.tiplabels<-append(req.tiplabels,i) }
      }
      
      if(length(req.tiplabels)>=2)
      {
        ## Pairwise combinations of tips
        tip.combinations<-as_tibble(t(combn(req.tiplabels,2)))
      }
    }
    if((length(mpair1)==1) & (length(mpair2)==1))
    {
      req.tiplabels1<-NULL
      for(i in 1:length(tiplabel))
      {
        if(tiplabel[i] %in% mpair1){req.tiplabels1<-append(req.tiplabels1,i) }
        if(tiplabel[i] %in% mpair2){req.tiplabels1<-append(req.tiplabels1,i) }
      }
      
      if(length(req.tiplabels1)>=2)
      {
        ## Pairwise combinations of tips
        tip.combinations1<-as_tibble(t(combn(req.tiplabels1,2)))
      }
    }  
    tip.combinations<-bind_rows(tip.combinations,tip.combinations1)
    if(nrow(tip.combinations)>=1)
    {
      names(tip.combinations) = c("tip1", "tip2")
      tip.combinations$name1<-info$label[tip.combinations$tip1]
      tip.combinations$name2<-info$label[tip.combinations$tip2]
      
      # For each row, get the most recent common ancestor of the pairwise tip combination
      tip.combinations$node <- as.integer(apply(tip.combinations, 1, 
                                                function(x) ape::getMRCA( 
                                                  tree@phylo, 
                                                  c(as.integer(x['tip1']), as.integer(x['tip2'])))))
      
      tip.combinations$mrca<-info$label[tip.combinations$node]
      tip.combinations$event<-info$events[tip.combinations$node]
      if(type=="sure")
      {
        if((tip.combinations$event=="duplication" & tip.combinations$mrca=="Clupeocephala") | (tip.combinations$event=="duplication" & tip.combinations$mrca=="Osteoglossocephalai"))
        {
          tree@data$events<-as.character(tree@data$events)
          FishWGD.nodes <- as.numeric(tip.combinations$node)
          tree@data$events[FishWGD.nodes]<-"FishWGD"
          tree@data$events<- factor(tree@data$events, levels=c("speciation", "duplication","FishWGD"))
          return(tree)
        }
        else{return(NA)}
      }
      
      if(type=="unsure")
      {
        if(!(tip.combinations$mrca=="Clupeocephala" | tip.combinations$mrca=="Osteoglossocephalai"))
        {
          tree@data$events<-as.character(tree@data$events)
          FishWGD.nodes <- as.numeric(tip.combinations$node)
          tree@data$events[FishWGD.nodes]<-"FishWGD"
          tree@data$events<- factor(tree@data$events, levels=c("speciation", "duplication","FishWGD"))
          return(tree)
        }
        else{return(NA)}
      }
    }
  }
}

## Function to identify candidate SSD trees  
## Identifying trees with no "Cluepeocephaha" or "Osteoglossocephalai" duplication 
candidate.unsureSSD.phylogeny<-function(tree)
{
  #Flag initialization 
  flag<-0
  
  # Reading the '@data' slot & identifying internal node clade labels and events
  gene.tree<-tree@phylo
  gene.data<-tree@data
  tiplabels<-gene.tree$tip.label
  #internalnodes<-internal.nodes(tree)
  internal.labels <-gene.data$label[which(!is.na(gene.data$events))]
   
  if(("Clupeocephala_d" %in% internal.labels) | ("Osteoglossocephalai_d" %in% internal.labels)){flag=1} ##If there is duplication at "Clupeocephala" or "Osteoglossocephalai" clade; we do not consider the trees
  
  if(flag==0)
  {
     return(tree)
  } 
  if(flag==1){return(NA)}
}
