# Libraries ----
library(ape)
library(dplyr)
library(phytools)
library(tidytree)

# Functions ----
# Map Selectome node (branch_number) to node in processed gene tree
retrieve_mrca_node <- function(tree, descendants) {
  # Get subset of child_genes_in_branch that are present in the tree
  tips <- subset(descendants, descendants %in% tree$tip.label)
  
  # Find MRCA of tips. If not found, set node to NA.
  mrca.node <- phytools::findMRCA(tree, tips)
  if ( is.null(mrca.node) ) { mrca.node <- NA }
  
  return(mrca.node)
}

# Main ----
# setwd("/Users/cbucao/Google_Drive/UNIL/Projects/Tina_Positive_Selection/Scripts")

# Import data ----
# Load gene trees
path <- "../Gene_trees"
files <- list.files(path)

trees <- lapply(files, function(tree.file) { read.tree(file.path(path, tree.file)) })
names(trees) <- files
# 6923 gene trees

# Useful for checking
trees.tbl <- lapply(trees, as_tibble)

# Load Selectome files
selectome <- read.table(file.path("..", "selectome_trees_output.tsv"), header=TRUE, sep="\t")
# 104185 branches across 6649 trees
# 2603 branches with selection (q < 0.05)

# Process Selectome results ----
selectome$child_genes_of_branch <- strsplit(selectome$child_genes_of_branch, ",")

mapped.selectome.results <- 
  selectome[,c("tree_file", "selectome_tree_ac", "branch_number",
               "child_genes_of_branch", "taxon", "taxon_id",
               "lrt", "pvalue", "qvalue", "selected")]

# Map Selectome node (branch_number) to corresponding node in processed gene tree
mapped.selectome.results$node <- 
  apply(mapped.selectome.results, MARGIN=1, function(row) {
    retrieve_mrca_node(trees[[row$tree_file]], row$child_genes_of_branch)
    })

# Remove nodes with no corresponding tips in the tree
mapped.selectome.results <-
  mapped.selectome.results[!is.na(mapped.selectome.results$node),]
# 64643 branches across 6448 trees
# 2113 branches with selection (q < 0.05)

# Retrieve corresponding parent nodes
mapped.selectome.results$parent <-
  apply(mapped.selectome.results, MARGIN=1, function(row) {
    p <- trees[[row$tree_file]]$edge[,1][trees[[row$tree_file]]$edge[,2]==row$node]
    if (length(p) == 0) { p <- row$node } # No parent node
    parent <- p
  }) %>%
  unlist()

# Retrieve taxon labels
mapped.selectome.results$label <-
  apply(mapped.selectome.results, MARGIN=1, function(row) {
    ntips <- length(trees[[row$tree_file]]$tip.label)
    label <- trees[[row$tree_file]]$node.label[row$node-ntips]
  })

# Important: Check if taxon labels of Selectome results and processed gene trees are the same
summary(mapped.selectome.results$taxon == mapped.selectome.results$label)
# 50632 true
# 14011 false

# Keep only rows with matching taxon labels
keep <- mapped.selectome.results$taxon == mapped.selectome.results$label
mapped.selectome.results.true <- mapped.selectome.results[keep,] 
# 50632 branches across 6396 trees
# 1514 branches with selection (q < 0.05)

# mapped.selectome.results.false <- mapped.selectome.results[!keep,]

final.mapped.selectome.results.true <- 
  mapped.selectome.results.true[,c("tree_file", "parent", "node", "label",
                                   "lrt", "pvalue", "qvalue", "selected")]

# Save image and export file ----
write.table(final.mapped.selectome.results.true, file.path("..", "mapped.selectome.results.tsv"), row.names=FALSE, sep="\t")

save.image("mapped_selectome_to_trees.Rdata")
