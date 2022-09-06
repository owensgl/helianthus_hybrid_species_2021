library(phytools)
library(tidyverse)
library(TreeTools)
library(castor)
library(textclean)

tree_names <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10.2kb.physical.subsample.treenames.trees.txt",
                       col_names = c("tree_id","tree"))

tree_names$tree <- mgsub(tree_names$tree, c("NHKG_16_01_R_196", "GO8_GB179","PAR_posas_01","ANN1283","ARG0295","DEB_1135","PET0662","PET0424","PET0495","664647_GIG"), 
                         c("ano", "des","par","ann","arg","deb","niv", "pet_fal","pet_pet","outgroup"))


pet_positions <- tibble()
for (i in 1:nrow(tree_names)){
  print(i)
  #Test for pet_pet
  tmp_tree <- ape::read.tree(text=tree_names$tree[i])
  tmp_tree <- ape::drop.tip(tmp_tree, "par")
  tmp_tree <- root(tmp_tree, outgroup = "outgroup", resolve.root = TRUE)
  tmp_tree <- ape::drop.tip(tmp_tree, "pet_fal")
  node_mrca <- findMRCA(tmp_tree, tips=c("ann","arg","pet_pet"),type="node")
  sub_tree <- extract.clade(tmp_tree, node_mrca)
  if(length(sub_tree$tip.label) == 3){
    #Test if ANN,ARG are sister
    node_mrca_1 <- findMRCA(sub_tree, tips=c("ann","arg"),type="node")
    sub_tree_1 <- extract.clade(sub_tree, node_mrca_1)
    if(length(sub_tree_1$tip.label) == 2){
      tmp <- tibble(tree = i, species = "pet_pet",pattern = "ann_arg")
      pet_positions <- rbind(pet_positions, tmp)
    }
    node_mrca_2 <- findMRCA(sub_tree, tips=c("ann","pet_pet"),type="node")
    sub_tree_2 <- extract.clade(sub_tree, node_mrca_2)
    if(length(sub_tree_2$tip.label) == 2){
      tmp <- tibble(tree = i, species = "pet_pet",pattern = "ann_pet")
      pet_positions <- rbind(pet_positions, tmp)
      
    }
    node_mrca_3 <- findMRCA(sub_tree, tips=c("arg","pet_pet"),type="node")
    sub_tree_3 <- extract.clade(sub_tree, node_mrca_3)
    if(length(sub_tree_3$tip.label) == 2){
      tmp <- tibble(tree = i, species = "pet_pet",pattern = "arg_pet")
      pet_positions <- rbind(pet_positions, tmp)
    }
    
  }
  
  #Test for pet_fal
  tmp_tree <- ape::read.tree(text=tree_names$tree[i])
  tmp_tree <- ape::drop.tip(tmp_tree, "par")
  tmp_tree <- root(tmp_tree, outgroup = "outgroup", resolve.root = TRUE)
  tmp_tree <- ape::drop.tip(tmp_tree, "pet_pet")
  node_mrca <- findMRCA(tmp_tree, tips=c("ann","arg","pet_fal"),type="node")
  sub_tree <- extract.clade(tmp_tree, node_mrca)
  if(length(sub_tree$tip.label) == 3){
    #Test if ANN,ARG are sister
    node_mrca_1 <- findMRCA(sub_tree, tips=c("ann","arg"),type="node")
    sub_tree_1 <- extract.clade(sub_tree, node_mrca_1)
    if(length(sub_tree_1$tip.label) == 2){
      tmp <- tibble(tree = i, species = "pet_fal",pattern = "ann_arg")
      pet_positions <- rbind(pet_positions, tmp)
    }
    node_mrca_2 <- findMRCA(sub_tree, tips=c("ann","pet_fal"),type="node")
    sub_tree_2 <- extract.clade(sub_tree, node_mrca_2)
    if(length(sub_tree_2$tip.label) == 2){
      tmp <- tibble(tree = i, species = "pet_fal",pattern = "ann_pet")
      pet_positions <- rbind(pet_positions, tmp)
      
    }
    node_mrca_3 <- findMRCA(sub_tree, tips=c("arg","pet_fal"),type="node")
    sub_tree_3 <- extract.clade(sub_tree, node_mrca_3)
    if(length(sub_tree_3$tip.label) == 2){
      tmp <- tibble(tree = i, species = "pet_fal",pattern = "arg_pet")
      pet_positions <- rbind(pet_positions, tmp)
    }
    
  }
}

plot(ape::read.tree(text=tree_names$tree[254]))
pet_positions %>%
  group_by(species, pattern) %>%
  summarize(count=n()) %>%
  ggplot(.,aes(x=species,y=count,fill=pattern)) +
  geom_bar(stat="identity",position="dodge")
 