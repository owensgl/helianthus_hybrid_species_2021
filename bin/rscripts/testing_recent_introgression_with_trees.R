library(phytools)
library(tidyverse)
library(TreeTools)
library(castor)
library(textclean)

#Test to see for trees where ann, arg and pet are monophyletic, whether its ann,arg or ann,pet
tree_names <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10.2kb.physical.treenames.trees.txt",
                       col_names = c("tree_id","tree"))

pet_pet <- c("PET0568","PET0495")
pet_fal <- c("PET0765","PET0424")
ann <- c("ANN1029","ANN1283")
arg <- c("ARG0143","ARG0295")
par <- c("PAR_3","PAR_posas_01")

removal_list <- c(pet_pet,pet_fal, ann,arg,par)
#Test for pet_pet
pet_positions <- tibble()
for (i in 1:nrow(tree_names)){
  print(i)
  tmp_tree <- ape::read.tree(text=tree_names$tree[i])
  #select the samples to use
  for (pet_sample in pet_pet){
    for (ann_sample in ann){
      for (arg_sample in arg){
        trio_name <- paste(ann_sample,arg_sample,pet_sample,sep=".")
        tmp_removal <- removal_list[which(removal_list != ann_sample & 
                                               removal_list != arg_sample &
                                            removal_list != pet_sample)]
        tmp_tree <- ape::read.tree(text=tree_names$tree[i])
        tmp_tree <- ape::drop.tip(tmp_tree, tmp_removal)
        tmp_tree <- root(tmp_tree, outgroup = "664647_GIG", resolve.root = TRUE)
        trio <- c(ann_sample,arg_sample,pet_sample)
        node_mrca <- findMRCA(tmp_tree, tips=trio,type="node")
        sub_tree <- extract.clade(tmp_tree, node_mrca)
        if(length(sub_tree$tip.label) == 3){
          #Test if ANN,ARG are sister
          node_mrca_1 <- findMRCA(sub_tree, tips=c(ann_sample,arg_sample),type="node")
          sub_tree_1 <- extract.clade(sub_tree, node_mrca_1)
          if(length(sub_tree_1$tip.label) == 2){
            tmp <- tibble(tree = i, species = "pet_pet",pattern = "ann_arg",trio=trio_name)
            pet_positions <- rbind(pet_positions, tmp)
          }
          node_mrca_2 <- findMRCA(sub_tree, tips=c(ann_sample,pet_sample),type="node")
          sub_tree_2 <- extract.clade(sub_tree, node_mrca_2)
          if(length(sub_tree_2$tip.label) == 2){
            tmp <- tibble(tree = i, species = "pet_pet",pattern = "ann_pet",trio=trio_name)
            pet_positions <- rbind(pet_positions, tmp)
            
          }
          node_mrca_3 <- findMRCA(sub_tree, tips=c(arg_sample,pet_sample),type="node")
          sub_tree_3 <- extract.clade(sub_tree, node_mrca_3)
          if(length(sub_tree_3$tip.label) == 2){
            tmp <- tibble(tree = i, species = "pet_pet",pattern = "arg_pet",trio=trio_name)
            pet_positions <- rbind(pet_positions, tmp)
          }
        }
      }
    }
  }
}
#Testing for pet_fal
for (i in 1:nrow(tree_names)){
  print(i)
  tmp_tree <- ape::read.tree(text=tree_names$tree[i])
  #select the samples to use
  for (pet_sample in pet_fal){
    for (ann_sample in ann){
      for (arg_sample in arg){
        trio_name <- paste(ann_sample,arg_sample,pet_sample,sep=".")
        tmp_removal <- removal_list[which(removal_list != ann_sample & 
                                            removal_list != arg_sample &
                                            removal_list != pet_sample)]
        tmp_tree <- ape::read.tree(text=tree_names$tree[i])
        tmp_tree <- ape::drop.tip(tmp_tree, tmp_removal)
        tmp_tree <- root(tmp_tree, outgroup = "664647_GIG", resolve.root = TRUE)
        trio <- c(ann_sample,arg_sample,pet_sample)
        node_mrca <- findMRCA(tmp_tree, tips=trio,type="node")
        sub_tree <- extract.clade(tmp_tree, node_mrca)
        if(length(sub_tree$tip.label) == 3){
          #Test if ANN,ARG are sister
          node_mrca_1 <- findMRCA(sub_tree, tips=c(ann_sample,arg_sample),type="node")
          sub_tree_1 <- extract.clade(sub_tree, node_mrca_1)
          if(length(sub_tree_1$tip.label) == 2){
            tmp <- tibble(tree = i, species = "pet_fal",pattern = "ann_arg",trio=trio_name)
            pet_positions <- rbind(pet_positions, tmp)
          }
          node_mrca_2 <- findMRCA(sub_tree, tips=c(ann_sample,pet_sample),type="node")
          sub_tree_2 <- extract.clade(sub_tree, node_mrca_2)
          if(length(sub_tree_2$tip.label) == 2){
            tmp <- tibble(tree = i, species = "pet_fal",pattern = "ann_pet",trio=trio_name)
            pet_positions <- rbind(pet_positions, tmp)
            
          }
          node_mrca_3 <- findMRCA(sub_tree, tips=c(arg_sample,pet_sample),type="node")
          sub_tree_3 <- extract.clade(sub_tree, node_mrca_3)
          if(length(sub_tree_3$tip.label) == 2){
            tmp <- tibble(tree = i, species = "pet_fal",pattern = "arg_pet",trio=trio_name)
            pet_positions <- rbind(pet_positions, tmp)
          }
        }
      }
    }
  }
}

pet_positions %>%
  group_by(species, pattern) %>%
  summarize(count=n()) %>%
  ggplot(.,aes(x=species,y=count,fill=pattern)) +
  geom_bar(stat="identity",position="dodge")  
