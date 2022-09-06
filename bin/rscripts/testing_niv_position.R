library(phytools)
library(tidyverse)
library(TreeTools)
library(castor)
library(textclean)

#Test to see for trees where niv, petpet and petfal are monophyletic and where they fall
tree_names <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10.2kb.physical.treenames.trees.txt",
                       col_names = c("tree_id","tree"))

pet_pet <- c("PET0568","PET0495")
pet_fal <- c("PET0765","PET0424")
niv <- c("PET0695","PET0662")
deb <- c("DEB_1837","DEB_1135")
ann <- c("ANN1029","ANN1283")
arg <- c("ARG0143","ARG0295")
par <- c("PAR_3","PAR_posas_01")

removal_list <- c(pet_pet,pet_fal, niv)
#Test for pet_pet
niv_positions <- tibble()
for (i in 1:nrow(tree_names)){
  print(i)
  tmp_tree <- ape::read.tree(text=tree_names$tree[i])
  #select the samples to use
  for (petpet_sample in pet_pet){
    for (petfal_sample in pet_fal){
      for (niv_sample in niv){
        trio_name <- paste(petpet_sample,petfal_sample,niv_sample,sep=".")
        tmp_removal <- removal_list[which(removal_list != petpet_sample & 
                                            removal_list != petfal_sample &
                                            removal_list != niv_sample)]
        tmp_tree <- ape::read.tree(text=tree_names$tree[i])
        tmp_tree <- ape::drop.tip(tmp_tree, tmp_removal)
        tmp_tree <- root(tmp_tree, outgroup = "664647_GIG", resolve.root = TRUE)
        trio <- c(petpet_sample,petfal_sample,niv_sample)
        node_mrca <- findMRCA(tmp_tree, tips=trio,type="node")
        sub_tree <- extract.clade(tmp_tree, node_mrca)
        if(length(sub_tree$tip.label) == 3){
          #Test if ANN,ARG are sister
          node_mrca_1 <- findMRCA(sub_tree, tips=c(petpet_sample,petfal_sample),type="node")
          sub_tree_1 <- extract.clade(sub_tree, node_mrca_1)
          if(length(sub_tree_1$tip.label) == 2){
            tmp <- tibble(tree = i,pattern = "petpet_petfal",trio=trio_name)
            niv_positions <- rbind(niv_positions, tmp)
          }
          node_mrca_2 <- findMRCA(sub_tree, tips=c(niv_sample,petpet_sample),type="node")
          sub_tree_2 <- extract.clade(sub_tree, node_mrca_2)
          if(length(sub_tree_2$tip.label) == 2){
            tmp <- tibble(tree = i,pattern = "petpet_niv",trio=trio_name)
            niv_positions <- rbind(niv_positions, tmp)
            
          }
          node_mrca_3 <- findMRCA(sub_tree, tips=c(niv_sample,petfal_sample),type="node")
          sub_tree_3 <- extract.clade(sub_tree, node_mrca_3)
          if(length(sub_tree_3$tip.label) == 2){
            tmp <- tibble(tree = i,pattern = "petfal_niv",trio=trio_name)
            niv_positions <- rbind(niv_positions, tmp)
          }
        }
      }
    }
  }
}

niv_positions %>%
  group_by(trio, pattern) %>%
  summarize(count=n()) %>%
  ggplot(.,aes(x=trio,y=count,fill=pattern)) +
  geom_bar(stat="identity",position="dodge")  


removal_list <- c(ann,arg, deb)
#Test for deb
deb_positions <- tibble()
for (i in 1:nrow(tree_names)){
  print(i)
  tmp_tree <- ape::read.tree(text=tree_names$tree[i])
  #select the samples to use
  for (ann_sample in ann){
    for (arg_sample in arg){
      for (deb_sample in deb){
        trio_name <- paste(ann_sample,arg_sample,deb_sample,sep=".")
        tmp_removal <- removal_list[which(removal_list != ann_sample & 
                                            removal_list != arg_sample &
                                            removal_list != deb_sample)]
        tmp_tree <- ape::read.tree(text=tree_names$tree[i])
        tmp_tree <- ape::drop.tip(tmp_tree, tmp_removal)
        tmp_tree <- root(tmp_tree, outgroup = "664647_GIG", resolve.root = TRUE)
        trio <- c(ann_sample,arg_sample,deb_sample)
        node_mrca <- findMRCA(tmp_tree, tips=trio,type="node")
        sub_tree <- extract.clade(tmp_tree, node_mrca)
        if(length(sub_tree$tip.label) == 3){
          #Test if ANN,ARG are sister
          node_mrca_1 <- findMRCA(sub_tree, tips=c(ann_sample,arg_sample),type="node")
          sub_tree_1 <- extract.clade(sub_tree, node_mrca_1)
          if(length(sub_tree_1$tip.label) == 2){
            tmp <- tibble(tree = i,pattern = "ann_arg",trio=trio_name)
            deb_positions <- rbind(deb_positions, tmp)
          }
          node_mrca_2 <- findMRCA(sub_tree, tips=c(ann_sample,deb_sample),type="node")
          sub_tree_2 <- extract.clade(sub_tree, node_mrca_2)
          if(length(sub_tree_2$tip.label) == 2){
            tmp <- tibble(tree = i,pattern = "ann_deb",trio=trio_name)
            deb_positions <- rbind(deb_positions, tmp)
            
          }
          node_mrca_3 <- findMRCA(sub_tree, tips=c(arg_sample,deb_sample),type="node")
          sub_tree_3 <- extract.clade(sub_tree, node_mrca_3)
          if(length(sub_tree_3$tip.label) == 2){
            tmp <- tibble(tree = i,pattern = "arg_deb",trio=trio_name)
            deb_positions <- rbind(deb_positions, tmp)
          }
        }
      }
    }
  }
}

deb_positions %>%
  group_by(trio, pattern) %>%
  summarize(count=n()) %>%
  ggplot(.,aes(x=trio,y=count,fill=pattern)) +
  geom_bar(stat="identity",position="dodge") 
