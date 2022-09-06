library(phytools)
library(tidyverse)
library(TreeTools)
library(castor)
library(textclean)

tree_names <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10.2kb.physical.subsample.treenames.trees.txt",
                       col_names = c("tree_id","tree"))

tree_names$tree <- mgsub(tree_names$tree, c("NHKG_16_01_R_196", "GO8_GB179","PAR_posas_01","ANN1283","ARG0295","DEB_1135","PET0662","PET0424","PET0495","664647_GIG"), 
                         c("ano", "des","par","ann","arg","deb","niv", "pet_fal","pet_pet","outgroup"))


niv_positions <- tibble()
for (i in 1:nrow(tree_names)){
  print(i)
  #Test for pet_pet
  tmp_tree <- ape::read.tree(text=tree_names$tree[i])
  tmp_tree <- ape::drop.tip(tmp_tree, c("par", "deb","ann","arg"))
  tmp_tree <- root(tmp_tree, outgroup = "outgroup", resolve.root = TRUE)
  node_mrca <- findMRCA(tmp_tree, tips=c("niv","pet_fal","pet_pet"),type="node")
  sub_tree <- extract.clade(tmp_tree, node_mrca)
  node_mrca_1 <- findMRCA(tmp_tree, tips=c("niv","ano","des"),type="node")
  sub_tree_1 <- extract.clade(tmp_tree, node_mrca_1)
  if(length(sub_tree$tip.label) == 3){
    tmp <- tibble(tree = i, species = "niv",pattern = "pet_group")
    niv_positions <- rbind(niv_positions, tmp)
  }else if(length(sub_tree_1$tip.label) == 3){
    tmp <- tibble(tree = i, species = "niv",pattern = "ano_des_group")
    niv_positions <- rbind(niv_positions, tmp)
  }else{
    tmp <- tibble(tree = i, species = "niv",pattern = "other")
    niv_positions <- rbind(niv_positions, tmp)
  }
}


pdf("plots/niveus_position.pdf",height=4,width=4)
niv_positions %>%
  mutate(species = "H. niveus") %>%
  mutate(pattern = case_when(pattern == "ano_des_group" ~ "niv+ano/des",
                             pattern == "pet_group" ~ "niv+pet",
                              TRUE ~ "other")) %>%
  group_by(species, pattern) %>% 
  summarize(count=n()) %>%
  ggplot(.,aes(x=species,y=count,fill=pattern)) +
  geom_bar(stat="identity",position="dodge") +
  theme_cowplot() +
  ylab("Tree count") +
  scale_fill_manual(values=c("#e63946","#a8dadc","#1d3557"),
                    name="Pattern") +
  xlab("")
dev.off()
