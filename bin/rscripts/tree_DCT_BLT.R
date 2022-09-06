library(phytools)
library(tidyverse)
library(TreeTools)
library(castor)
library(textclean)

tree_names <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10.2kb.physical.subsample.treenames.trees.txt",
                       col_names = c("tree_id","tree"))

tree_names$tree <- mgsub(tree_names$tree, c("NHKG_16_01_R_196", "GO8_GB179","PAR_posas_01","ANN1283","ARG0295","DEB_1135","PET0662","PET0424","PET0495","664647_GIG"), 
                      c("ano", "des","par","ann","arg","deb","niv", "pet_fal","pet_pet","outgroup"))

#Use dsuite to select sets of samples to set
dsuite_results <- read_tsv("../d/gowens22.mappable.variant.snps.filtered.dp4.missing80.tree_1_tree.txt") %>%
  select(P1,P2,P3)
total_results_rows = nrow(dsuite_results) * nrow(tree_names)
test_results <- matrix(nrow = total_results_rows, ncol = 8)
counter = 1
for (x in 1:nrow(dsuite_results)){
  print(x)
  sample_1 <- dsuite_results$P1[x]
  sample_2 <- dsuite_results$P2[x]
  sample_3 <- dsuite_results$P3[x]
  for (i in 1:nrow(tree_names)){
    tmp_tree <- ape::read.tree(text=tree_names$tree[i])
    samples_to_find <- c(sample_1,sample_2,sample_3)
    tmp_tree <- Preorder(tmp_tree)
    sub_tree <- get_subtree_with_tips(tmp_tree,samples_to_find)
    tmp_tree <- Preorder(sub_tree$subtree)
    
    test_12 <- Subtree(tmp_tree,findMRCA(tmp_tree,c(sample_1,sample_2)))
    if (test_12$Nnode == 1){
      result <- "12"
    }
    test_13 <- Subtree(tmp_tree,findMRCA(tmp_tree,c(sample_1,sample_3)))
    if (test_13$Nnode == 1){
      result <- "13"
    }
    test_23 <- Subtree(tmp_tree,findMRCA(tmp_tree,c(sample_2,sample_3)))
    if (test_23$Nnode == 1){
      result <- "23"
    }
    
    tree_matrix <- cophenetic.phylo(tmp_tree)

    
    dist_12 <- tree_matrix[sample_1,sample_2]
    dist_13 <- tree_matrix[sample_1,sample_3]
    dist_23 <- tree_matrix[sample_2,sample_3]
    
    
    test_results[counter, ] <- c(sample_1, sample_2, sample_3, tree_names$tree_id[i],
                                 dist_12,dist_13,dist_23,result)
    
    
    
    counter <- counter+1
  }

}

#DCT from 155 drosophila paper
test_tibble <- as_tibble(test_results)
colnames(test_tibble) <- c("P1","P2","P3","tree","dist_12","dist_13","dist_23","topology")
test_tibble %>%
  group_by(P1,P2,P3, topology) %>%
  count() %>%
  pivot_wider(names_from = topology,values_from = n) -> test_counts

write_tsv(test_tibble,"gowens22.mappable.dp4.missing80.10.2kb.physical.subsample.topology.txt")
test_counts$dct_pvalue <- NA
for (i in 1:nrow(test_counts)){
  test_counts$dct_pvalue[i] <- chisq.test(c(test_counts$`13`[i],test_counts$`23`[i]))$p.value
}
test_counts$dct_pvalue_hoch <- p.adjust(test_counts$dct_pvalue, method="hochberg")
write_tsv(test_counts,"gowens22.mappable.dp4.missing80.10.2kb.physical.subsample.dct.results.txt")

#BLT from 155 drosophila paper
test_tibble %>%
  filter(topology != "12") %>%
  mutate(dtype = topology, d = case_when(topology == 23 ~ dist_23,
                                            topology == 13 ~ dist_13)) %>%
  mutate(d = as.numeric(d)) -> test_blt

blt_result <- tibble()
for ( x in 1:nrow(dsuite_results)){
  print(x)
  sample_1 <- dsuite_results$P1[x]
  sample_2 <- dsuite_results$P2[x]
  sample_3 <- dsuite_results$P3[x]
  test_blt %>%
    filter(P1 == sample_1,
           P2 == sample_2,
           P3 == sample_3) -> tmp
  tmp %>%
    group_by(dtype) %>%
    summarize(mean_d = mean(d)) %>%
    pivot_wider(names_from = dtype,values_from = mean_d) %>%
    mutate(dif = `13` - `23`) %>% pull(dif) -> dif
  
  a <- tmp %>% filter(dtype == "13") %>% pull(d)
  b <- tmp %>% filter(dtype == "23") %>% pull(d)
  result <- wilcox.test(a,b) 
  result_tmp <- tibble(P1 = sample_1, P2 = sample_2, P3 = sample_3,
                       blt_pvalue = result$p.value, 
                       blt_dif = dif)
  blt_result <- rbind(blt_result, result_tmp)
}
blt_result$blt_pvalue_hoch <- p.adjust(blt_result$blt_pvalue, method="hochberg")
write_tsv(blt_result,"gowens22.mappable.dp4.missing80.10.2kb.physical.subsample.blt.results.txt")


inner_join(blt_result, test_counts) %>% 
  group_by(P1,P2,P3) %>%
  summarize(blt_result = case_when(blt_pvalue_hoch > 0.05 ~ "nonsig",
                                blt_dif > 0 ~ "P2-P3",
                                blt_dif < 0 ~ "P1-P3"),
         dct_result = case_when(dct_pvalue_hoch > 0.05 ~ "nonsig",
                                `13` > `23` ~ "P1-P3",
                                `13` < `23` ~ "P2-P3")) %>% View()
  pivot_longer(-c(P1,P2,P3),names_to="test",values_to="result") %>%
  mutate(trio = paste(P1,P2,P3,sep=".")) %>%
  ggplot(.,aes(x=test,y=trio,fill=result)) +
  geom_tile()
  
  

introgression_summary_dct <- tibble()
for (i in 1:nrow(test_counts)){
  print(i)
  P1 <- test_counts$P1[i]
  P2 <- test_counts$P2[i]
  P3 <- test_counts$P3[i]
  p <- test_counts$dct_pvalue_hoch[i]
  `13` <- test_counts$`13`
  `23` <- test_counts$`23`
  
  if (p < 0.05 & `13` < `23`){
    tmp <- tibble(spe1 = P2, spe2 = P3, signal=1)
    introgression_summary_dct <- rbind(introgression_summary_dct, tmp)
    tmp <- tibble(spe1 = P3, spe2 = P2, signal=1)
    introgression_summary_dct <- rbind(introgression_summary_dct, tmp)
    tmp <- tibble(spe1 = P1, spe2 = P3, signal=0)
    introgression_summary_dct <- rbind(introgression_summary_dct, tmp)
    tmp <- tibble(spe1 = P3, spe2 = P1, signal=0)
    introgression_summary_dct <- rbind(introgression_summary_dct, tmp)
  }else if(p < 0.05 & `13` > `23`){
    tmp <- tibble(spe1 = P1, spe2 = P3, signal=1)
    introgression_summary_dct <- rbind(introgression_summary_dct, tmp)
    tmp <- tibble(spe1 = P3, spe2 = P1, signal=1)
    introgression_summary_dct <- rbind(introgression_summary_dct, tmp)
    tmp <- tibble(spe1 = P2, spe2 = P3, signal=0)
    introgression_summary_dct <- rbind(introgression_summary_dct, tmp)
    tmp <- tibble(spe1 = P3, spe2 = P2, signal=0)
    introgression_summary_dct <- rbind(introgression_summary_dct, tmp)
  }else{
    tmp <- tibble(spe1 = P1, spe2 = P3, signal=0)
    introgression_summary_dct <- rbind(introgression_summary_dct, tmp)
    tmp <- tibble(spe1 = P3, spe2 = P1, signal=0)
    introgression_summary_dct <- rbind(introgression_summary_dct, tmp)
    tmp <- tibble(spe1 = P2, spe2 = P3, signal=0)
    introgression_summary_dct <- rbind(introgression_summary_dct, tmp)
    tmp <- tibble(spe1 = P3, spe2 = P2, signal=0)
    introgression_summary_dct <- rbind(introgression_summary_dct, tmp)
  }

}
  
  

introgression_summary_blt <- tibble()
for (i in 1:nrow(blt_result)){
  print(i)
  P1 <- blt_result$P1[i]
  P2 <- blt_result$P2[i]
  P3 <- blt_result$P3[i]
  p <- blt_result$blt_pvalue_hoch[i]
  blt_dif <- blt_result$blt_dif
  if (p < 0.05 & blt_dif > 0){
    tmp <- tibble(spe1 = P2, spe2 = P3, signal=1)
    introgression_summary_blt <- rbind(introgression_summary_blt, tmp)
    tmp <- tibble(spe1 = P3, spe2 = P2, signal=1)
    introgression_summary_blt <- rbind(introgression_summary_blt, tmp)
    tmp <- tibble(spe1 = P1, spe2 = P3, signal=0)
    introgression_summary_blt <- rbind(introgression_summary_blt, tmp)
    tmp <- tibble(spe1 = P3, spe2 = P1, signal=0)
    introgression_summary_blt <- rbind(introgression_summary_blt, tmp)
  }else if(p < 0.05 & blt_dif < 0){
    tmp <- tibble(spe1 = P1, spe2 = P3, signal=1)
    introgression_summary_blt <- rbind(introgression_summary_blt, tmp)
    tmp <- tibble(spe1 = P3, spe2 = P1, signal=1)
    introgression_summary_blt <- rbind(introgression_summary_blt, tmp)
    tmp <- tibble(spe1 = P2, spe2 = P3, signal=0)
    introgression_summary_blt <- rbind(introgression_summary_blt, tmp)
    tmp <- tibble(spe1 = P3, spe2 = P2, signal=0)
    introgression_summary_blt <- rbind(introgression_summary_blt, tmp)
  }else{
    tmp <- tibble(spe1 = P1, spe2 = P3, signal=0)
    introgression_summary_blt <- rbind(introgression_summary_blt, tmp)
    tmp <- tibble(spe1 = P3, spe2 = P1, signal=0)
    introgression_summary_blt <- rbind(introgression_summary_blt, tmp)
    tmp <- tibble(spe1 = P2, spe2 = P3, signal=0)
    introgression_summary_blt <- rbind(introgression_summary_blt, tmp)
    tmp <- tibble(spe1 = P3, spe2 = P2, signal=0)
    introgression_summary_blt <- rbind(introgression_summary_blt, tmp)
  }
  
}

species_order <- tibble(spe = c("par","arg","ann","ano","des","niv","deb","pet_fal","pet_pet"),
                        rank = 1:9)
introgression_summary_blt %>%
  group_by(spe1,spe2) %>%
  summarize(introgression_signal = mean(signal), tests=n(), proportion=paste0(tests*introgression_signal,"/", tests)) %>%
  
  inner_join(species_order %>% rename(spe1 = spe,rank1=rank)) %>%
  inner_join(species_order %>% rename(spe2 = spe,rank2=rank)) %>%
  filter(rank1 > rank2) %>%
  ggplot(.,aes(x=fct_reorder(spe1,-rank1),y=fct_reorder(spe2,-rank2),fill=introgression_signal)) +
  geom_tile(color="black",size=2) +
  scale_fill_viridis_c() +
  theme_cowplot() +
  scale_fill_viridis_c(name="BLT introgression\nsignal") +
  theme_cowplot() +
  xlab("") +
  ylab("") +
  geom_label(aes(label=proportion),fill="white")


introgression_summary_dct %>%
  group_by(spe1,spe2) %>%
  summarize(introgression_signal = mean(signal), tests=n(), proportion=paste0(tests*introgression_signal,"/", tests)) %>%
  inner_join(species_order %>% rename(spe1 = spe,rank1=rank)) %>%
  inner_join(species_order %>% rename(spe2 = spe,rank2=rank)) %>%
  filter(rank1 > rank2) %>%
  ggplot(.,aes(x=fct_reorder(spe1,-rank1),y=fct_reorder(spe2,-rank2),fill=introgression_signal)) +
  geom_tile(color="black",size=2) +
  scale_fill_viridis_c() +
  theme_cowplot() +
  scale_fill_viridis_c(name="DCT introgression\nsignal") +
  theme_cowplot() +
  xlab("") +
  ylab("") +
  geom_label(aes(label=proportion),fill="white")
  
test_counts %>%
  mutate(sig = case_when(dct_pvalue_hoch < 0.05 ~ "sig",
                         TRUE ~ "non-sig")) %>%
  group_by(sig) %>%
  summarize(n=n()) %>%
  mutate(freq = n / sum(n)) %>%
  ggplot(.,aes(x="Dstat",y=freq,fill=sig)) +
  geom_bar(stat="identity") +
  theme_cowplot() +
  scale_fill_brewer(palette = "Set1",name="Significant?") +
  xlab("") + ylab("Proportion of trios")
blt_result %>%
  mutate(sig = case_when(blt_pvalue_hoch < 0.05 ~ "sig",
                         TRUE ~ "non-sig")) %>%
  group_by(sig) %>%
  summarize(n=n()) %>%
  mutate(freq = n / sum(n)) %>%
  ggplot(.,aes(x="Dstat",y=freq,fill=sig)) +
  geom_bar(stat="identity") +
  theme_cowplot() +
  scale_fill_brewer(palette = "Set1",name="Significant?") +
  xlab("") + ylab("Proportion of trios")
