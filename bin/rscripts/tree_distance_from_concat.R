#This is for calculating the distance in weird tree space between the concatenated tree and everyone else.
library("treespace")
library("adegenet")
library("adegraphics")
library("rgl")
library(tidyverse)
library(cowplot)

distdex<-function(i,j,n){ #given row, column, and n, return index
  if(i==j){0
  }else if(i > j){
    n*(j-1) - j*(j-1)/2 + i-j
  }else{
    n*(i-1) - i*(i-1)/2 + j-i  
  }
}
distdex(100,1,3426)
#all_trees <- read.tree("../vcf/gowens22.mappable.dp4.missing80.10kb.startingwithconcat.trees",keep.multi = T)
tree_names <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10kb.treenames.txt",
                       col_names=c("name"))
all_trees_rooted <- root(all_trees,outgroup="664647_GIG",resolve.root = TRUE)
#tree_pca <- treespace(all_trees_rooted,nf=10,processors = 10)
saveRDS(tree_pca,file="gowens22.mappable.dp4.missing80.10kb.startingwithconcat.treespace.Robj")
dist_from_concat <- tibble(dist=tree_pca$D[1:3425], id=tree_names$name)

dist_from_concat %>%
  separate(id, c("chr","start","end"),"-",convert=T) %>%
  ggplot(.,aes(x=start,y=dist)) +
  geom_point() +
  facet_wrap(~chr)

#Add recombination rate
genetic_map <- read_tsv("/home/owens/ref/Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt") %>% 
  mutate(rate = lead(cM)-cM) %>% 
  rename(start = pos) %>% mutate(end = start + 999999) %>%
  mutate(rate = case_when(is.na(rate) ~ cM - lag(cM),
                          rate < 0 ~ cM - lag(cM),
                          TRUE ~ rate))
dist_from_concat %>%
  separate(id, c("chr","start","end"),"-",convert=T) %>%
  select(chr,start, end) %>% unique() -> tree_windows

tree_windows$cm_rate <- NA
for (i in 1:nrow(tree_windows)){
  print(i)
  chosen_start <- tree_windows$start[i]
  chosen_end <-tree_windows$end[i]
  chosen_chr <- tree_windows$chr[i]
  start_rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
    filter(start <= chosen_start, end >= chosen_start) %>%
    pull(rate)
  end_rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
    filter(start <= chosen_end, end >= chosen_end) %>%
    pull(rate)
  mean_rate <- mean(start_rate,end_rate)
  tree_windows$cm_rate[i] <- mean_rate
}  
tree_windows %>%
  mutate(quantile_rank = ntile(tree_windows$cm_rate,5)) %>%
  mutate(id = paste(chr,start,end,sep="-"))-> tree_windows

dist_from_concat %>%
  inner_join(tree_windows) %>%
  ggplot(.,aes(x=as.factor(quantile_rank),y=dist)) + geom_boxplot()

dist_from_concat %>%
  inner_join(tree_windows) %>%
  lm(data=., dist ~ cm_rate) %>%
  summary()
