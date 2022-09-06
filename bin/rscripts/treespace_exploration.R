library("treespace")
library("adegenet")
library("adegraphics")
library("rgl")
library(tidyverse)
library(cowplot)
library(ggtree)

distdex<-function(i,j,n){ #given row, column, and n, return index
  if(i==j){0
  }else if(i > j){
    n*(j-1) - j*(j-1)/2 + i-j
  }else{
    n*(i-1) - i*(i-1)/2 + j-i  
  }
}
all_trees <- read.tree("../vcf/gowens22.mappable.dp4.missing80.10kb.trees",keep.multi = T)
#all_trees <- read.tree("../vcf/tmp.trees",keep.multi = T)
tree_names <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10kb.treenames.txt",
                       col_names=c("name"))
names(all_trees) <- tree_names$name

all_trees_rooted <- root(all_trees,outgroup="664647_GIG",resolve.root = TRUE)
tree_pca <- treespace(all_trees_rooted,nf=10,processors = 10)
#saveRDS(tree_pca,file="gowens22.mappable.dp4.missing80.10kb.treespace.Robj")

as_tibble(tree_pca$pco$li,rownames = "loci") %>%
  separate(loci, c("chr","start","end"),"-",convert=T) %>%
  ggplot(.,aes(x=A1,y=A2,color=chr)) + 
  geom_point() +
  theme_cowplot() +
  facet_wrap(~chr)

as_tibble(tree_pca$pco$li,rownames = "loci") %>%
  separate(loci, c("chr","start","end"),"-",convert=T) %>%
  ggplot(.,aes(y=A5,x=start,color=chr)) + 
  geom_point() +
  theme_cowplot() +
  facet_wrap(~chr)

plotGroves(tree_pca$pco, lab.cex=1.5,xax=4,yax=3)
tree_pca.groves <- findGroves(tree_pca)
plotGrovesD3(tree_pca.groves,xax=1,yax=2)
as_tibble(tree_pca.groves$groups) %>%
  cbind(tree_names) %>%
  separate(name, c("chr","start","end"),"-",convert=T) %>%
  filter(chr == "Ha412HOChr17") %>%
  ggplot(.,aes(x=start,xend=end, yend=value, y=value,color=value)) + 
  geom_segment(size=3) +
  theme_cowplot()

as_tibble(tree_pca.groves$groups) %>%
  cbind(tree_names) %>%
  separate(name, c("chr","start","end"),"-",convert=T) %>%
  arrange(chr, start)

medtrees <- medTree(all_trees_rooted, tree_pca.groves$groups)

med.trees <- lapply(medtrees, function(e) ladderize(e$trees[[1]]))
par(mfrow=c(2,3))
for(i in 1:length(med.trees)) plot(med.trees[[i]], main=paste("cluster",i),cex=1.5)


######
ggdensitree(all_trees, alpha=.3, colour='steelblue') + 
  geom_tiplab(size=3) + xlim(0, 45)


###################
tree_pca.ann <- treespace(all_trees_rooted,nf=5,emphasise.tips=c("ANN1283","ANN1029"),emphasise.weight=3)
plotGroves(tree_pca.ann$pco, lab.cex=1.5,xax=1,yax=2)
tree_pca.ann.groves <- findGroves(tree_pca.ann, nclust=10)
plotGrovesD3(tree_pca.ann.groves,xax=1,yax=2)

medtrees <- medTree(all_trees_rooted, tree_pca.ann.groves$groups)

med.trees <- lapply(medtrees, function(e) ladderize(e$trees[[1]]))
par(mfrow=c(2,5))
for(i in 1:length(med.trees)) plot(med.trees[[i]], main=paste("cluster",i),cex=1.5)
plot(med.trees[[10]])
