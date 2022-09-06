source("rscripts/treemix_plotting_tidy_functions.R")
#Load up treemix_plotting_tidy_functions.R
library(tidyverse)
library(patchwork)
library(cowplot)
treemix_trees <- list()

treemix <- read_treemix(paste0("/media/drive_5_usb/Childs/helianthus_phylogeny/treemix/owens22.mappable.variant.snps.filtered.dp4.missing80.pruned.0"))
treemix_tree <- plot_treemix(treemix) + ylab("") + xlab("Drift parameter") +
  theme_cowplot() + 
  theme(axis.ticks.y = element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_gradient2("Migration\nweight", high = "red", mid = "yellow") +
  scale_x_continuous(expand = c(0, .01)) + ggtitle("Mig:0")
treemix_trees[[11]] <- treemix_tree
for (i in seq(1,10)){
  treemix <- read_treemix(paste0("/media/drive_5_usb/Childs/helianthus_phylogeny/treemix/owens22.mappable.variant.snps.filtered.dp4.missing80.pruned.",i))
  treemix_tree <- plot_treemix(treemix) + ylab("") + xlab("Drift parameter") +
    theme_cowplot() + 
    theme(axis.ticks.y = element_blank(),
          axis.text.y=element_blank()) +
    scale_colour_gradient2("Migration\nweight", high = "red", mid = "yellow") +
    scale_x_continuous(expand = c(0, .01)) + ggtitle(paste0("Mig:",i))
  treemix_trees[[i]] <- treemix_tree
}

#Treemix likelihoods scores
likelihoods <- tibble(migration_nodes=numeric(),likelihood=numeric())
likelihood_chosen <- read_delim(paste0("/media/owens/Childs/helianthus_phylogeny/treemix/owens22.mappable.variant.snps.filtered.dp4.missing80.pruned.0.llik"),
                                delim=" ")[[7]]
tmp <- tibble(migration_nodes=0,likelihood=likelihood_chosen)
likelihoods <- rbind(likelihoods, tmp)
for (i in seq(1,10)){
  likelihood_chosen <- read_delim(paste0("/media/owens/Childs/helianthus_phylogeny/treemix/owens22.mappable.variant.snps.filtered.dp4.missing80.pruned.",i,".llik"),
                                  delim=" ")[[7]]
  tmp <- tibble(migration_nodes=i,likelihood=likelihood_chosen)
  likelihoods <- rbind(likelihoods, tmp)
}
likelihoods %>%
  mutate(likelihood_gain = likelihood - lag(likelihood)) %>%
  filter(migration_nodes > 0) %>%
  ggplot(.,aes(x=migration_nodes,y=likelihood_gain)) + geom_line() +
  theme_cowplot()+ 
  ylab("Likelihood difference") + xlab("Migration edges") +
  scale_x_continuous(breaks=seq(1,10)) + 
  geom_hline(yintercept = 0,linetype="dotted")

treemix_trees[[1]]


treemix <- read_treemix(paste0("/media/owens/Childs/helianthus_phylogeny/treemix/test"))
treemix_tree <- plot_treemix(treemix) + ylab("") + xlab("Drift parameter") +
  theme_cowplot() + 
  theme(axis.ticks.y = element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_gradient2("Migration\nweight", high = "red", mid = "yellow") +
  scale_x_continuous(expand = c(0, .01)) + ggtitle("Mig:0")
treemix_tree
