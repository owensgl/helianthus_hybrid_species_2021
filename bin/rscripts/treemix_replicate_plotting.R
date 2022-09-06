#Compare the likelihood of treemix replicates 
source("rscripts/treemix_plotting_tidy_functions.R")
#Load up treemix_plotting_tidy_functions.R
library(tidyverse)
library(patchwork)
library(cowplot)
library(ggbeeswarm)

reps <- read_tsv("../treemix_noss_reps/allreps.txt",col_names=F) %>%
  pull(X1)
reps <- unique(reps)

likelihoods <- tibble()
for (i in seq(0,10)){
  for (rep in reps){
    likelihood_chosen <- read_delim(paste0("/media/drive_5_usb/Childs/helianthus_phylogeny/treemix_noss_reps/regular/owens22.mappable.variant.snps.filtered.dp4.missing80.pruned.",rep,".",i,".llik"),
                                    delim=" ")[[7]]
    tmp <- tibble(migration_nodes=i,likelihood=likelihood_chosen,rep=rep)
    likelihoods <- rbind(likelihoods, tmp)
  }
}

likelihoods <- unique(likelihoods)
likelihoods %>%
  write_tsv("TreeMix_likelihoods.txt")
likelihoods %>%
  ggplot(.,aes(x=as.factor(migration_nodes),y=likelihood)) + geom_jitter(width=0.2) + theme_cowplot()
pdf("plots/treemix_likelihoods.v1.pdf")
likelihoods %>%
  #filter(migration_nodes >= 4) %>%
  ggplot(.,aes(x=as.factor(migration_nodes),y=likelihood)) + 
  geom_quasirandom(method='pseudorandom',alpha=.2)+ 
  theme_cowplot() +
  xlab("Migration nodes") +
  ylab("Likelihood")
  #facet_wrap(~migration_nodes,scales="free" )
dev.off()

max_likelihoods <- likelihoods %>%
  group_by(migration_nodes) %>%
  summarize(max_lik = max(likelihood)) 
likelihood_change <- tibble()
for (i in 1:10){
  l_dif = max_likelihoods$max_lik[i+1] - max_likelihoods$max_lik[i]
  likelihood_change <- rbind(likelihood_change, tibble(migration_nodes = i, l_dif = l_dif))
}
pdf("plots/treemix_likelihood_improvement.v1.pdf")
likelihood_change %>%
  ggplot(.,aes(x=as.factor(migration_nodes),y=l_dif)) + 
  geom_point() + 
  theme_cowplot() +
  ylab("Likelihood improvement") +
  xlab("Migration nodes")
dev.off()

top_reps <- likelihoods %>%
  inner_join(max_likelihoods) %>%
  filter(likelihood >= (max_lik-10)) %>%
  group_by(migration_nodes) %>%
  distinct(likelihood) %>%
  slice_head(n=5) 

top_reps <- likelihoods %>%
  inner_join(top_reps) %>%
  group_by(migration_nodes,likelihood) %>%
  slice_head(n=1) %>% 
  arrange(migration_nodes, desc(likelihood))

pdf("plots/treemix_top_5.pdf",height=4,width=4)
for (n in 1:nrow(top_reps)){
  i <- top_reps$migration_nodes[n]
  rep_n <- top_reps$rep[n]
  lik <- top_reps$likelihood[n]

  treemix <- read_treemix(paste0("/media/drive_5_usb/Childs/helianthus_phylogeny/treemix_noss_reps/regular/owens22.mappable.variant.snps.filtered.dp4.missing80.pruned.",rep_n,".",i))
  treemix_tree <- plot_treemix(treemix) + ylab("") + xlab("Drift parameter") +
    theme_cowplot() + 
    theme(axis.ticks.y = element_blank(),
          axis.text.y=element_blank()) +
    scale_colour_gradient2("Migration\nweight", high = "red", mid = "yellow") +
    scale_x_continuous(expand = c(0, .01)) + ggtitle(paste0("Mig:",i, " ","Likelihood:",lik))
  print(treemix_tree)
}
dev.off()
all_plots <- list()
p <- 1
for (n in 1:nrow(top_reps)){
  i <- top_reps$migration_nodes[n]
  rep_n <- top_reps$rep[n]
  lik <- top_reps$likelihood[n]
  
  treemix <- read_treemix(paste0("/media/drive_5_usb/Childs/helianthus_phylogeny/treemix_noss_reps/regular/owens22.mappable.variant.snps.filtered.dp4.missing80.pruned.",rep_n,".",i))
  treemix_tree <- plot_treemix(treemix) + ylab("") + xlab("Drift parameter") +
    theme_cowplot() + 
    theme(axis.ticks.y = element_blank(),
          axis.text.y=element_blank()) +
    scale_colour_gradient2("Migration\nweight", high = "red", mid = "yellow") +
    scale_x_continuous(expand = c(0, .01)) + ggtitle(paste0("Mig:",i, " ","Likelihood:",lik))
  all_plots[[p]] <- treemix_tree
  p <- p+ 1
}
#43 plots in total

pdf("plots/treemix_top_5_ordered.pdf",height=10,width=10)

all_plots[[1]] + all_plots[[2]] + all_plots[[3]] + all_plots[[4]] + 
  all_plots[[5]] + all_plots[[6]] + all_plots[[7]] + all_plots[[8]] +
  all_plots[[9]] 
all_plots[[10]] + all_plots[[11]] + all_plots[[12]] +
  all_plots[[13]] + all_plots[[14]] + all_plots[[15]] + 
  all_plots[[16]] + all_plots[[17]] + all_plots[[18]] 

all_plots[[19]] + all_plots[[20]] + all_plots[[21]] +
  all_plots[[22]] + all_plots[[23]] + all_plots[[24]] + 
  all_plots[[25]] + all_plots[[26]] + all_plots[[27]] 

all_plots[[28]] + all_plots[[29]] + all_plots[[30]] +
  all_plots[[31]] + all_plots[[32]] + all_plots[[33]] + 
  all_plots[[34]] + all_plots[[35]] + all_plots[[36]] 

all_plots[[37]] + all_plots[[38]] + all_plots[[39]] +
  all_plots[[40]] + all_plots[[41]] + all_plots[[42]] + 
  all_plots[[43]]  
dev.off()
top_5_reps <- likelihoods %>%
  inner_join(max_likelihoods) %>%
  #filter(likelihood == max_lik) %>%
  group_by(migration_nodes) %>%
  distinct(likelihood,.keep_all=T) %>%
  arrange(desc(likelihood)) %>%
  slice_head(n=5) 

top_5_reps_3 <- top_5_reps %>% 
  filter(migration_nodes == 3) 

treemix_trees_3nodes <- list()
for (i in seq(0,4)){
  rep_n <- top_5_reps_3$rep[i+1]
  treemix <- read_treemix(paste0("/media/owens/Childs/helianthus_phylogeny/treemix_noss_reps/owens22.mappable.variant.snps.filtered.dp4.missing80.pruned.",rep_n,".",3))
  treemix_tree <- plot_treemix(treemix) + ylab("") + xlab("Drift parameter") +
    theme_cowplot() + 
    theme(axis.ticks.y = element_blank(),
          axis.text.y=element_blank()) +
    scale_colour_gradient2("Migration\nweight", high = "red", mid = "yellow") +
    scale_x_continuous(expand = c(0, .01)) + ggtitle(paste0("Mig:",i))
  treemix_trees_3nodes[[i+1]] <- treemix_tree
}
treemix_trees_3nodes[1]
