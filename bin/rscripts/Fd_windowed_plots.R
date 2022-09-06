library(tidyverse)
library(cowplot)
dstats <- read_tsv("/media/owens/Childs/helianthus_phylogeny/d/deb_pet_pet_ann_localFstats_test_50_25.txt")

dstats %>%
  mutate(chr_n = gsub("Ha412HOChr","",chr)) %>%
  #filter(chr_n == 13) %>%
  ggplot(.,aes(x=windowStart,y=f_d)) +
  geom_point(aes(color=f_d)) +
  facet_wrap(~chr_n) +
  geom_hline(yintercept=0,linetype="dotted") +
  theme_cowplot() 
