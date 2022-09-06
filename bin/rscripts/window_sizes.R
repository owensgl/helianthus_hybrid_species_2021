#Plot of window size

library(tidyverse)
library(cowplot)
trees <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10kb.windowsize.txt",
                  col_names = c("chr","start","end","length"))

trees %>%
  mutate(middle = (start+end)/2) %>%
  ggplot(.,aes(length/1000000)) + geom_histogram() +
  theme_cowplot() +
  xlab("Mb")

trees %>%
  mutate(middle = (start+end)/2) %>%
  ggplot(.,aes(x=middle/1000000,y=length)) +
  facet_wrap(~chr,ncol=1) + 
  geom_point(aes(color=length/1000000)) + 
  scale_color_viridis_c(name="window size\n(Mb)") +
  theme_cowplot() +
  xlab("Mb")
