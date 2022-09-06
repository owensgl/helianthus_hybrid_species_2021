library(tidyverse)
library(cowplot)

data <- read_tsv("../vcf/gowens22.mappable.variant.snps.filtered.dp4.missing80.commonderived.txt.gz")

data %>% 
  filter(sample == "ANN1283") %>%
  mutate(common_derived_present =case_when(common_derived_present == 0.5 ~ 1,
                                           TRUE ~ common_derived_present)) %>% 
  group_by(sample,chr) %>%
  mutate(mean_derived = rollapply(common_derived_present, 50, mean, fill = NA, align = "right")) %>%
  ggplot(.,aes(x=pos/1000000,y=mean_derived)) +
  geom_line() +
  facet_wrap(~chr) +
  theme_cowplot() +
  xlab("MBp") +
  ylab("Proportion_derived")
  
#dbinom(30, size=50, prob=0.799)
data %>% 
  mutate(common_derived_present =case_when(common_derived_present == 0.5 ~ 1,
                                           TRUE ~ common_derived_present)) %>% 
  group_by(sample) %>%
  summarize(mean = mean(common_derived_present))

data %>% 
  filter(sample == "DEB_1837") %>%
  #filter(chr == "Ha412HOChr07") %>%
  mutate(common_derived_present =case_when(common_derived_present == 0.5 ~ 1,
                                           TRUE ~ common_derived_present)) %>% 
  group_by(sample,chr) %>%
  mutate(mean_derived = rollapply(common_derived_present, 50, mean, fill = NA, align = "right")) %>%
  ggplot(.,aes(x=pos/1000000,y=mean_derived)) +
  geom_line() +
  facet_wrap(~chr) +
  theme_cowplot() +
  xlab("MBp") +
  ylab("Proportion_derived")


data %>%
  filter(chr == "Ha412HOChr07") %>%
  filter(pos > 130000000, pos < 140000000) %>%
  ggplot(.,aes(x=as.factor(pos),y=sample,fill=common_derived_present)) + 
  geom_tile() +
  scale_fill_viridis_c() +
  ggtitle("Chr07, 130-140MBp")

data %>%
  filter(chr == "Ha412HOChr11") %>%
  filter(pos > 140000000, pos < 160000000) %>%
  ggplot(.,aes(x=as.factor(pos),y=sample,fill=common_derived_present)) + 
  geom_tile() +
  scale_fill_viridis_c() +
  ggtitle("Chr11, 140-160MBp")

data %>%
  filter(chr == "Ha412HOChr13") %>%
  filter(pos > 165000000, pos < 175000000) %>%
  ggplot(.,aes(x=as.factor(pos),y=sample,fill=common_derived_present)) + 
  geom_tile() +
  scale_fill_viridis_c() +
  ggtitle("Chr13, 165-175MBp")

data %>%
  filter(chr == "Ha412HOChr17") %>%
  filter(pos > 60000000, pos < 80000000) %>%
  ggplot(.,aes(x=as.factor(pos),y=sample,fill=common_derived_present)) + 
  geom_tile() +
  scale_fill_viridis_c() +
  ggtitle("Chr07, 130-140MBp")
