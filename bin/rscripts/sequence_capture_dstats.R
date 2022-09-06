#This is for plotting D stats from the sequence capture section
library(tidyverse)
library(cowplot)
library(ggbeeswarm)

data <- read_tsv("../sequence_capture/genes.dsuite_tree.txt")

data %>%
  rename(P1 = P2, P2 = P1) %>%
  mutate(Dstatistic = -Dstatistic) -> data_flip

rbind(data, data_flip) -> complete_data

large_perennials <- c("H_salicifolius","H_maximiliani","H_giganteus",
                      "H_verticillatus","H_grosseserratus","H_nuttalli",
                      "H_divaricatus","H_microcephalus","H_cusickii",
                      "H_arizonensis","H_laciniatus")
south_perennials <- c("H_longifolius","H_carnosus","H_radula",
                      "H_atrorubens","H_silphioides","H_heterophyllus",
                      "H_angustifolius","H_floridanus")
other_perennials <- c("H_mollis","H_occidentalis","H_gracilentus","H_agrestis")
outgroup_perennials <- c("H_porteri")
annuus_group <- c("H_annuus","H_argophyllus")
exilis_group <- c("H_exilis")
petiolaris_group <- c("H_niveus","H_petiolaris","H_debilis","H_praecox")

#Exilis compared to annuus/arg against perennials
complete_data %>%
  filter(P1 %in% annuus_group) %>%
  filter(!P2 %in% annuus_group) %>%
  filter(P2 %in% exilis_group) %>%
  filter(!P3 %in% petiolaris_group) %>%
  mutate(ID = row_number()) %>% 
  mutate(group_3 = case_when(P3 %in% large_perennials ~ "large_perennials",
                             P3 %in% south_perennials ~ "south_perennials",
                             P3 %in% other_perennials ~ "other_perennials",
                             P3 %in% outgroup_perennials ~ "outgroup_perennials")) %>% 
  mutate(sig = case_when(`p-value` < 0.05 ~ "significant",
                         TRUE ~ "non-significant")) %>%
  ggplot(.,aes(x=group_3,y=Dstatistic)) + 
  geom_jitter(aes(color=P1,shape=sig),width=0.2)  +
  theme_cowplot() +
  geom_hline(yintercept=0,linetype="dotted")

#petiolaris clade compared to annuus/arg against perennials

complete_data %>%
  filter(P1 %in% annuus_group) %>%
  filter(!P2 %in% annuus_group) %>%
  filter(P2 %in% petiolaris_group) %>%
  filter(!P3 %in% exilis_group) %>%
  mutate(ID = row_number()) %>% 
  mutate(group_3 = case_when(P3 %in% large_perennials ~ "large_perennials",
                             P3 %in% south_perennials ~ "south_perennials",
                             P3 %in% other_perennials ~ "other_perennials",
                             P3 %in% outgroup_perennials ~ "outgroup_perennials")) %>% 
  mutate(sig = case_when(`p-value` < 0.05 ~ "significant",
                         TRUE ~ "non-significant")) %>% 
  ggplot(.,aes(x=group_3,y=Dstatistic)) + 
  geom_beeswarm(aes(color=P1),width=1,dodge.width=1)  +
  theme_cowplot() +
  geom_hline(yintercept=0,linetype="dotted") +
  geom_boxplot(aes(color=P1),outlier.shape=NA,fill=NA)
  
complete_data %>%
  filter(P1 %in% exilis_group) %>%
  filter(!P2 %in% annuus_group) %>%
  filter(P2 %in% petiolaris_group) %>%
  filter(!P3 %in% annuus_group) %>%
  mutate(ID = row_number()) %>% 
  mutate(group_3 = case_when(P3 %in% large_perennials ~ "large_perennials",
                             P3 %in% south_perennials ~ "south_perennials",
                             P3 %in% other_perennials ~ "other_perennials",
                             P3 %in% outgroup_perennials ~ "outgroup_perennials")) %>% 
  mutate(sig = case_when(`p-value` < 0.05 ~ "significant",
                         TRUE ~ "non-significant")) %>% 
ggplot(.,aes(x=group_3,y=Dstatistic)) + 
  geom_beeswarm(aes(color=group_3,shape=sig),width=0.2)  +
  theme_cowplot() +
  geom_hline(yintercept=0,linetype="dotted")
  