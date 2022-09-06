library(tidyverse)

data <- read_tsv("../d/gowens22.mappable.variant.snps.filtered.dp4.missing80_BBAA.txt")


data %>%
  filter(P1 == "ann",P2 == "par", P3 == "pet_fal")  %>%
  select(P1, BBAA, ABBA,BABA) %>%
  pivot_longer(-P1, names_to = "type",values_to="count") %>%
  arrange(desc(type)) %>%
  ggplot(.,aes(x=fct_rev(type),y=count,fill=type)) +
  geom_bar(stat="identity") +
  theme_cowplot() +
  xlab("") +
  scale_fill_manual(values=c("#4ED86A","#FBD030","#367EFB")) +
  theme(legend.position = "none")


data %>%
  filter(P1 == "ano",P2 == "pet_fal", P3 == "ann")  %>%
  select(P1, BBAA, ABBA,BABA) %>%
  pivot_longer(-P1, names_to = "type",values_to="count") %>%
  arrange(desc(type)) %>%
  ggplot(.,aes(x=fct_rev(type),y=count,fill=type)) +
  geom_bar(stat="identity") +
  theme_cowplot() +
  xlab("") +
  scale_fill_manual(values=c("#e63946","#a8dadc","#1d3557")) +
  theme(legend.position = "none")
ggsave("plots/gowens22.mappable.variant.snps.filtered.dp4.missing80_BBAA.ano.v1.pdf")

data %>%
  filter(P1 == "des",P2 == "pet_fal", P3 == "ann")  %>%
  select(P1, BBAA, ABBA,BABA) %>%
  rename(ABBA = BABA, BABA = ABBA) %>%
  rename(ABBA = BBAA, BBAA = ABBA) %>%
  
  pivot_longer(-P1, names_to = "type",values_to="count") %>%
  arrange(desc(type)) %>%
  ggplot(.,aes(x=fct_rev(type),y=count,fill=type)) +
  geom_bar(stat="identity") +
  theme_cowplot() +
  xlab("") +
  scale_fill_manual(values=c("#4ED86A","#FBD030","#367EFB")) +
  theme(legend.position = "none")
