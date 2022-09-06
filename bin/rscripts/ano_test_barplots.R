library(tidyverse)
library(RcppRoll)
data <- read_tsv("../d/gowens22.mappable.variant.snps.filtered.dp4.missing80_BBAA.txt")
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
  theme(legend.position = "none") +
  ylab("Site count")
ggsave("plots/gowens22.mappable.variant.snps.filtered.dp4.missing80_BBAA.ano.v1.pdf",width=12,height=4)

twisst <- read_tsv("../twisst/ann.petfal.ano.weights.txt",comment="#",
                   col_names=c("ann_pet","ann_ano","pet_ano"),skip=4)
tree_ids <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10.2kb.physical.treenames.trees.txt",
                     col_names=c("tree_id","tree"))

chisq.test(c(1078,720))$p.value
twisst %>%
  cbind(tree_ids) %>%
  select(-tree) %>%
  pivot_longer(-tree_id, names_to = "topology",values_to="weight") %>%
  group_by(topology) %>%
  summarize(weight_count=sum(weight)/32)  %>%
  arrange(desc(topology)) %>%
  ggplot(.,aes(x=fct_relevel(topology, c("pet_ano","ann_ano","ann_pet")),y=weight_count,fill=topology)) +
  geom_bar(stat="identity") +
  theme_cowplot() +
  xlab("") +
  scale_fill_manual(values=c("#a8dadc","#e63946","#1d3557")) +
  theme(legend.position = "none") +
  ylab("Topology counts")
ggsave("plots/ann.petfal.ano.weights.twisst.v1.pdf",width=12,height=4)

twisst %>%
  cbind(tree_ids) %>%
  select(-tree) %>%
  as_tibble() %>%
  separate(tree_id,c("chr","start","end"),"-",convert=T) %>%
  arrange(chr,start) %>%
  pivot_longer(-c(chr,start,end), names_to = "topology",values_to="weight") %>%
  mutate(percent = weight/32) -> twisst_long

twisst_long %>%
  group_by(chr,topology) %>%
  mutate(roll_mean_percent = roll_mean(percent, n = 5, align = "right", fill = NA),
         roll_mean_pos = roll_mean(start, n = 5, align = "right", fill = NA)) %>%
  filter(chr == "Ha412HOChr03") %>% 
  ggplot(.,aes(x=roll_mean_pos, y=roll_mean_percent,color=topology,fill=topology)) +
  geom_area()

###DCT plot
dct_results <- read_tsv("gowens22.mappable.dp4.missing80.10.2kb.physical.subsample.dct.results.txt")
dct_results %>%
  filter(P1 == "ano",P2 == "pet_fal",P3 == "ann") %>%
  select(P1, `12`,`13`,`23`) %>%
  pivot_longer(-P1, names_to = "topology",values_to="count") %>%
  ggplot(.,aes(x=topology,y=count,fill=topology)) +
  geom_bar(stat="identity") +
  theme_cowplot() +
  xlab("") +
  scale_fill_manual(values=c("#e63946","#a8dadc","#1d3557")) +
  theme(legend.position = "none")
  
#BLT TEST RESULTS
test_tibble <- read_tsv("gowens22.mappable.dp4.missing80.10.2kb.physical.subsample.topology.txt")

test_tibble %>%
  filter(topology != "12") %>%
  mutate(dtype = topology, d = case_when(topology == 23 ~ dist_23,
                                         topology == 13 ~ dist_13)) %>%
  mutate(d = as.numeric(d)) -> test_blt

test_blt %>%
  filter(P1 == "ano",P2 == "pet_fal",P3=="ann") %>%
  select(P1,dist_13, dist_23) %>%
  pivot_longer(-P1, names_to = "type",values_to= "distance") %>%
  ggplot(.,aes(x=type,y=distance,fill=type)) +
  geom_boxplot(alpha=0.7) +
  theme_cowplot() +
  scale_fill_manual(values=c("#a8dadc","#e63946"),
                    name="Topology") +
  ylab("Genetic\ndistance") +
  xlab("") +
  coord_cartesian(ylim=c(0,0.08))
ggsave("plots/gowens22.mappable.dp4.missing80.10.2kb.physical.subsample.ano.blt.pdf")
