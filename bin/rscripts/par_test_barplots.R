library(tidyverse)
library(RcppRoll)
data <- read_tsv("../d/gowens22.mappable.variant.snps.filtered.dp4.missing80_BBAA.txt")
data %>%
  filter(P1 == "ann",P2 == "par", P3 == "pet_fal")  %>%
  select(P1, BBAA, ABBA,BABA) %>%
  pivot_longer(-P1, names_to = "type",values_to="count") %>%
  arrange(desc(type)) %>%
  ggplot(.,aes(x=fct_relevel(type,c("ABBA","BBAA","BABA")),y=count,fill=type)) +
  geom_bar(stat="identity") +
  theme_cowplot() +
  xlab("") +
  scale_fill_manual(values=c("#1d3557","#e63946","#a8dadc")) +
  theme(legend.position = "none") +
  ylab("Site count")
ggsave("plots/gowens22.mappable.variant.snps.filtered.dp4.missing80_BBAA.par.v1.pdf",width=12,height=4)

twisst <- read_tsv("../twisst/ann.petfal.par.weights.txt",comment="#",
                   col_names=c("ann_pet","ann_par","pet_par"),skip=4)
tree_ids <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10.2kb.physical.treenames.trees.txt",
                     col_names=c("tree_id","tree"))

chisq.test(c(1054,1039))$p.value
twisst %>%
  cbind(tree_ids) %>%
  select(-tree) %>%
  pivot_longer(-tree_id, names_to = "topology",values_to="weight") %>%
  group_by(topology) %>%
  summarize(weight_count=sum(weight)/32)  %>%
  arrange(desc(topology)) %>%
  ggplot(.,aes(x=fct_relevel(topology, c("pet_par","ann_par","ann_pet")),y=weight_count,fill=topology)) +
  geom_bar(stat="identity") +
  theme_cowplot() +
  xlab("") +
  scale_fill_manual(values=c("#a8dadc","#e63946","#1d3557")) +
  theme(legend.position = "none") +
  ylab("Topology counts")
ggsave("plots/ann.petfal.par.weights.twisst.v1.pdf",width=12,height=4)

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
  filter(chr == "Ha412HOChr15") %>% 
  ggplot(.,aes(x=roll_mean_pos, y=roll_mean_percent,color=topology,fill=topology)) +
  geom_area()


#BLT TEST RESULTS
test_tibble <- read_tsv("gowens22.mappable.dp4.missing80.10.2kb.physical.subsample.topology.txt")

test_tibble %>%
  filter(topology != "12") %>%
  mutate(dtype = topology, d = case_when(topology == 23 ~ dist_23,
                                         topology == 13 ~ dist_13)) %>%
  mutate(d = as.numeric(d)) -> test_blt

test_blt %>%
  filter(P1 == "ann",P2 == "par",P3=="pet_fal") %>%
  select(P1,dist_13, dist_23) %>%
  pivot_longer(-P1, names_to = "type",values_to= "distance") %>%
  ggplot(.,aes(x=type,y=distance,fill=type)) +
  geom_boxplot(alpha=0.7) +
  theme_cowplot() +
  scale_fill_manual(values=c("#e63946","#1d3557"),
                    name="Topology") +
  ylab("Genetic\ndistance") +
  xlab("")
ggsave("plots/gowens22.mappable.dp4.missing80.10.2kb.physical.subsample.par.blt.pdf")


#quibl results
quibl <- read_tsv("../quibl/gowens22.mappable.dp4.missing80.10kb.subsample.quibl.formatted.txt")

quibl$triplet <- mgsub(quibl$triplet , c("NHKG_16_01_R_196", "GO8_GB179","PAR_posas_01","ANN1283","ARG0295","DEB_1135","PET0662","PET0424","PET0495","664647_GIG"), 
                       c("ano", "des","par","ann","arg","deb","niv", "petfal","petpet","outgroup"))

quibl$outgroup <- mgsub(quibl$outgroup , c("NHKG_16_01_R_196", "GO8_GB179","PAR_posas_01","ANN1283","ARG0295","DEB_1135","PET0662","PET0424","PET0495","664647_GIG"), 
                        c("ano", "des","par","ann","arg","deb","niv", "petfal","petpet","outgroup"))
quibl %>% 
  filter(triplet == "ann_par_petfal")
