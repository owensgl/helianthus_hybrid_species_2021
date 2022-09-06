library(tidyverse)

quibl <- read_tsv("../quibl/gowens22.mappable.dp4.missing80.10kb.subsample.quibl.formatted.txt")

quibl$triplet <- mgsub(quibl$triplet , c("NHKG_16_01_R_196", "GO8_GB179","PAR_posas_01","ANN1283","ARG0295","DEB_1135","PET0662","PET0424","PET0495","664647_GIG"), 
                      c("ano", "des","par","ann","arg","deb","niv", "petfal","petpet","outgroup"))

quibl$outgroup <- mgsub(quibl$outgroup , c("NHKG_16_01_R_196", "GO8_GB179","PAR_posas_01","ANN1283","ARG0295","DEB_1135","PET0662","PET0424","PET0495","664647_GIG"), 
                       c("ano", "des","par","ann","arg","deb","niv", "petfal","petpet","outgroup"))
quibl %>%
  separate(triplet,c("P1","P2","P3"), "_") %>%
  filter(Species_tree == "no") %>%
  mutate(sig = case_when(`BIC2-1` < -30 ~ 1,
                         TRUE ~ 0)) %>%
  group_by(P1, P2, P3) %>%
  summarize(sig = max(sig)) %>%
  group_by(sig) %>%
  summarize(n=n()) %>%
  mutate(freq = n / sum(n)) %>%
  ggplot(.,aes(x="Dstat",y=freq,fill=sig)) +
  geom_bar(stat="identity") +
  theme_cowplot() +
  scale_fill_brewer(palette = "Set1",name="Significant?") +
  xlab("") + ylab("Proportion of trios")

quibl %>%
  filter(grepl("niv",triplet)) %>%
  filter(grepl("des|ano", triplet)) %>%
  filter(outgroup == "niv")
