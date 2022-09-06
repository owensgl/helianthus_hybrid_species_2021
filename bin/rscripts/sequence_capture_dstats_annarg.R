library(tidyverse)
library(PNWColors)
library(patchwork)
colors <- pnw_palette("Bay",5)
dtree <- read_tsv("../sequence_capture/genes.dsuite_tree.txt")

rev_dtree <- dtree %>%
  rename(P1 = P2, P2 = P1) %>%
  mutate(Dstatistic = -Dstatistic)

all_annuals <- c("H_annuus","H_argophyllus", "H_debilis","H_exilis","H_niveus","H_petiolaris","H_praecox")
other_annuals <- c("H_debilis","H_exilis","H_niveus","H_petiolaris","H_praecox")
large_perennials <- c("H_salicifolius","H_maximiliani","H_giganteus",
                      "H_verticillatus","H_grosseserratus","H_nuttalli",
                      "H_divaricatus","H_microcephalus","H_cusickii",
                      "H_arizonensis","H_laciniatus")
south_perennials <- c("H_longifolius","H_carnosus","H_radula",
                      "H_atrorubens","H_silphioides","H_heterophyllus",
                      "H_angustifolius","H_floridanus")
other_perennials <- c("H_mollis","H_occidentalis","H_gracilentus","H_agrestis")
outgroup_perennials <- c("H_porteri")
other_perennials <- c("H_agrestis","H_arizonensis","H_atrorubens","H_carnosus","H_cusickii","H_angustifolius","H_divaricatus","H_giganteus",
                      "H_gracilentus","H_grosseserratus","H_heterophyllus","H_laciniatus","H_longifolius","H_maximiliani","H_microcephalus",
                      "H_mollis","H_nuttalli","H_occidentalis", "H_radula","H_salicifolius","H_silphioides","H_verticillatus","H_floridanus")

plot_1 <- dtree %>%
  rbind(rev_dtree) %>%
  filter(P1 == "H_annuus" | P1 == "H_argophyllus" ) %>%
  filter(P2 != "H_argophyllus" & P2 != "H_annuus") %>%
  filter(P2 %in% other_annuals) %>%
  mutate(P1 = gsub("H_","",P1)) %>%
  mutate(P2 = gsub("H_","",P2)) %>%
  ggplot(.,aes(x=P2,y=Dstatistic,color=P1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(aes(yintercept=0),linetype="dashed") +
  theme_cowplot()  +
  geom_point(position = position_jitterdodge(),
             alpha=0.6) +
  ylab("D") +
  ggtitle("P1-Annual-Perennial") +
  scale_color_manual(values=c("#e63946","#1d3557")) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))


plot_2 <- dtree %>%
  rbind(rev_dtree) %>%
  filter(P1 == "H_maximiliani") %>%
  filter(P2 %in% other_perennials) %>%
  filter(P3 %in% all_annuals) %>%
  mutate(P1 = gsub("H_","",P1)) %>%
  mutate(P2 = gsub("H_","",P2)) %>%
  mutate(P3 = gsub("H_","",P3)) %>%
  ggplot(.,aes(x=P3,y=Dstatistic)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(aes(yintercept=0),linetype="dashed") +
  theme_cowplot() +
  geom_jitter(width=0.2,alpha=0.5) +
  ylab("D") +
  ggtitle("Maximiliani-Perennial-Annual") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

pdf("plots/sequence_caption_candidates.v1.pdf",
    height=8,width=6)
plot_1 / plot_2 
dev.off()

