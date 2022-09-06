library(tidyverse)
library(PNWColors)
library(patchwork)
tree_ids <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10.2kb.physical.treenames.trees.txt",
                     col_names=c("tree_id","tree"))

twisst_ano <- read_tsv("../twisst/anngroup.petgroup.ano.weights.txt",comment="#",
                       col_names=c("ann_pet","ann_hyb","pet_hyb"),skip=4)
twisst_ano %>%
  cbind(tree_ids) %>%
  select(-tree) %>%
  as_tibble() %>%
  separate(tree_id,c("chr","start","end"),"-",convert=T) %>%
  arrange(chr,start) %>%
  pivot_longer(-c(chr,start,end), names_to = "topology",values_to="weight") %>%
  mutate(percent = weight/256) -> twisst_ano_long

twisst_des <- read_tsv("../twisst/anngroup.petgroup.des.weights.txt",comment="#",
                       col_names=c("ann_pet","ann_hyb","pet_hyb"),skip=4)
twisst_des %>%
  cbind(tree_ids) %>%
  select(-tree) %>%
  as_tibble() %>%
  separate(tree_id,c("chr","start","end"),"-",convert=T) %>%
  arrange(chr,start) %>%
  pivot_longer(-c(chr,start,end), names_to = "topology",values_to="weight") %>%
  mutate(percent = weight/256) -> twisst_des_long

twisst_par <- read_tsv("../twisst/anngroup.petgroup.par.weights.txt",comment="#",
                       col_names=c("ann_pet","ann_hyb","pet_hyb"),skip=4)
twisst_par %>%
  cbind(tree_ids) %>%
  select(-tree) %>%
  as_tibble() %>%
  separate(tree_id,c("chr","start","end"),"-",convert=T) %>%
  arrange(chr,start) %>%
  pivot_longer(-c(chr,start,end), names_to = "topology",values_to="weight") %>%
  mutate(percent = weight/256) -> twisst_par_long


twisst_ano_plot <- twisst_ano_long %>%
  group_by(chr,topology) %>%
  mutate(next_end = roll_max(start, n = 2, align = "left", fill = NA)) %>%
  group_by(chr, start,next_end) %>%
  summarize(max_per = max(percent), type = topology[which(percent == max_per)]) %>%
  mutate(type = case_when(max_per < 0.5 ~ "unknown",
                          TRUE ~ type)) %>%
  mutate(chr = gsub("Ha412HO","",chr)) %>%
  mutate(species = "Ano")

twisst_des_plot <- twisst_des_long %>%
  group_by(chr,topology) %>%
  mutate(next_end = roll_max(start, n = 2, align = "left", fill = NA)) %>%
  group_by(chr, start,next_end) %>%
  summarize(max_per = max(percent), type = topology[which(percent == max_per)]) %>%
  mutate(type = case_when(max_per < 0.5 ~ "unknown",
                          TRUE ~ type)) %>%
  mutate(chr = gsub("Ha412HO","",chr)) %>%
  mutate(species = "Des")

twisst_par_plot <- twisst_par_long %>%
  group_by(chr,topology) %>%
  mutate(next_end = roll_max(start, n = 2, align = "left", fill = NA)) %>%
  group_by(chr, start,next_end) %>%
  summarize(max_per = max(percent), type = topology[which(percent == max_per)]) %>%
  mutate(type = case_when(max_per < 0.5 ~ "unknown",
                          TRUE ~ type)) %>%
  mutate(chr = gsub("Ha412HO","",chr)) %>%
  mutate(species = "Par")

#pdf("plots/twisst_example_smoothed.pdf",height=12,width=24)
twisst_long %>%
  group_by(chr,topology) %>%
  mutate(roll_mean_percent = roll_mean(percent, n = 2, align = "right", fill = NA),
         roll_mean_pos = roll_mean(start, n = 2, align = "right", fill = NA)) %>%
  mutate(chr = gsub("Ha412HO","",chr)) %>%
  ggplot(.,aes(x=roll_mean_pos/1000000, y=roll_mean_percent,color=topology,fill=topology)) +
  geom_area() +
  theme_cowplot() +
  facet_wrap(~chr) +
  xlab("MBp") 

pdf("plots/anngroup.petgroup.hybrid.twisstmajority.v1.pdf",height=4,width=12)
colors <- pnw_palette("Bay",3)
rbind(twisst_ano_plot, twisst_des_plot) %>%
  rbind(., twisst_par_plot) %>%
  ggplot(.,aes(x=start/1000000, xend=next_end/1000000, y=species,yend=species,color=type)) +
  geom_segment(size=3) +
  theme_cowplot() +
  facet_wrap(~chr) +
  xlab("MBp") + 
  scale_color_manual(values=c(colors,"grey"),
                     name="Dominant topology",
                     labels=c("ANN-HYB",
                              "ANN-PET",
                              "PET-HYB",
                              "UNCERTAIN")) +
  ylab("Species") 
dev.off()
### Correlation between measures using ANN-HYB

twisst_par_annhyb <- twisst_par_long %>%
  mutate(species = "par") %>%
  filter(topology == "ann_hyb") 

twisst_ano_annhyb <- twisst_ano_long %>%
  mutate(species = "ano") %>%
  filter(topology == "ann_hyb") 

twisst_des_annhyb <- twisst_des_long %>%
  mutate(species = "des") %>%
  filter(topology == "ann_hyb") 

tmp1 <- rbind(twisst_ano_annhyb,twisst_des_annhyb, twisst_par_annhyb ) %>%
  select(chr,start, percent,species) %>%
  pivot_wider(names_from = species,values_from=percent)
des_par_cor <- cor(tmp1$des, tmp1$par,method="pearson") 

ano_des_cor <- cor(tmp1$des, tmp1$ano,method="pearson") 

ano_par_cor <- cor(tmp1$par, tmp1$ano,method="pearson") 

plot_1 <- rbind(twisst_ano_annhyb,twisst_des_annhyb, twisst_par_annhyb ) %>%
  select(chr,start, percent,species) %>%
  pivot_wider(names_from = species,values_from=percent) %>%
  ggplot(.,aes(x=par,y=des)) + geom_density_2d_filled(alpha = 0.5) + geom_smooth(method="lm") +
  theme(legend.position = "none") +
  ylab("Des") +
  xlab("Par") +
  ggtitle("Proportion ANN-HYB, r2 = 0.04")

rbind(twisst_ano_annhyb,twisst_des_annhyb, twisst_par_annhyb ) %>%
  select(chr,start, percent,species) %>%
  pivot_wider(names_from = species,values_from=percent) %>%
  lm(ano ~ par,data=.) %>% summary()

plot_2 <- rbind(twisst_ano_annhyb,twisst_des_annhyb, twisst_par_annhyb ) %>%
  select(chr,start, percent,species) %>%
  pivot_wider(names_from = species,values_from=percent) %>%
  ggplot(.,aes(x=par,y=ano)) + geom_density_2d_filled(alpha = 0.5) + geom_smooth(method="lm") +
  theme(legend.position = "none") +
  ylab("Ano") +
  xlab("Par") +
  ggtitle("Proportion ANN-HYB, r2 = 0.04")

rbind(twisst_ano_annhyb,twisst_des_annhyb, twisst_par_annhyb ) %>%
  select(chr,start, percent,species) %>%
  pivot_wider(names_from = species,values_from=percent) %>%
  lm(ano ~ des,data=.) %>% summary()
plot_3 <- rbind(twisst_ano_annhyb,twisst_des_annhyb, twisst_par_annhyb ) %>%
  select(chr,start, percent,species) %>%
  pivot_wider(names_from = species,values_from=percent) %>%
  ggplot(.,aes(x=des,y=ano)) + geom_density_2d_filled(alpha = 0.5) + geom_smooth(method="lm") +
  theme(legend.position = "none") +
  ylab("Ano") +
  xlab("Des") +
  ggtitle("Proportion ANN-HYB, r2 = 0.49")

pdf("plots/anngroup.petgroup.hybrid.twisstcor.v1.pdf",height=4,width=12)
plot_1 + plot_2 + plot_3
dev.off()


### Correlation between measures using PET-ann

twisst_par_annhyb <- twisst_par_long %>%
  mutate(species = "par") %>%
  filter(topology == "pet_hyb") 

twisst_ano_annhyb <- twisst_ano_long %>%
  mutate(species = "ano") %>%
  filter(topology == "pet_hyb") 

twisst_des_annhyb <- twisst_des_long %>%
  mutate(species = "des") %>%
  filter(topology == "pet_hyb") 

tmp2 <- rbind(twisst_par_annhyb,twisst_ano_annhyb, twisst_des_annhyb ) %>%
  select(chr,start, percent,species) %>%
  pivot_wider(names_from = species,values_from=percent)
des_par_cor <- cor(tmp2$des, tmp2$par,method="pearson") 

ano_des_cor <- cor(tmp2$des, tmp2$ano,method="pearson") 

ano_par_cor <- cor(tmp2$par, tmp2$ano,method="pearson") 

plot_1 <- rbind(twisst_ano_annhyb,twisst_des_annhyb, twisst_par_annhyb ) %>%
  select(chr,start, percent,species) %>%
  pivot_wider(names_from = species,values_from=percent) %>%
  ggplot(.,aes(x=par,y=des)) + geom_density_2d_filled(alpha = 0.5) + geom_smooth(method="lm") +
  theme(legend.position = "none") +
  ylab("Des") +
  xlab("Par") +
  ggtitle("Proportion ANN-PET, r2 = 0.04")

rbind(twisst_ano_annhyb,twisst_des_annhyb, twisst_par_annhyb ) %>%
  select(chr,start, percent,species) %>%
  pivot_wider(names_from = species,values_from=percent) %>%
  lm(ano ~ par,data=.) %>% summary()

plot_2 <- rbind(twisst_ano_annhyb,twisst_des_annhyb, twisst_par_annhyb ) %>%
  select(chr,start, percent,species) %>%
  pivot_wider(names_from = species,values_from=percent) %>%
  ggplot(.,aes(x=par,y=ano)) + geom_density_2d_filled(alpha = 0.5) + geom_smooth(method="lm") +
  theme(legend.position = "none") +
  ylab("Ano") +
  xlab("Par") +
  ggtitle("Proportion ANN-PET, r2 = 0.04")

rbind(twisst_ano_annhyb,twisst_des_annhyb, twisst_par_annhyb ) %>%
  select(chr,start, percent,species) %>%
  pivot_wider(names_from = species,values_from=percent) %>%
  lm(ano ~ des,data=.) %>% summary()
plot_3 <- rbind(twisst_ano_annhyb,twisst_des_annhyb, twisst_par_annhyb ) %>%
  select(chr,start, percent,species) %>%
  pivot_wider(names_from = species,values_from=percent) %>%
  ggplot(.,aes(x=des,y=ano)) + geom_density_2d_filled(alpha = 0.5) + geom_smooth(method="lm") +
  theme(legend.position = "none") +
  ylab("Ano") +
  xlab("Des") +
  ggtitle("Proportion ANN-PET, r2 = 0.50")

pdf("plots/anngroup.petgroup.hybrid.twisstcor.v3.pdf",height=4,width=12)
plot_1 + plot_2 + plot_3
dev.off()
