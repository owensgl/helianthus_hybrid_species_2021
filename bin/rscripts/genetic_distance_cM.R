library(tidyverse)
library(rstatix)

dist_files <- list.files("../vcf/seqdist/")
species_list <- read_tsv("species_id.txt")
all_distances <- tibble()
for (i in dist_files){
  tmp <- read_tsv(paste0("../vcf/seqdist/",i)) %>%
    select(Header_1,Header_2,Distance) %>%
    mutate(region = gsub(".seqdist","",i)) %>%
    separate(region, c("chr","start","end"),"-",convert=T)
  tmp2 <- tmp %>%
    rename(Header_2 = Header_1, Header_1 = Header_2) 
  
  all_distances <- rbind(all_distances,tmp)
  all_distances <- rbind(all_distances,tmp2)
}
write_tsv(all_distances, "all_seqdist.txt.gz")
all_distances <- read_tsv("all_seqdist.txt.gz")

perennials <- c("664647_GIG","DIV_1956","DEC_1895","GRO_2043")

genetic_map <- read_tsv("/home/owens/ref/Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt") %>% 
  mutate(rate = lead(cM)-cM) %>% 
  rename(start = pos) %>% mutate(end = start + 999999) %>%
  mutate(rate = case_when(is.na(rate) ~ cM - lag(cM),
                          rate < 0 ~ cM - lag(cM),
                          TRUE ~ rate))

outgroup_dist <- all_distances %>%
  filter(Header_2 %in% perennials,
         !Header_1 %in% perennials) %>%
  inner_join(species_list %>% rename(Header_1 = sample)) %>%
  group_by(species,chr,start,end,) %>%
  summarize(mean_dist = mean(Distance)) 
outgroup_dist %>%
  group_by(chr,start,end) %>%
  summarize() -> windows
  
windows$cm_rate <- NA
for (i in 1:nrow(windows)){
  print(i)
  chosen_start <- windows$start[i]
  chosen_end <-windows$end[i]
  chosen_chr <- windows$chr[i]
  start_rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
    filter(start <= chosen_start, end >= chosen_start) %>%
    pull(rate)
  end_rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
    filter(start <= chosen_end, end >= chosen_end) %>%
    pull(rate)
  mean_rate <- mean(start_rate,end_rate)
  windows$cm_rate[i] <- mean_rate
}
windows %>%
  ungroup() %>%
  mutate(quantile_rank = ntile(windows$cm_rate,5)) -> windows
outgroup_dist %>%
  inner_join(windows) %>%
  ggplot(.,aes(x=as.factor(quantile_rank),y=mean_dist,fill=species)) +
  geom_boxplot() +
  facet_wrap(~species) +
  ylab("Distance from perennials") +
  xlab("cM quintile")

outgroup_dist %>%
  inner_join(windows) %>%
  ggplot(aes(x=cm_rate,y=mean_dist,color=species)) +
  geom_smooth(method="lm") +
  theme_cowplot() +
  ylab("Distance from perennials")
  
all_distances %>%
  filter(Header_2 %in% perennials,
         !Header_1 %in% perennials) %>%
  inner_join(species_list %>% rename(Header_1 = sample)) %>%
  group_by(species,chr,start,end,) %>%
  summarize(mean_dist = mean(Distance)) %>%
  filter(chr == "Ha412HOChr15") %>%
  ggplot(.,aes(x=start,y=mean_dist,color=species)) +
  geom_line() +
  theme_cowplot()
  
  ggplot(.,aes(x=start,y=Distance,color=Header_1)) + 
  geom_line()
all_distances %>%
  inner_join(species_list %>% rename(Header_1 = sample, species_1 = species)) %>%
  filter(species_1 == "ano") %>%
  rename(sample=Header_2) %>%
  inner_join(species_list) %>%
  group_by(species,start,end) %>%
  summarize(mean_dist = mean(Distance)) %>%
  filter(species == "ann" | species == "pet_fal") %>%
  pivot_wider(names_from = species,values_from=mean_dist) %>%
  mutate(dist_dif = ann - pet_fal,
         type=case_when(dist_dif < 0 ~ "ann_like",
                        dist_dif > 0 ~ "pet_like"),
         constant=1) %>%
  select(start,dist_dif,type) %>%
  pivot_wider(values_from = dist_dif,names_from=type,values_fill=0) %>%
  mutate(         ymax_ann = case_when(ann_like > 0 ~ ann_like,
                                   TRUE ~ 0),
                  ymin_ann = case_when(ann_like < 0 ~ ann_like,
                                   TRUE ~ 0),
                  ymax_pet = case_when(pet_like > 0 ~ pet_like,
                                       TRUE ~ 0),
                  ymin_pet = case_when(pet_like < 0 ~ pet_like,
                                       TRUE ~ 0)) %>%
ggplot(.,aes(x=start)) +
  geom_ribbon(aes(ymax = ymax_ann, ymin = ymin_ann), colour = NA,fill="red") +
  geom_ribbon(aes(ymax = ymax_pet, ymin = ymin_pet), colour = NA,fill="blue") +
  geom_hline(yintercept=0,linetype="dotted") +
  theme_cowplot()
####Distance to root
pers <- c("664647_GIG","DEC_1895","DIV_1956","GRO_2043")
species_id <- read_tsv("../meta/species_ids.txt") %>%
  rename(Header_2 = sample)

all_distances %>%
  filter(Header_1 %in% pers) %>%
  filter(!Header_2 %in% pers) %>%
  inner_join(species_id) %>%
  group_by(species, chr, start, end) %>%
  summarize(outgroup_dist = mean(Distance)) -> outgroup_dist
  
all_distances %>%
  filter(Header_1 == "ANN1029" | Header_1 == "ANN1283") %>%
  inner_join(species_id) %>%
  filter(species != "ann") %>%
  inner_join(outgroup_dist) %>%
  mutate(normalize_dist = Distance/outgroup_dist) %>%
  filter(normalize_dist < 3) %>%
  ggplot(.,aes(x=species,y=normalize_dist,color=species)) +
  geom_jitter(width=0.2,alpha=0.2) +
  geom_boxplot(outlier.shape=NA) +
  geom_hline(aes(yintercept=1),linetype="dotted") +
  theme_cowplot() +
  ggtitle("Genetic distance from annuus, normalized against perennials")
  
all_distances %>%
  filter(Header_1 == "ANN1029" | Header_1 == "ANN1283") %>%
  inner_join(species_id) %>%
  filter(species != "ann") %>%
  inner_join(outgroup_dist) %>%
  mutate(normalize_dist = Distance/outgroup_dist) %>%
  group_by(species) %>%
  summarize(mean_dist = mean(normalize_dist,na.rm=T),
            sd_dist = sd(normalize_dist,na.rm=T)) %>%
  ggplot(.,aes(x=fct_reorder(species,mean_dist),color=species)) +
  geom_point(aes(y=mean_dist)) + 
  geom_errorbar(aes(ymin=mean_dist-sd_dist,ymax=mean_dist+sd_dist)) +
  theme_cowplot() +
  geom_hline(aes(yintercept=1),linetype="dotted") +
  ggtitle("Genetic distance from annuus, normalized against perennials")
  

