library(tidyverse)
library(cowplot)
dstat <- read_tsv("../d/gowens22.mappable.variant.snps.filtered.dp4.missing80.vcf.gz_BBAA.txt")

dstat$pvalue_hoch <- p.adjust(dstat$`p-value`, method="hochberg")
introgression_summary <- tibble()
for (i in 1:nrow(dstat)){
  print(i)
  P1 <- dstat$P1[i]
  P2 <- dstat$P2[i]
  P3 <- dstat$P3[i]
  p <- dstat$pvalue_hoch[i]
  if (p < 0.05){
    tmp <- tibble(spe1 = P2, spe2 = P3, signal=1)
    introgression_summary <- rbind(introgression_summary, tmp)
    tmp <- tibble(spe1 = P3, spe2 = P2, signal=1)
    introgression_summary <- rbind(introgression_summary, tmp)
  }else{
    tmp <- tibble(spe1 = P2, spe2 = P3, signal=0)
    introgression_summary <- rbind(introgression_summary, tmp)
    tmp <- tibble(spe1 = P3, spe2 = P2, signal=0)
    introgression_summary <- rbind(introgression_summary, tmp)
  }
  tmp <- tibble(spe1 = P1, spe2 = P3, signal=0)
  introgression_summary <- rbind(introgression_summary, tmp)
  tmp <- tibble(spe1 = P3, spe2 = P1, signal=0)
  introgression_summary <- rbind(introgression_summary, tmp)
}
species_order <- tibble(spe = c("par","arg","ann","ano","des","niv","deb","pet_fal","pet_pet"),
                        rank = 1:9)
introgression_summary %>%
  group_by(spe1,spe2) %>%
  summarize(introgression_signal = mean(signal), tests=n(), proportion=paste0(tests*introgression_signal,"/", tests)) %>%
  inner_join(species_order %>% rename(spe1 = spe,rank1=rank)) %>%
  inner_join(species_order %>% rename(spe2 = spe,rank2=rank)) %>%
  filter(rank1 > rank2) %>%
  ggplot(.,aes(x=fct_reorder(spe1,-rank1),y=fct_reorder(spe2,-rank2),fill=introgression_signal)) +
  geom_tile(color="black",size=2) +
  scale_fill_viridis_c(name="Dstat Introgression\nsignal") +
  theme_cowplot() +
  xlab("") +
  ylab("") +
  geom_label(aes(label=proportion),fill="white")
  

dstat %>%
  mutate(sig = case_when(pvalue_hoch < 0.05 ~ "sig",
                         TRUE ~ "non-sig")) %>%
  group_by(sig) %>%
  summarize(n=n()) %>%
  mutate(freq = n / sum(n)) %>%
  ggplot(.,aes(x="Dstat",y=freq,fill=sig)) +
  geom_bar(stat="identity") +
  theme_cowplot() +
  scale_fill_brewer(palette = "Set1",name="Significant?") +
  xlab("") + ylab("Proportion of trios")
  
