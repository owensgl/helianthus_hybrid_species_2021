library(tidyverse)
library(cowplot)

data <- read_tsv("../vcf/gowens22.mappable.variant.snps.stats.subset.txt",
         col_names=c("QD","DP","MQ","MQRankSum","FS","ReadPosSum","SOR"),na=".")

limits <- tibble(stat=c("DP","FS","QD","SOR","MQ"),
                 cutoff=c(1000,50,2,4,30 ))
pdf("plots/vcf_quality_cutoffs.v1.pdf")
data %>%
  mutate(locus=row_number()) %>%
  pivot_longer(-locus, names_to = "stat", values_to = "value") %>%
  ggplot(.,aes(value)) +
  geom_density() +
  facet_wrap(~stat,scales="free") +
  theme_cowplot() +
  geom_vline(data=limits,aes(xintercept=cutoff),
             linetype="dotted")
dev.off()
