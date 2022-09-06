library(tidyverse)
library(cowplot)

tree_pattern <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10kb.geneconcord.cf.stat_tree",comment = "#")
gene_concordance <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10kb.geneconcord.cf.stat",comment = "#")

tree_pattern <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10.2kbphysical.concordance.cf.stat_tree",comment = "#")
gene_concordance <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10.2kbphysical.concordance.cf.stat",comment = "#")

chisq = function(DF1, DF2, N){
  tryCatch({
    # converts percentages to counts, runs chisq, gets pvalue
    chisq.test(c(round(DF1*N)/100, round(DF2*N)/100))$p.value
  },
  error = function(err) {
    # errors come if you give chisq two zeros
    # but here we're sure that there's no difference
    return(1.0)
  })
}

gene_concordance %>% 
  group_by(ID) %>%
  mutate(gEF_p = chisq(gDF1, gDF2, gN)) %>% 
  filter(gEF_p < 0.05)
