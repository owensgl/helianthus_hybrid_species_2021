library(tidyverse)
library(ape)
library(rstatix)
library(cowplot)
tree_pattern <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10kb.geneconcord.cf.stat_tree",comment = "#")
tree_names <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10kb.treenames.trees.txt",col_names = c("loci","tree")) %>%
  mutate(TreeID = row_number())

tree_size <- read_tsv("../vcf/gowens22.mappable.dp4.missing80.10kb.windowsize.txt",
                      col_names = c("chr","start","end","length"))
tree_names <- tree_names %>% 
  mutate(size = tree_size$length)
tree_names <- tree_names %>%
  mutate(size_quantile_rank = ntile(tree_names$size,5))

genetic_map <- read_tsv("/home/owens/ref/Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt") %>% 
  mutate(rate = lead(cM)-cM) %>% 
  rename(start = pos) %>% mutate(end = start + 999999) %>%
  mutate(rate = case_when(is.na(rate) ~ cM - lag(cM),
                          rate < 0 ~ cM - lag(cM),
                          TRUE ~ rate))
tree_names %>%
  separate(loci, c("chr","start","end"),"-",convert=T) %>%
  select(chr,start, end) %>% unique() -> tree_windows

tree_windows$cm_rate <- NA
for (i in 1:nrow(tree_windows)){
  print(i)
  chosen_start <- tree_windows$start[i]
  chosen_end <-tree_windows$end[i]
  chosen_chr <- tree_windows$chr[i]
  start_rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
    filter(start <= chosen_start, end >= chosen_start) %>%
    pull(rate)
  end_rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
    filter(start <= chosen_end, end >= chosen_end) %>%
    pull(rate)
  mean_rate <- mean(start_rate,end_rate)
  tree_windows$cm_rate[i] <- mean_rate
}
tree_windows %>%
  mutate(quantile_rank = ntile(tree_windows$cm_rate,5)) %>%
  mutate(loci = paste(chr,start,end,sep="-")) %>%
  inner_join(tree_names)-> tree_windows



tree_pattern %>%
  filter(ID == "25") %>%
  inner_join(tree_names) %>%
  group_by(ID) %>%
  summarize(p_gC = sum(gC)/n(),
            p_gD1 = sum(gD1)/n(),
            p_gD2 = sum(gD2)/n()) %>%
  dplyr::select(ID, p_gC,p_gD1,p_gD2) %>%
  pivot_longer(-ID, names_to = "proportion",values_to="value") %>%
  ggplot(.,aes(x=ID,y=value,fill=proportion)) +
  geom_bar(stat="identity",position="dodge")

tree_pattern %>%
  filter(ID == "25") %>%
  inner_join(tree_names) %>%
  inner_join(tree_windows %>% mutate(loci = paste(chr,start,end,sep="-"))) %>%
  group_by(ID,quantile_rank) %>%
  summarize(p_gC = sum(gC)/n(),
            p_gD1 = sum(gD1)/n(),
            p_gD2 = sum(gD2)/n()) %>%
  dplyr::select(ID,quantile_rank, p_gC,p_gD1,p_gD2) %>%
  pivot_longer(-c(ID,quantile_rank), names_to = "proportion",values_to="value") %>%
  ggplot(.,aes(x=quantile_rank,y=value,fill=proportion)) +
  geom_bar(stat="identity",position="dodge")



tree_pattern %>%
  inner_join(tree_names) %>%
  inner_join(tree_windows %>% mutate(loci = paste(chr,start,end,sep="-"))) %>%
  group_by(ID) %>%
  mutate(total_gC = sum(gC)/n()) %>%
  filter(total_gC < 0.90) %>%
  #filter(ID == "31") %>%
  group_by(ID,quantile_rank) %>%
  summarize(p_gC = sum(gC)/n(),
            p_gD1 = sum(gD1)/n(),
            p_gD2 = sum(gD2)/n(),
            remainder = 1-p_gC - p_gD1 - p_gD2) %>%
  pivot_longer(-c(ID,quantile_rank), names_to = "type",values_to="value") %>%
  filter(type == "p_gC") %>%
  ggplot(.,aes(x=quantile_rank,y=value,color=as.factor(ID),group=ID)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  #facet_grid(~ID,scales="free_y") +
  xlab("cM quantile") +
  ylab("Window proportion") +
  theme_cowplot()

tree_pattern %>%
  inner_join(tree_names) %>%
  inner_join(tree_windows %>% mutate(loci = paste(chr,start,end,sep="-"))) %>%
  group_by(ID) %>%
  mutate(total_gC = sum(gC)/n()) %>%
  filter(total_gC < 0.90) %>%
  #filter(ID == "31") %>%
  group_by(ID,size_quantile_rank) %>%
  summarize(p_gC = sum(gC)/n(),
            p_gD1 = sum(gD1)/n(),
            p_gD2 = sum(gD2)/n(),
            remainder = 1-p_gC - p_gD1 - p_gD2) %>%
  pivot_longer(-c(ID,size_quantile_rank), names_to = "type",values_to="value") %>%
  filter(type == "p_gC") %>%
  ggplot(.,aes(x=size_quantile_rank,y=value,color=as.factor(ID),group=ID)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  #facet_grid(type~ID,scales="free_y") +
  xlab("Information quantile") +
  ylab("Window proportion") +
  theme_cowplot()

tree_names %>%
  inner_join(tree_windows %>% mutate(loci = paste(chr,start,end,sep="-"))) %>%
  ggplot(.,aes(x=size,y=cm_rate)) +
  geom_point() +
  geom_smooth(method="lm") +
  scale_x_log10() +
  scale_y_log10()

cor.test(tree_names %>%
  inner_join(tree_windows %>% mutate(loci = paste(chr,start,end,sep="-"))) %>%
    select(cm_rate) %>% pull(),
  tree_names %>%
    inner_join(tree_windows %>% mutate(loci = paste(chr,start,end,sep="-"))) %>%
    select(size)  %>% pull(),
  method="spearman"
)
tree_names %>%
  inner_join(tree_windows %>% mutate(loci = paste(chr,start,end,sep="-"))) %>%
  group_by(size_quantile_rank, quantile_rank) %>%
  summarize(count=n()) %>%
  ggplot(.,aes(x=size_quantile_rank,y=quantile_rank)) +
  geom_tile(aes(fill=count)) +
  scale_fill_viridis_c()

tree_pattern %>%
  inner_join(tree_names) %>%
  inner_join(tree_windows %>% mutate(loci = paste(chr,start,end,sep="-"))) %>%
  group_by(ID) %>%
  mutate(total_gC = sum(gC)/n()) %>%
  filter(total_gC < 0.90) %>%
  #filter(ID == "31") %>%
  group_by(ID,size_quantile_rank,quantile_rank) %>%
  summarize(p_gC = sum(gC)/n()) %>%
  pivot_longer(-c(ID,size_quantile_rank,quantile_rank), names_to = "type",values_to="value") %>%
  ggplot(.,aes(x=quantile_rank,y=value,color=as.factor(ID),group=ID)) +
  geom_point() +
  geom_smooth(method="lm",) +
  facet_grid(size_quantile_rank~ID,scales="free_y") +
  xlab("cM quantile") +
  ylab("Window proportion") +
  theme_cowplot()

#Trying logistic regression
library(car)
tree_pattern %>%
  inner_join(tree_names) %>%
  inner_join(tree_windows %>% mutate(loci = paste(chr,start,end,sep="-"))) %>%
  group_by(ID) %>%
  mutate(total_gC = sum(gC)/n()) %>%
  filter(total_gC < 0.90) %>%
  select(ID) %>% unique() -> testable_nodes


glm.fit <- tree_pattern %>%
  inner_join(tree_names) %>%
  inner_join(tree_windows %>% mutate(loci = paste(chr,start,end,sep="-"))) %>%
  group_by(ID) %>%
  mutate(total_gC = sum(gC)/n()) %>%
  filter(total_gC < 0.90) %>%
  mutate(ID = as.character(ID)) %>%
  select(cm_rate,gC,ID) %>%
  glm(gC ~ ID * cm_rate , data = ., family = binomial) 

y <- summary(glm.fit)
z <- anova(glm.fit,test="LRT")

glm.fit <- tree_pattern %>%
  inner_join(tree_names) %>%
  inner_join(tree_windows %>% mutate(loci = paste(chr,start,end,sep="-"))) %>%
  group_by(ID) %>%
  mutate(total_gC = sum(gC)/n()) %>%
  filter(total_gC < 0.90) %>%
  mutate(ID = as.character(ID)) %>%
  select(cm_rate,gC,ID,size) %>%
  glm(gC ~ ID * size , data = ., family = binomial) 

y <- summary(glm.fit)
z <- anova(glm.fit,test="LRT")


results_glm <- tibble()
for (node in testable_nodes$ID){
  

glm.fit <- tree_pattern %>%
  inner_join(tree_names) %>%
  inner_join(tree_windows %>% mutate(loci = paste(chr,start,end,sep="-"))) %>%
  group_by(ID) %>%
  mutate(total_gC = sum(gC)/n()) %>%
  filter(total_gC < 0.90) %>%
  filter(ID == node) %>%
  select(cm_rate,gC,size) %>%
  glm(gC ~ cm_rate, data = ., family = binomial) 

y <- summary(glm.fit)
est <- y$coefficients[2,1]
p_1 <- y$coefficients[2,4]

glm.fit <- tree_pattern %>%
  inner_join(tree_names) %>%
  inner_join(tree_windows %>% mutate(loci = paste(chr,start,end,sep="-"))) %>%
  group_by(ID) %>%
  mutate(total_gC = sum(gC)/n()) %>%
  filter(total_gC < 0.90) %>%
  filter(ID == node) %>%
  select(cm_rate,gC,size) %>%
  glm(gC ~ size + cm_rate, data = ., family = binomial) 

x <- summary(glm.fit)
z <- anova(glm.fit,test="LRT")
p_2 <- z$`Pr(>Chi)`[3]
results_glm <- rbind(results_glm, tibble(node=node, est=est,p_single=p_1,p_seq=p_2))
}
write_tsv(results_glm, "gowens22.mappable.dp4.missing80.10kb.geneconcord.cf.glmspeciestree.txt" )
