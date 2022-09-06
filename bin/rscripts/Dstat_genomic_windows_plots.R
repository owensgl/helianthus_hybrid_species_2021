library(tidyverse)
library(cowplot)
library(zoo)
library(patchwork)

# genetic_map <- read_tsv("/media/drive_5_usb/Fuchs/ref/Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt") %>% 
#   mutate(rate = lead(cM)-cM) %>% 
#   rename(start = pos) %>% mutate(end = start + 999999) %>%
#   mutate(rate = case_when(is.na(rate) ~ cM - lag(cM),
#                           rate < 0 ~ cM - lag(cM),
#                           TRUE ~ rate))


# data$cm_rate <- NA
# for (i in 1:nrow(data)){
#   print(i)
#   chosen_start <- data$windowStart[i]
#   chosen_end <-data$windowEnd[i]
#   chosen_chr <- data$chr[i]
#   start_rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
#     filter(start <= chosen_start, end >= chosen_start) %>%
#     pull(rate)
#   end_rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
#     filter(start <= chosen_end, end >= chosen_end) %>%
#     pull(rate)
#   mean_rate <- mean(start_rate,end_rate)
#   data$cm_rate[i] <- mean_rate
# }
pdf("plots/local_Dstat.v2.pdf")
data <- read_tsv("../d/ann_pet_pet_ano_localFstats__50_25.txt")

p1.1 <- data %>%
  ggplot(aes(D)) +
  geom_histogram(fill="#a8dadc") +
  theme_cowplot() +
  theme(axis.text = element_text(size=8)) +
  coord_cartesian(xlim=c(-1,1))


p1.2 <- data %>%
  mutate(chr = gsub("Ha412HO","",chr)) %>%
  mutate(mean_D = rollapply(D, 10, mean, fill = NA, align = "right")) %>%
  mutate(middle = (windowStart + windowEnd)/2) %>%
  ggplot(.,aes(x=middle/1000000,y=D)) +
  geom_point(alpha=0.2,color="#a8dadc") +
  geom_line(aes(x=middle/1000000,y=mean_D)) +
  facet_wrap(~chr) +
  geom_hline(aes(yintercept=0)) +
  scale_color_viridis_d() +
  theme_cowplot() +
  xlab("MBp") +
  ggtitle("Ann-PetPet-Ano") +
  theme(axis.text = element_text(size=8))

p1.2 + inset_element(p1.1, left = 0.4, bottom = 0, right = 1, top = 0.2)


data <- read_tsv("../d/ann_pet_fal_ano_localFstats__50_25.txt")

p2.1 <- data %>%
  ggplot(aes(D)) +
  geom_histogram(fill="#a8dadc") +
  theme_cowplot() +  
  theme(axis.text = element_text(size=8)) +
  coord_cartesian(xlim=c(-1,1))


p2.2 <- data %>%
  mutate(chr = gsub("Ha412HO","",chr)) %>%
  mutate(mean_D = rollapply(D, 10, mean, fill = NA, align = "right")) %>%
  mutate(middle = (windowStart + windowEnd)/2) %>%
  ggplot(.,aes(x=middle/1000000,y=D)) +
  geom_point(alpha=0.2,color="#a8dadc") +
  geom_line(aes(x=middle/1000000,y=mean_D)) +
  facet_wrap(~chr) +
  geom_hline(aes(yintercept=0)) +
  scale_color_viridis_d() +
  theme_cowplot() +
  xlab("MBp") +
  ggtitle("Ann-PetFal-Ano") +
  theme(axis.text = element_text(size=8))



p2.2 + inset_element(p2.1, left = 0.4, bottom = 0, right = 1, top = 0.2)

data <- read_tsv("../d/ann_niv_ano_localFstats__50_25.txt")

p7.1 <- data %>%
  ggplot(aes(D)) +
  geom_histogram(fill="#a8dadc") +
  theme_cowplot() +
  theme(axis.text = element_text(size=8)) +
  coord_cartesian(xlim=c(-1,1))


p7.2 <- data %>%
  mutate(chr = gsub("Ha412HO","",chr)) %>%
  mutate(mean_D = rollapply(D, 10, mean, fill = NA, align = "right")) %>%
  mutate(middle = (windowStart + windowEnd)/2) %>%
  ggplot(.,aes(x=middle/1000000,y=D)) +
  geom_point(alpha=0.2,color="#a8dadc") +
  geom_line(aes(x=middle/1000000,y=mean_D)) +
  facet_wrap(~chr) +
  geom_hline(aes(yintercept=0)) +
  scale_color_viridis_d() +
  theme_cowplot() +
  xlab("MBp") +
  ggtitle("Ann-Niv-Ano") +
  theme(axis.text = element_text(size=8))

p7.2 + inset_element(p1.1, left = 0.4, bottom = 0, right = 1, top = 0.2)

data <- read_tsv("../d/ann_pet_pet_des_localFstats__50_25.txt")

p3.1 <- data %>%
  ggplot(aes(D)) +
  geom_histogram(fill="#a8dadc") +
  theme_cowplot() +
  theme(axis.text = element_text(size=8)) +
  coord_cartesian(xlim=c(-1,1))


p3.2 <- data %>%
  mutate(chr = gsub("Ha412HO","",chr)) %>%
  mutate(mean_D = rollapply(D, 10, mean, fill = NA, align = "right")) %>%
  mutate(middle = (windowStart + windowEnd)/2) %>%
  ggplot(.,aes(x=middle/1000000,y=D)) +
  geom_point(alpha=0.2,color="#a8dadc") +
  geom_line(aes(x=middle/1000000,y=mean_D)) +
  facet_wrap(~chr) +
  geom_hline(aes(yintercept=0)) +
  scale_color_viridis_d() +
  theme_cowplot() +
  xlab("MBp") +
  ggtitle("Ann-PetPet-Des") +
  theme(axis.text = element_text(size=8))


p3.2 + inset_element(p3.1, left = 0.4, bottom = 0, right = 1, top = 0.2)

data <- read_tsv("../d/ann_pet_fal_des_localFstats__50_25.txt")

p4.1 <- data %>%
  ggplot(aes(D)) +
  geom_histogram(fill="#a8dadc") +
  theme_cowplot() +
  theme(axis.text = element_text(size=8)) +
  coord_cartesian(xlim=c(-1,1))


p4.2 <- data %>%
  mutate(chr = gsub("Ha412HO","",chr)) %>%
  mutate(mean_D = rollapply(D, 10, mean, fill = NA, align = "right")) %>%
  mutate(middle = (windowStart + windowEnd)/2) %>%
  ggplot(.,aes(x=middle/1000000,y=D)) +
  geom_point(alpha=0.2,color="#a8dadc") +
  geom_line(aes(x=middle/1000000,y=mean_D)) +
  facet_wrap(~chr) +
  geom_hline(aes(yintercept=0)) +
  scale_color_viridis_d() +
  theme_cowplot() +
  xlab("MBp") +
  ggtitle("Ann-PetFal-Des") +  
  theme(axis.text = element_text(size=8))


p4.2 + inset_element(p4.1, left = 0.4, bottom = 0, right = 1, top = 0.2)

data <- read_tsv("../d/ann_niv_des_localFstats__50_25.txt")

p7.1 <- data %>%
  ggplot(aes(D)) +
  geom_histogram(fill="#a8dadc") +
  theme_cowplot() +
  theme(axis.text = element_text(size=8)) +
  coord_cartesian(xlim=c(-1,1))


p8.2 <- data %>%
  mutate(chr = gsub("Ha412HO","",chr)) %>%
  mutate(mean_D = rollapply(D, 10, mean, fill = NA, align = "right")) %>%
  mutate(middle = (windowStart + windowEnd)/2) %>%
  ggplot(.,aes(x=middle/1000000,y=D)) +
  geom_point(alpha=0.2,color="#a8dadc") +
  geom_line(aes(x=middle/1000000,y=mean_D)) +
  facet_wrap(~chr) +
  geom_hline(aes(yintercept=0)) +
  scale_color_viridis_d() +
  theme_cowplot() +
  xlab("MBp") +
  ggtitle("Ann-Niv-Des") +
  theme(axis.text = element_text(size=8))

p8.2 + inset_element(p1.1, left = 0.4, bottom = 0, right = 1, top = 0.2)

data <- read_tsv("../d/ann_pet_pet_par_localFstats__50_25.txt")

p5.1 <- data %>%
  ggplot(aes(D)) +
  geom_histogram(fill="#a8dadc") +
  theme_cowplot() +
  theme(axis.text = element_text(size=8)) +
  coord_cartesian(xlim=c(-1,1))


p5.2 <- data %>%
  mutate(chr = gsub("Ha412HO","",chr)) %>%
  mutate(mean_D = rollapply(D, 10, mean, fill = NA, align = "right")) %>%
  mutate(middle = (windowStart + windowEnd)/2) %>%
  ggplot(.,aes(x=middle/1000000,y=D)) +
  geom_point(alpha=0.2,color="#a8dadc") +
  geom_line(aes(x=middle/1000000,y=mean_D)) +
  facet_wrap(~chr) +
  geom_hline(aes(yintercept=0)) +
  scale_color_viridis_d() +
  theme_cowplot() +
  xlab("MBp") +
  ggtitle("Ann-PetPet-Par") +
  theme(axis.text = element_text(size=8))


p5.2 + inset_element(p5.1, left = 0.4, bottom = 0, right = 1, top = 0.2)

data <- read_tsv("../d/ann_pet_fal_par_localFstats__50_25.txt")

p6.1 <- data %>%
  ggplot(aes(D)) +
  geom_histogram(fill="#a8dadc") +
  theme_cowplot() +
  theme(axis.text = element_text(size=8)) +
  coord_cartesian(xlim=c(-1,1))


p6.2 <- data %>%
  mutate(chr = gsub("Ha412HO","",chr)) %>%
  mutate(mean_D = rollapply(D, 10, mean, fill = NA, align = "right")) %>%
  mutate(middle = (windowStart + windowEnd)/2) %>%
  ggplot(.,aes(x=middle/1000000,y=D)) +
  geom_point(alpha=0.2,color="#a8dadc") +
  geom_line(aes(x=middle/1000000,y=mean_D)) +
  facet_wrap(~chr) +
  geom_hline(aes(yintercept=0)) +
  scale_color_viridis_d() +
  theme_cowplot() +
  xlab("MBp") +
  ggtitle("Ann-PetFal-Par") +
  theme(axis.text = element_text(size=8))


p6.2 + inset_element(p6.1, left = 0.4, bottom = 0, right = 1, top = 0.2)

data <- read_tsv("../d/ann_niv_par_localFstats__50_25.txt")

p9.1 <- data %>%
  ggplot(aes(D)) +
  geom_histogram(fill="#a8dadc") +
  theme_cowplot() +
  theme(axis.text = element_text(size=8)) +
  coord_cartesian(xlim=c(-1,1))


p9.2 <- data %>%
  mutate(chr = gsub("Ha412HO","",chr)) %>%
  mutate(mean_D = rollapply(D, 10, mean, fill = NA, align = "right")) %>%
  mutate(middle = (windowStart + windowEnd)/2) %>%
  ggplot(.,aes(x=middle/1000000,y=D)) +
  geom_point(alpha=0.2,color="#a8dadc") +
  geom_line(aes(x=middle/1000000,y=mean_D)) +
  facet_wrap(~chr) +
  geom_hline(aes(yintercept=0)) +
  scale_color_viridis_d() +
  theme_cowplot() +
  xlab("MBp") +
  ggtitle("Ann-Niv-Par") +
  theme(axis.text = element_text(size=8))

p9.2 + inset_element(p1.1, left = 0.4, bottom = 0, right = 1, top = 0.2)



dev.off()

