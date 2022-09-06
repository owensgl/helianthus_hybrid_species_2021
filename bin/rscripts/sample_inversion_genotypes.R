#Check inversion status of samples

ann <- read_tsv("/home/owens/bin/ha412_lostruct/Ha412HO_inv.v3.pcasites.annuus.inversiongenotypes.txt")

ann %>%
  filter(sample == "ANN1029" | sample == "ANN1283") %>%
  select(sample,triangle_genotype,sv_name) %>% View()

pet <- read_tsv("/home/owens/bin/ha412_lostruct/Ha412HO_inv.v3.pcasites.petiolaris.inversiongenotypes.txt")

pet %>%
  filter(sample == "PET0568" | sample == "PET0495" | sample == "PET0765" | sample == "PET0424") %>%
  select(sample,triangle_genotype,sv_name) %>% 
  pivot_wider(names_from=sv_name,values_from=triangle_genotype) %>% View()

arg <- read_tsv("/home/owens/bin/ha412_lostruct/Ha412HO_inv.v3.pcasites.argophyllus.inversiongenotypes.txt")

arg %>%
  filter(sample == "ARG0143" | sample == "ARG0295") %>%
  select(sample,triangle_genotype,sv_name) %>% 
  pivot_wider(names_from=sv_name,values_from=triangle_genotype) %>% View()
