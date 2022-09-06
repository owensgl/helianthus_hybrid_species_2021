#Find genomic location of capture data

genes <- read_tsv("../sequence_capture/genes.txt",col_names = "gene_number")
all_gene_locations <- tibble()
for (gene_number in genes$gene_number){
  read_tsv(paste0("../sequence_capture/gene_locations/",gene_number,".blastout.txt"),
           col_names = c("species","chr","id_per","length","mis1","mis2","mis3",
                         "length2","start","end","p","evalue")) %>%
    filter(evalue > 150) %>%
    arrange(desc(evalue)) %>%
    head(1) %>%
    select(chr,start,end) %>%
    mutate(gene_number = gene_number) -> tmp
  all_gene_locations <- rbind(all_gene_locations,tmp)
}

write_tsv(all_gene_locations,"../sequence_capture/gene_locations.txt")
