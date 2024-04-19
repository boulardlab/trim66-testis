library(tidyverse)

root <- "/scratch/tabaro/scrna-seq/testis/data/biomart"

macaque_mapping <- file.path(root, "macaque-one2one.txt")
mouse_mapping <- file.path(root, "mouse-one2one.txt")

human2macaque <- read_delim(macaque_mapping, col_names = F) %>%
  select(-X5) %>%
  rename_with(~c("human_gene_id", "human_gene_name",
                 "macaque_gene_id", "macaque_gene_name"))

human2mouse <- read_delim(mouse_mapping, col_names = F) %>%
  select(-X5) %>%
  rename_with(~c("human_gene_id", "human_gene_name",
                 "mouse_gene_id", "mouse_gene_name"))

one2one <- inner_join(human2macaque, human2mouse) %>%
  filter(!is.na(human_gene_name)) %>%
  mutate(
    macaque_gene_name = if_else(is.na(macaque_gene_name),
                                human_gene_name,
                                macaque_gene_name))

write_tsv(one2one, file.path(root, "human2macaque2mouse.txt"), na = ".")
