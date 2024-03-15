library(tidyverse)
library(ggpubr)

sample_sheet_trim66_gfp_total <- read_csv("trim66-gfp/sample_sheet_trim66-gfp_totalRNA-seq.csv")  %>% select(name, genotype)
sample_sheet_trim66_gfp_es <- read_csv("trim66-gfp/sample_sheet_trim66-elongatedSpermatids_Ago2021_noeGFP2_noeGFP3_noeGFP6.csv")  %>% select(name, genotype)

sample_sheet <- bind_rows(
  sample_sheet_trim66_gfp_es %>% mutate(library = "ES total RNA-seq"),
  sample_sheet_trim66_gfp_total %>% mutate(library = "RS total RNA-seq")
  )

ls <- list.files(
  "alignments/",
  pattern = "*.coverage.txt",
  recursive = TRUE,
  full.names = TRUE)

samtools_coverage <- ls %>%
  map(read_delim) %>%
  set_names(basename(ls)) %>%
  map(dplyr::select, numreads) %>%
  map2(.y = ls, ~mutate(.x, sample = sub(".*_(t?MBAO[0-9]+|eGFP[0-9]+)_.*", "\\1", basename(.y)))) %>%
  map(~left_join(.x, sample_sheet, by = c(sample = "name"))) %>%
  reduce(bind_rows)

ls <- list.files(
  "trim66-gfp/alignments/",
  pattern = "*.alignment_metrics.txt",
  recursive = TRUE,
  full.names = TRUE)

picard_alignment_stats <- ls %>%
  map(read_delim, skip = 6, n_max = 3) %>%
  map(filter, CATEGORY == "PAIR") %>%
  map(select, PF_READS_ALIGNED) %>%
  map2(.y = ls, ~mutate(.x, sample = sub(".*_(t?MBAO[0-9]+|eGFP[0-9]+)_.*", "\\1", basename(.y)))) %>%
  reduce(bind_rows) %>%
  rename(library_size = "PF_READS_ALIGNED")

p <- samtools_coverage %>%
  inner_join(picard_alignment_stats) %>%
  mutate(
    scaling_factor = library_size / 1e6,
    rpm = numreads / scaling_factor,
    fpkm = rpm / (875 / 1e3), # 875 is the length of GFP
    library = factor(
      library,
      levels = c("RS total RNA-seq", "ES total RNA-seq"))
  ) %>%
  select(-c(scaling_factor, rpm)) %>%
  ggplot(aes(genotype, fpkm)) +
  geom_boxplot(outlier.size = NA) +
  geom_jitter(width = 0.25) +
  stat_compare_means(comparisons = list(c("KO", "WT"))) +
  facet_grid(
    cols = vars(library),
    labeller = labeller(
      library = c(`ES total RNA-seq` = "Elongating spermatids",
                  `RS total RNA-seq` = "Round spermatids"))) +
  scale_x_discrete(name = "Genotype") +
  scale_y_continuous(name = "FPKM") +
  ggtitle("GFP signal in Trim66 mutants") +
  theme_bw(base_size = 16) +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(size = 16))
ggsave("trim66-gfp/gfp-fpkm.pdf", plot = p, width = 8, height = 8)

