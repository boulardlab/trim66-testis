brary(tidyverse)
library(rstatix)
library(ggpubr)

dat <- tibble(
  sample = seq(17204, 17209),
  npeaks = c(39, 58, 46, 59, 56, 82),
  genotype = c("WT", "KO", "KO", "WT", "WT", "KO")
)


 
dat %>%
  mutate(genotype = as.factor(genotype)) %>%
  ggplot(aes(x = genotype, y = npeaks)) +
  geom_dotplot(
    binaxis='y', 
    stackdir='center',
    # stackratio=1.5, 
    dotsize=0.75
    ) +
  stat_summary(
    fun.data=mean_sdl,
    fun.args=list(mult=1),
    # geom="crossbar",
    size=1,
    color="red",
    ) +
  stat_compare_means(
    method = "t.test", 
    comparisons = list(c("WT", "KO")),
    label = "p.format",
    size = 5,
    bracket.size = 1.2
    ) +
  xlab("Genotype") +
  ylab("Number of peaks") +
  theme_classic(base_size = 24) + 
  theme(
    aspect.ratio = 2
  )
ggsave("~/Pictures/trim66-sperm/npeaks-dotplot.pdf", height = 8, width = 4)
