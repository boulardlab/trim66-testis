library(ggplot2)

dat_path <- as.character(snakemake@input[[1]])
dat <- read.table(dat_path, header = TRUE)

print(dat)

p <- ggplot(dat, aes(log2FoldChange, -log10(padj))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05)) +
    geom_vline(xintercept = 0) +
    xlim(-3, 3) +
    theme_bw() +
    theme(aspect.ratio = 1)

ggsave(snakemake@output[[1]], p)
