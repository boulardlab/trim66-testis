library(data.table)
library(grDevices)

s1 <- as.character(snakemake@wildcards[["internal_sample"]])
s2 <- as.character(snakemake@wildcards[["external_sample"]])

dat <- fread(as.character(snakemake@input), header = FALSE)
dat <- dat[V4 + V5 > 0, ]

tau <- cor(
    dat$V4,
    dat$V5,
    method = "spearman"
)

writeLines(
    text = sprintf("%s\t%s\t%.6f", s1, s2, tau),
    con = as.character(snakemake@output[["corr"]]),
    sep = "\n"
)

# p <- ggplot(dat, aes(x = V4, y = V5)) +
#     geom_point() +
#     geom_abline(
#         slope = 1,
#         intercept = 0,
#         color = "dodgerblue",
#         size = 1.2
#     ) +
#     geom_smooth(method = "lm") +
#     xlab(s1) +
#     ylab(s2) +
#     theme_bw(base_size = 18) +
#     theme(aspect.ratio = 1)
# ggsave(
#     filename = as.character(snakemake@output[["plot"]]),
#     plot = p,
#     width = 8,
#     height = 8
# )

pdf(as.character(snakemake@output[["plot"]]), width = 8, height = 8)
palette <- heat.colors(30)
smoothScatter(
    log2(dat$V4),
    log2(dat$V5),
    xlab = s1,
    ylab = s2,
    pch = 10,
    col = "white",
    colramp = colorRampPalette(palette)
)
dev.off()
