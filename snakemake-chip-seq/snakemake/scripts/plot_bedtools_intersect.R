intersect_results <- as.character(snakemake@input[[1]])
plot_path <- as.character(snakemake@output)

dat <- read.delim(intersect_results, header = FALSE)


pdf(
    plot_path, height = 5, width = 5
)
hist(
    dat[, 4],
    main = "Number of samples over consensus peaks",
    xlab = "Number of samples"
    )
dev.off()
