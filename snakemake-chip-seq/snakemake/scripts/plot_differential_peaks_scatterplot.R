sample_sheet <- as.data.frame(snakemake@params[["sample_sheet"]])
contrast_column <- as.character(snakemake@params[["contrast_column"]])
contrast_levels <- as.character(snakemake@params[["contrast_levels"]])
samples_labels <- sample_sheet[, contrast_column]

matrix_path <- as.character(snakemake@input[[1]])
mat <- read.csv(matrix_path, header = TRUE, check.names = FALSE)
rownames(mat) <- sprintf("%s:%d-%d", mat$chr, mat$start, mat$end)
print(head(mat))
coords <- mat[, c("chr", "start", "end")]
mat[, c("chr", "start", "end")] <- NULL
mat <- mat[, match(colnames(mat), sample_sheet$sample_name)]

samples_label <- factor(samples_labels, levels = contrast_levels)
samples_label <- as.numeric(samples_label)

g1 <- as.matrix(mat[, samples_label == 1])
g2 <- as.matrix(mat[, samples_label == 2])

dat <- data.frame(
    m1 = apply(g1, 1, mean),
    m2 = apply(g2, 1, mean))

fit <- lm(m2 ~ m1, data = dat)
r2 <- summary(fit)$r.squared

pdf(
    snakemake@output[[1]],
    height = 7,
    width = 7
)
par(
    pty = "s",
    cex = 1.4
)
plot(
    apply(g1, 1, mean),
    apply(g2, 1, mean),
    xlab = contrast_levels[1],
    ylab = contrast_levels[2],
    pch = 19,
    sub = bquote("R"^2 * " = " * .(round(r2, 4))),
    main = "Average normalized peak intensity values"
)
abline(coef = c(0, 1))
abline(coef = fit$coefficients, col = "blue", lwd = 2)
dev.off()