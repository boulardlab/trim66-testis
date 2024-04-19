sample_sheet <- as.data.frame(snakemake@params[["sample_sheet"]])
contrast_column <- as.character(snakemake@params[["contrast_column"]])
contrast_levels <- as.character(snakemake@params[["contrast_levels"]])
samples_labels <- sample_sheet[, contrast_column]

matrix_path <- as.character(snakemake@input[[1]])
mat <- read.csv(matrix_path, header = TRUE, check.names = FALSE)
rownames(mat) <- sprintf("%s:%d-%d", mat$chr, mat$start, mat$end)
coords <- mat[, c("chr", "start", "end")]
mat[, c("chr", "start", "end")] <- NULL
mat <- mat[, match(colnames(mat), sample_sheet$sample_name)]

samples_label <- factor(samples_labels, levels = contrast_levels)
samples_label <- as.numeric(samples_label)

g1 <- as.matrix(mat[, samples_label == 1])
g2 <- as.matrix(mat[, samples_label == 2])

nc1 <- ncol(g1)
nc2 <- ncol(g2)

mat <- cbind(g1, g2)

lfc <- apply(
    mat, 1,
    function(current_row) {
        median1 <- median(current_row[1:nc1] + 1)
        median2 <- median(current_row[(nc1 + 1):(nc1 + nc2)] + 1)
        log2(median1) - log2(median2)
    })

pvalue <- apply(
    mat, 1,
    function(current_row) {
        dat1 <- current_row[1:nc1]
        dat2 <- current_row[(nc1 + 1):(nc1 + 1 + nc2)]
        wilcox.test(dat1, dat2)$p.value
        })
qvalue <- p.adjust(pvalue, "fdr")

output <- data.frame(
    log2FoldChange = lfc,
    pvalue = pvalue,
    padj = qvalue
)
output <- cbind(coords, output)
write.table(
    output,
    snakemake@output[[1]],
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
    )
