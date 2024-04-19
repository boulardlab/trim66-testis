library(pheatmap)
library(RColorBrewer)

input_files <- as.character(snakemake@input)
dat <- lapply(input_files, scan, sep = "\t", what = "character")
dat <- as.data.frame(do.call(rbind, dat), stringsAsFactors = FALSE)
colnames(dat) <- c("s1", "s2", "cor")

mat <- reshape(dat, idvar = "s1", timevar = "s2", direction = "wide")
rownames(mat) <- mat[, 1]
mat <- mat[, -1]
mat <- apply(mat, c(1, 2), as.numeric)
mat <- as.matrix(mat)

pheatmap(
    mat,
    filename = as.character(snakemake@output),
    cellwidth = 20,
    cellheight = 20,
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
)
