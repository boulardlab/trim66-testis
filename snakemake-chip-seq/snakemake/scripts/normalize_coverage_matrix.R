mat_path <- as.character(snakemake@input[["matrix"]])
output_path <- as.character(snakemake@output[[1]])
frag_size_path <- as.character(snakemake@input[["fragment_size"]])
nmapped_path <- as.character(snakemake@input[["number_of_mapped_reads"]])
effective_genome_size <- as.numeric(snakemake@params[["effective_genome_size"]])
sample_names <- as.character(snakemake@params[["sample_names"]])

fragsize <- sapply(frag_size_path, scan)
nmapped <- sapply(nmapped_path, scan)

mat <- read.delim(mat_path, header = FALSE)

coords <- mat[, c(1, 2, 3)]
dat <- as.matrix(mat[, seq(4, ncol(mat))])

scale_factor <- nmapped * fragsize / effective_genome_size
scale_factor <- 1 / scale_factor

# divide each column for its scaling factor
dat <- sweep(dat, 2, scale_factor, FUN = "*")

dat <- as.data.frame(dat)
ret <- cbind(coords, dat)
colnames(ret) <- c("chr", "start", "end", sample_names)

write.csv(
    ret,
    output_path,
    row.names = FALSE,
    quote = FALSE
    )