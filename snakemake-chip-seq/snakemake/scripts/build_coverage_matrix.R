
input_files <- as.character(snakemake@input)
sample_sheet <- as.data.frame(snakemake@params[["sample_sheet"]])
output_file <- as.character(snakemake@output)

all_bed <- lapply(seq_along(input_files), function(i) {
    p <- input_files[i]

    bn <- sub(".bed", "", basename(p))
    i <- grep(bn, sample_sheet[, "filename"])
    sample <- sample_sheet[i, "sample_name"]

    bed <- read.delim(p, header = FALSE)
    bed <- bed[, c(1, 2, 3, 4)]
    colnames(bed) <- c("chr", "start", "end", as.character(sample))
    
    return(bed)
})

ret <- Reduce(merge, all_bed)

write.csv(
    ret,
    output_file,
    row.names = FALSE,
    quote = FALSE
    )