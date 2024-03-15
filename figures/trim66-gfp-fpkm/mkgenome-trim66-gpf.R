library(Biostrings)

dir.create("references/trim66-gfp")
download.file("https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/latest/mm10.fa.gz",
              destfile = "references/trim66-gfp/mm10.fa.gz")
genome_sequence <- readDNAStringSet("references/trim66-gfp/mm10.fa")

gfp <- readDNAStringSet("trim66-gfp/gfp.fasta")
gfp_rc <- reverseComplement(gfp)

length(genome_sequence[["chr7"]])
# 145441459

length(gfp_rc[[1]])
# 875
# expected: 145441459 + 875 = 145442334

chr7 <- genome_sequence[["chr7"]]
upstream <- subseq(chr7, start = 1, end = 109484668)
downstream <- subseq(chr7, start = 109484669, end = length(chr7))
result <- paste0(upstream, gfp_rc, downstream)

nchar(result)
# 145442334

genome_sequence[["chr7"]] <- result

writeXStringSet(genome_sequence, "references/trim66-gfp/mm10-gfp.fa")
