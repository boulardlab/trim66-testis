library(rtracklayer)
library(Gviz)

scheme <- getScheme(name = "default")

# Palette: https://www.ibm.com/design/language/color/
scheme$GdObject$background.title <- "transparent"
scheme$GdObject$col.axis <- "black"
scheme$GdObject$fontcolor.title <- "black"
scheme$GdObject$min.distance <- 6
scheme$GdObject$cex.title <- 0.7
scheme$GdObject$cex.axis <- 0.6
# scheme$GdObject$rotation.title <- 360

scheme$DataTrack$col.grid <- '#e0e0e0'
scheme$DataTrack$col.histogram <- '#4589ff'
scheme$DataTrack$col.boxplotFrame <- NA
scheme$DataTrack$fill.histogram <- '#a6c8ff'
scheme$DataTrack$col.border.title <- NA
scheme$DataTrack$lwd.border.title <- 0
scheme$DataTrack$jitter.x <- TRUE

scheme$GeneRegionTrack$col <- NA
scheme$GeneRegionTrack$col.border.title <- NA
scheme$GeneRegionTrack$lwd.border.title <- 0
scheme$GeneRegionTrack$fill <- '#6f6f6f'
scheme$GeneRegionTrack$gfp <- "#42be65"
scheme$GeneRegionTrack$default <- scheme$GeneRegionTrack$fill

scheme$GenomeAxisTrack$add35 <- TRUE
scheme$GenomeAxisTrack$add53 <- TRUE
scheme$GenomeAxisTrack$col <- '#6f6f6f'
scheme$GenomeAxisTrack$col.range <- NA
scheme$GenomeAxisTrack$fill.range <- '#c1c7cd'

scheme$AnnotationTrack$col <- NA

addScheme(scheme, "myScheme")
options(Gviz.scheme = "myScheme", ucscChromosomeNames=FALSE)

if (! "gr" %in% ls()) {
  gr <- import("data/references/MGI.gff3")
  seqlevelsStyle(gr) <- "UCSC"
  genome(gr) <- "mm10"

  chr <- as.character(unique(seqnames(gr)))
  gen <- genome(gr)
}

grr <- import("data/references/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff.gz")
seqlevelsStyle(grr) <- "UCSC"
sele <- subset(grr, seqnames == "chr7" & start > 109430000 & end < 109520000)
RegTrack <- AnnotationTrack(sele, stacking = "dense", showOverplotting = TRUE,
                            name = "Ensembl regbuild", rotation.title=360,
                            shape = "box",
                            `Promoter Flanking Region`="#1192e8",
                            `CTCF Binding Site`="#3ddbd9",
                            `Enhancer`="#ffafd2",
                            `Promoter`="#82cfff",
                            `Open chromatin`="#e8daff",
                            `TF binding site`="#78a9ff")
feature(RegTrack) <- sele$feature_type



# trim66 <- subset(gr, gene_id == "MGI:2152406" & type == "gene")
# export(trim66, con = "data/references/Trim66_MGI.gtf")
# awk -F "\t" -v OFS="\t" '$3~"exon"{split($9, a, ";"); print $1,$4,$5,$6,$7,a[1]}' Trim66_MGI.gtf | sed -r 's/ID "(.+)"/\1/g' > Trim66_MGI_exons.bed
# sort -k2,2n Trim66_MGI_exons.bed  > Trim66_MGI_exons_sorted.bed
# SLOP=20; bedtools slop -b $SLOP -i Trim66_MGI_exons_sorted.bed -g GRCm38.primary_assembly.genome.chromsize  > Trim66_MGI_exons_slop${SLOP}.bed
# bedtools complement -g GRCm38.primary_assembly.genome.chromsize -i Trim66_MGI_exons_slop${SLOP}.bed > Trim66_MGI_exons_complement.bed
# introns <- import("data/references/Trim66_MGI_complement.bed")

gTrack <- GenomeAxisTrack(gr)

# myGeneModels <- subset(as.data.frame(gr),
#                        subset = !is.na(exon_id),
#                        select = c("seqnames", "start", "end", "width", "strand",
#                                   "transcript_name", "gene_id", "gene_name",
#                                   "exon_id"))
myGeneModels <- subset(as.data.frame(gr),
                       subset = !is.na(exon_id),
                       select = c("seqnames", "start", "end", "width", "strand",
                                  "transcript_id", "gene_id",  "exon_id"))
myGeneModels[myGeneModels$transcript_id == "ENSMUST00000033339", "transcript_id"]  <- "Trim66-201"
myGeneModels[myGeneModels$transcript_id == "ENSMUST00000106739", "transcript_id"]  <- "Trim66-202"
myGeneModels[myGeneModels$transcript_id == "ENSMUST00000106741", "transcript_id"]  <- "Trim66-203"
myGeneModels[myGeneModels$transcript_id == "ENSMUST00000137704", "transcript_id"]  <- "Trim66-204"

# colnames(myGeneModels) <- c("chromosome", "start", "end", "width", "strand",
#                             "transcript", "gene", "symbol", "exon")
colnames(myGeneModels) <- c("chromosome", "start", "end", "width", "strand",
                            "transcript", "gene", "exon")

myGeneModels$color <- "default"
myGeneModels[myGeneModels$exon == "ENSMUSE00000205644", "color"] <- "gfp"
grTrack <- GeneRegionTrack(myGeneModels,
                           genome = "mm10",
                           chromosome = "chr7",
                           feature=as.vector(myGeneModels$color),
                           name = "MGI")

ls <- list.files("data/rna-seq/alignments/star/TRIM66_Dec2020/",
                 pattern = "Aligned.sortedByCoord.out.bam$", full.names = TRUE)
dTracks <- sapply(ls, function(p) {
  n <- gsub(".*(MBAO[0-9]+).*$", "\\1", p)
  AlignmentsTrack(p, isPaired = TRUE, name=n, ylim = c(0,100))
  # DataTrack(p, genome = "mm10", type = c("g", "l"), name = n, window = -1,  ylim = c(0,100))
})


# introns <-subset(gr, gene_id == "MGI:2152406" & type == "exon")
# introns <- reduce(introns)
# introns <- gaps(introns)
# intronsTrack <- AnnotationTrack(introns, shape = "box", name = "Trim66 introns")


pdf("analysis/pictures/trim66-gfp_insertion.pdf", height = 20, width = 8)
# plotTracks(c(gTrack, grTrack, intronsTrack, RegTrack, dTracks),
# sizes = c(1, 1, 1, 0.5, rep(5, length(dTracks))),
plotTracks(c(gTrack, grTrack, RegTrack, dTracks),
           sizes = c(1, 1, 0.5, rep(5, length(dTracks))),
           chromosome = "chr7",
           from = 109448000,
           to = 109517000,
           type = c("coverage", "sashimi"),
           transcriptAnnotation = "transcript",
           sashimiScore = 3,
           sashimiStrand = "-",
           lwd.sashimiMax = 6,
           # sashimiFilter = introns,
           # sashimiFilterTolerance = 5L,
           sashimiHeight = .9,
           minSashimiHeight = .8)
dev.off()


