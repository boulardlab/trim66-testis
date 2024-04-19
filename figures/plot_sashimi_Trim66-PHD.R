
# snakemake --directory trim66-sashimi/ --configfile trim66-sashimi/config.yaml --profile profile/slurm --snakefile 3t-seq/workflow/Snakefile \
# results/alignments/star/TRIM66_Oct2020/HVHTJBGXG_mbaj01_20s003046-1-1_Mielnicka_lane120s003046_sequence.Aligned.sortedByCoord.out.bam \
# results/alignments/star/TRIM66_Oct2020/HVHTJBGXG_mbaj02_20s003047-1-1_Mielnicka_lane120s003047_sequence.Aligned.sortedByCoord.out.bam \
# results/alignments/star/TRIM66_Oct2020/HVHTJBGXG_mbaj03_20s003048-1-1_Mielnicka_lane120s003048_sequence.Aligned.sortedByCoord.out.bam \
# results/alignments/star/TRIM66_Oct2020/HVHTJBGXG_mbaj04_20s003049-1-1_Mielnicka_lane120s003049_sequence.Aligned.sortedByCoord.out.bam \
# results/alignments/star/TRIM66_Oct2020/HVHTJBGXG_mbaj05_20s003050-1-1_Mielnicka_lane120s003050_sequence.Aligned.sortedByCoord.out.bam \
# results/alignments/star/TRIM66_Oct2020/HVHTJBGXG_mbaj06_20s003051-1-1_Mielnicka_lane120s003051_sequence.Aligned.sortedByCoord.out.bam \
# results/alignments/star/TRIM66_Oct2020/HVHTJBGXG_mbaj07_20s003052-1-1_Mielnicka_lane120s003052_sequence.Aligned.sortedByCoord.out.bam \
# results/alignments/star/TRIM66_Oct2020/HVHTJBGXG_mbaj08_20s003053-1-1_Mielnicka_lane120s003053_sequence.Aligned.sortedByCoord.out.bam \
# results/alignments/star/TRIM66_Oct2020/HVHTJBGXG_mbaj09_20s003054-1-1_Mielnicka_lane120s003054_sequence.Aligned.sortedByCoord.out.bam \
# results/alignments/star/TRIM66_Oct2020/HVHTJBGXG_mbaj10_20s003055-1-1_Mielnicka_lane120s003055_sequence.Aligned.sortedByCoord.out.bam \
# results/alignments/star/TRIM66_Oct2020/HVHTJBGXG_mbaj11_20s003056-1-1_Mielnicka_lane120s003056_sequence.Aligned.sortedByCoord.out.bam \
# results/alignments/star/TRIM66_Oct2020/HVHTJBGXG_mbaj12_20s003057-1-1_Mielnicka_lane120s003057_sequence.Aligned.sortedByCoord.out.bam

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
options(Gviz.scheme = "myScheme", ucscChromosomeNames = FALSE)

if (!"gr" %in% ls()) {
  gr <- import("references/MGI.2021.gff3")
  seqlevelsStyle(gr) <- "UCSC"
  genome(gr) <- "mm10"

  chr <- as.character(unique(seqnames(gr)))
  gen <- genome(gr)
}

grr <- import("references/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff.gz")
seqlevelsStyle(grr) <- "UCSC"
sele <- subset(grr, seqnames == "chr7" & start > 109430000 & end < 109520000)
RegTrack <- AnnotationTrack(sele, stacking = "dense", showOverplotting = TRUE,
                            name = "Ensembl regbuild", rotation.title = 360,
                            shape = "box",
                            `Promoter Flanking Region` = "#1192e8",
                            `CTCF Binding Site` = "#3ddbd9",
                            `Enhancer` = "#ffafd2",
                            `Promoter` = "#82cfff",
                            `Open chromatin` = "#e8daff",
                            `TF binding site` = "#78a9ff")
feature(RegTrack) <- sele$feature_type

gTrack <- GenomeAxisTrack(gr)

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

# myGeneModels$color <- "default"
# myGeneModels[myGeneModels$exon == "ENSMUSE00000205644", "color"] <- "gfp"
grTrack <- GeneRegionTrack(myGeneModels,
                           genome = "mm10",
                           chromosome = "chr7",
                           # feature=as.vector(myGeneModels$color),
                           name = "MGI")

ls <- list.files("trim66-sashimi/results/alignments/star_markdup/TRIM66_Oct2020/",
                 # pattern = "Aligned.sortedByCoord.out.bam$", full.names = TRUE)
                 pattern = ".bam$", full.names = TRUE)
dTracks <- sapply(ls, function(p) {
  n <- gsub(".*(mbaj[0-9]+).*$", "\\1", p)
  AlignmentsTrack(p, isPaired = TRUE, name = n, ylim = c(0,100))
})

introns <-subset(gr, gene_id == "MGI:2152406" & type == "exon")
introns <- reduce(introns)
introns <- gaps(introns)
intronsTrack <- AnnotationTrack(introns, shape = "box", name = "Trim66 introns")

pdf("trim66-sashimi/trim66-phd-null.pdf", height = 22, width = 8)
# plotTracks(c(gTrack, grTrack, intronsTrack, RegTrack, dTracks),
# sizes = c(1, 1, 1, 0.5, rep(5, length(dTracks))),
plotTracks(c(gTrack, grTrack, RegTrack, dTracks),
           sizes = c(1.25, 1, 0.5, rep(5, length(dTracks))),
           chromosome = "chr7",
           from = 109448000,
           to = 109517000,
           type = c("coverage", "sashimi"),
           transcriptAnnotation = "transcript",
           sashimiScore = 3,
           sashimiStrand = -1,
           lwd.sashimiMax = 4,
           # sashimiFilter = introns,
           # sashimiFilterTolerance = 5L,
           sashimiHeight = .9,
           minSashimiHeight = .8)
dev.off()
