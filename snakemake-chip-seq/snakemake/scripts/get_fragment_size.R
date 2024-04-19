##########
# This script uses cross-correlation plots to estimate the fragment length of
# single-end experiments.
# Descostes December 2019
##########

library(csaw)

#############
## PARAMS
#############

bamFile <- as.character(snakemake@input$bamFile)
outFile <- as.character(snakemake@output$report)
outFigure <- as.character(snakemake@output$figure)
max.delay <- as.numeric(snakemake@params$max_delay)

cat("Input BAM: ", bamFile, "\n")
cat("Output txt: ", outFile, "\n")
cat("Output png: ", outFigure, "\n")
cat("max delay: ", max.delay, "\n")

#############
## MAIN
#############

param <- readParam(minq = 20)
dedup.on <- reform(param, dedup = TRUE)

x <- correlateReads(bamFile, max.delay, param = dedup.on)

png(filename = outFigure)
plot(0:max.delay, x, type = "l", ylab = "CCF", xlab = "Delay (bp)")
dev.off()

extensionSize <- maximizeCcf(x)

write(extensionSize, file =  outFile, ncolumns = 1)
Sys.sleep(10)
