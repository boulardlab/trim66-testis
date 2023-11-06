library(grid)
library(gridBase)
library(gridExtra)

input_files <- as.character(snakemake@input)

sample_sheet <- as.data.frame(snakemake@params[["sample_sheet"]])
contrast_column <- as.character(snakemake@params[["contrast_column"]])
contrast_levels <- as.character(snakemake@params[["contrast_levels"]])

output_file <- as.character(snakemake@output)


samples_labels <- sample_sheet[, contrast_column]

samples_label <- factor(samples_labels, levels = contrast_levels)
samples_label <- as.numeric(samples_label)
names(samples_label) <- as.character(sample_sheet$sample_name)
samples_label[samples_label == 1] <- "firebrick"
samples_label[samples_label == 2] <- "forestgreen"

npeaks <- sapply(seq_along(input_files), function(i) {
    p <- input_files[i]

    bn <- basename(dirname(p))
    print(bn)
    i <- grep(bn, sample_sheet[, "filename"])
    sample <- sample_sheet[i, "sample_name"]

    npeaks <- scan(p, what = numeric())

    return(setNames(npeaks, as.character(sample)))
})

pdf(
    output_file,
    width = 10,
    height = 5
    )
par(mfcol = c(1, 2))
barplot(
    rev(npeaks),
    col = samples_label,
    main = "Number of peaks",
    horiz = TRUE,
    las = 1
)
legend(
    "bottomright",
    fill = setNames(c("firebrick", "forestgreen"), levels(samples_labels)),
    legend = contrast_levels,
    bty = "n"
)

plot.new()
vps <- baseViewports()
pushViewport(vps$figure)
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params = list(cex = 1.6)),
  colhead = list(fg_params = list(cex = 1.6)),
  rowhead = list(fg_params = list(cex = 1.6)))
df <- as.data.frame(npeaks)
colnames(df) <- "Number of peaks"
tgb <- tableGrob(df, theme = mytheme)
grid.draw(tgb)
upViewport()
dev.off()
