library(Seurat)
library(SingleCellExperiment)
library(glmGamPoi)
library(sctransform)
library(logger)

args <- commandArgs(trailingOnly = TRUE)
log_info("args: ", args)
do_subsetting <- as.numeric(args[1])
log_info("do_subsetting: ", do_subsetting)

root    <- "/scratch/tabaro/scrna-seq/testis/data"
con     <- file.path(root, "GSE142585/GSE142585_MergedMonkeyTestis5_DGE.txt")
genes   <- file.path(root, "GSE142585/GSE142585_MergedMonkeyTestis5_Genes.txt")
coldata <- file.path(root, "GSE142585/GSE142585_MergedMonkeyTestis5_PerCellAttributes.txt")

ori <- read.table(con)
log_info("Loaded ", con)


if (do_subsetting) {
  one2one <- read.delim(file.path(root, "biomart/human2macaque2mouse.txt"))
  log_info("Loaded mapping")

  k <- rownames(ori)[rownames(ori) %in% one2one$macaque_gene_name]
  m <- ori[k,]
  rownames(m) <- one2one$human_gene_id[match(k, one2one$macaque_gene_name)]
  ori <- m
  outfile <- file.path(root, "GSE142585/seurat-monkey-one2one.rds")
  log_info("Done mapping")

  min_features <- 600

} else {
  outfile <- file.path(root, "GSE142585/seurat-monkey-original.rds")
  min_features <- 0
}
log_info("min features: ", min_features)

m2 <- CreateSeuratObject(
  counts = ori,
  project = basename(dirname(con)),
  meta.data = read.table(coldata),
  names.field = 1,
  names.delim = "_",
  min.features = min_features)
log_info("Seurat object created")

Idents(m2) <- "CellType"

m2 <- SCTransform(m2,
                  method = "glmGamPoi",
                  vst.flavor = "v2",
                  verbose = FALSE)
log_info("SCT done")

m2 <- RunPCA(m2)
log_info("PCA done")

saveRDS(m2, outfile)
log_info("All done")

