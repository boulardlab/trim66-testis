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
con     <- file.path(root, "GSE142585/GSE142585_MergedHumanTestis4_DGE.txt")
genes   <- file.path(root, "GSE142585/GSE142585_MergedHumanTestis4_Genes.txt")
coldata <- file.path(root, "GSE142585/GSE142585_MergedHumanTestis4_PerCellAttributes.txt")

ori <- read.table(con)
log_info("Loaded ", con)

if (do_subsetting) {
  one2one <- read.delim(file.path(root, "biomart/human2macaque2mouse.txt"))
  log_info("Loaded mapping")

  k <- rownames(ori)[rownames(ori) %in% one2one$human_gene_name]
  ori <- ori[k,]
  rownames(ori) <- one2one$human_gene_id[match(k, one2one$human_gene_name)]
  outfile <- file.path(root, "GSE142585/seurat-human-one2one.rds")
  log_info("Done mapping")

  min_features <- 600
} else {
  outfile <- file.path(root, "GSE142585/seurat-human-original.rds")
  min_features <- 0
}

m2 <- CreateSeuratObject(
  counts = ori,
  meta.data = read.table(coldata),
  project = basename(dirname(con)),
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
