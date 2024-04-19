library(tidyverse)
library(Seurat)
library(SeuratObject)

root <- "/scratch/tabaro/scrna-seq/testis/data"

# Human
gse142585_h <- readRDS(file.path(root, "GSE142585/seurat-human-one2one.rds"))

# Mouse
gse112393   <- readRDS(file.path(root, "GSE112393/seurat-one2one.rds"))

# Monkey
gse142585_m <- readRDS(file.path(root, "GSE142585/seurat-monkey-one2one.rds"))

all_dat <- list(
  gse142585_h,
  gse112393,
  gse142585_m
)

features <- SelectIntegrationFeatures(object.list = all_dat,
                                      selection.method = "vst",
                                      assay = rep("SCT", length(all_dat)),
                                      nfeatures = 2000)

all_dat <- PrepSCTIntegration(all_dat,
                              assay = "SCT",
                              anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = all_dat,
                                  anchor.features = features,
                                  normalization.method = "SCT",
                                  reduction = "cca",
                                  l2.norm = TRUE,
                                  nn.method = "annoy",
                                  n.trees = 100)

integrated <- IntegrateData(anchorset = anchors, 
                            normalization.method = "SCT")

saveRDS(integrated, file.path(root, "pmid32504559-integrated.rds"))

