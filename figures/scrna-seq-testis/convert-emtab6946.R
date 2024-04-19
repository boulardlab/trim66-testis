library(Seurat)
library(glmGamPoi)
library(sctransform)

#cell_metadata <- read.table("/scratch/tabaro/scrna-seq/testis/data/E-MTAB-6946/cell_metadata.txt")
#write.table(cell_metadata, row.names = F, file = "/scratch/tabaro/scrna-seq/testis/data/E-MTAB-6946/cell_metadata_clean.txt")

mtx   <- "/scratch/tabaro/scrna-seq/testis/data/E-MTAB-6946/raw_counts.mtx"
genes <- "/scratch/tabaro/scrna-seq/testis/data/E-MTAB-6946/genes.tsv"
cells <- "/scratch/tabaro/scrna-seq/testis/data/E-MTAB-6946/cell_metadata_clean.txt"

m <- ReadMtx(
  mtx = mtx,
  features = genes,
  cells = cells,
  cell.column = 2,
  cell.sep = " ",
  skip.cell = 1,
  skip.feature = 1)

cells.meta <- read.table(cells, header = TRUE)
cell_names <- paste(cells.meta$Library,
                    cells.meta$Sample,
                    cells.meta$Barcode,
                    cells.meta$AnnotatedClusters,
                    sep = "|")
rownames(cells.meta) <- cell_names
colnames(m) <- cell_names

m1 <- CreateSeuratObject(
  counts = m,
  project = basename(dirname(mtx)),
  meta.data = cells.meta,
  names.field = 4,
  names.delim = "|")
Idents(m1) <- "AnnotatedClusters"

m1 <- SCTransform(m1,
                  method = "glmGamPoi",
                  vst.flavor = "v2",
                  verbose = FALSE)
m1 <- RunPCA(m1)
m1 <- RunUMAP(m1, dims = 1:50)
m1 <- RunTSNE(m1, pca =FALSE, perplexity=350)
# m1 <- FindVariableFeatures(m1, selection.method = "vst")
saveRDS(m1, "/scratch/tabaro/scrna-seq/testis/data/E-MTAB-6946/seurat.rds")

