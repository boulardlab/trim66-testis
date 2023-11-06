library(Seurat)
library(tidyverse)
library(ggplot2)
library(viridis)
library(ggrepel)

rename_celltype <- function(x) {
  # Fix mouse CellType
  ### Cell type:
  # 1-InnateLymph;
  # 2-Macrophage;
  # 3-Endothelial;
  # 4-Myoid;
  # 5-Leydig;
  # 6-Sertoli;
  # 7-Unknown;
  # 8-SPG;
  # 9-Scytes;
  # 10-STids;
  # 11-Elongating
  x[x == 1]  <- "InnateLymph"
  x[x == 2]  <- "Machrophage"
  x[x == 3]  <- "Endothelial"
  x[x == 4]  <- "Myoid"
  x[x == 5]  <- "ImmLeydig"
  x[x == 6]  <- "Sertoli"
  x[x == 7]  <- "Unknown"
  x[x == 8]  <- "Spermatogonia"
  x[x == 9]  <- "Spermatocyte"
  x[x == 10] <- "RoundSpermatid"
  x[x == 11] <- "Elongating"
  return(x)
}

integrated <- readRDS("/scratch/tabaro/scrna-seq/testis/data/pmid32504559-integrated.rds") %>%
  RunPCA() %>%
  RunUMAP(dims = 1:12)
integrated$CellType <- rename_celltype(integrated$CellType)
integrated$organism <- umap$organism
Idents(integrated) <- "CellType"

m <- GetAssayData(integrated, assay = "SCT", slot = "data")
trim66 <- m["ENSG00000166436",]

umap <- Embeddings(integrated, reduction = "umap") %>%
  as_tibble(rownames = "cell") %>%
  mutate(
    orig.ident = integrated$orig.ident,
    organism   = sub("(Human|Monkey).*", "\\1", orig.ident),
    organism   = ifelse(!grepl("Human|Monkey", organism), "Mouse", organism),
    ident      = integrated$CellType,
    ident      = rename_celltype(ident),
    Trim66     = trim66
  ) %>%
  filter(organism != "Monkey")

centroids <- umap %>%
  group_by(ident, organism) %>%
  summarize(X = mean(UMAP_1),
            Y = mean(UMAP_2))

p1 <- ggplot() +
  geom_point(data = umap, mapping = aes(UMAP_1, UMAP_2, color = ident), size = 0.125) +
  scale_color_viridis_d(option = "turbo") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  facet_grid(~ organism) +
  geom_label_repel(data = centroids,
                  mapping = aes(X, Y, label = ident),
                  color = "black",
                  size = 3,
                  max.overlaps = 10,
                  force = 5) +
  theme_classic(base_size = 18) +
  theme(aspect.ratio = 1,
        strip.background = element_blank())

p2 <- ggplot() +
  geom_point(data = umap, mapping = aes(UMAP_1, UMAP_2, color = Trim66), size = 0.125) +
  scale_color_viridis_c(direction = -1, option = "magma") +
  facet_grid(~ organism) +
  geom_label_repel(data = centroids,
                   mapping = aes(X, Y, label = ident),
                   color = "black",
                   size = 3,
                   max.overlaps = 10,
                   force = 5) +
  theme_classic(base_size = 18) +
  theme(aspect.ratio = 1,
        strip.background = element_blank())

library(patchwork)
p <- p1 / p2
ggsave(p, file = "/scratch/tabaro/scrna-seq/testis/data/pmid32504559-umaps.pdf",
       height = 12, width = 14)
