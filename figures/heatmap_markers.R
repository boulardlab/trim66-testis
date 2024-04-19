set.seed(123)

library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(DESeq2)
library(rtracklayer)

## LOAD MARKERS

markers_original <- read_csv("analysis/pictures/publication/markes_with_cluster.csv", col_names = F) %>%
  dplyr::rename(gene_name = "X1",
                cluster = "X2")
markers <- read_tsv("analysis/pictures/publication/MGIBatchReport_20211110_101621.txt")
markers_with_cluster_annot <- markers %>%
  dplyr::filter(`MGI Gene/Marker ID` != "No associated gene") %>%
  inner_join(markers_original, by = c(Input = "gene_name"))

markers_to_plot <- markers_with_cluster_annot %>% pull(`MGI Gene/Marker ID`)
cluster_annot <- markers_with_cluster_annot %>% pull(cluster, name = "MGI Gene/Marker ID")

## LOAD DDS

dds_paths <- c(
  "analysis/rdata/deseq2/TRIM66_Apr2021/dds.rds",
  "analysis/rdata/deseq2/TRIM66_elongatedSpermatids_Ago2021_noeGFP2_noeGFP3_noeGFP6/dds.rds"
)

all_dds <- lapply(dds_paths, readRDS) %>%
  set_names(basename(dirname(dds_paths)))

## NORMALIZATION

mean_centering <- function(mat) {
  tmat <- t(mat)
  t(scale(tmat, center = apply(tmat, 2, mean), scale = apply(tmat, 2, sd)))
}

all_fpkm <- map(all_dds, fpkm)

k <- map(all_fpkm, rownames) %>%
  purrr::reduce(intersect)

z <- map(all_fpkm, function(x) x[k,]) %>%
  purrr::reduce(cbind) %>%
  mean_centering
colnames(z) <- sub(".*_((?:tMBAO|eGFP)[0-9]+)_.*", "\\1", colnames(z))

## SUBSET MARKERS AND BUILD rowAnnotation OBJECT

markers_mat <- z[rownames(z) %in% markers_to_plot,]
cluster_annot <- cluster_annot[rownames(markers_mat)]

left_annot <- rowAnnotation(Cluster = cluster_annot,
                            col = list(Cluster = c(C11 = "#55ffe6",
                                                   C12 = "#778beb")),
                            show_annotation_name = FALSE)

## BUILD HeatmapAnnotation

get_samples <- function(dds, g) colData(dds) %>%
  as_tibble(rownames = "rn") %>%
  filter(genotype == g) %>%
  dplyr::select(rn, name)

wt_rs <- get_samples(all_dds[[1]], "WT")
wt_es <- get_samples(all_dds[[2]], "WT")

ko_rs <- get_samples(all_dds[[1]], "KO")
ko_es <- get_samples(all_dds[[2]], "KO")

wt_samples <- bind_rows(wt_rs, wt_es)
ko_samples <- bind_rows(ko_rs, ko_es)

genotype <- vector("character", ncol(z))
genotype[colnames(z) %in% wt_samples$name] <- "WT"
genotype[colnames(z) %in% ko_samples$name] <- "KO"

top_annot <- HeatmapAnnotation(Stage = ifelse(grepl("eGFP", colnames(z)), "ES", "RS"),
                               Genotype = genotype,
                               col = list(Stage = c(ES = "#e66868",
                                                    RS = "#eb8686"),
                                          Genotype = c(KO = "#303a52",
                                                       WT = "#596174")))


################
## HISTONES
################

mgi <- import("/g/boulard/Francesco/projects/trim66/data/references/MGI_mod.gff3",
              colnames=c("Name", "description", "gene_id", "mgi_type", "type"),
              feature.type = c("gene", "pseudogene"))
names(mgi) <- mgi$gene_id

histones <- subset(mcols(mgi), (grepl("histone", description) &
                                  !grepl("deacetylase|linker|binding|chaperone|methyl|homolog|kinase|regulator|aminotransferase|reader|factor", description) &
                                  grepl("protein", mgi_type)))  %>%
  as_tibble() %>%
  mutate(variant = sub("^(H[1-4][A-Z]?).*", "\\1", description)) %>%
  filter(grepl("^H[1-4]", variant))

protamines <- subset(mcols(mgi), gene_id %in% c("MGI:97765", "MGI:97766", "MGI:106601")) %>%
  as_tibble

transition_proteins <- subset(mcols(mgi), gene_id %in% c("MGI:98784", "MGI:98785")) %>%
  as_tibble

ids <- c(histones$gene_id,
         transition_proteins$gene_id,
         protamines$gene_id)

## SUBSET FOR HISTONES/TRANSITION PROTEINS/PROTAMINES
m <- z[rownames(z) %in% ids,]

## DO HIERARCHICAL CLUSTERING ON THE ENTIRE MATRIX
hc <- hclust(dist(t(z[rownames(z) %in% c(ids, markers_to_plot),])))

histone_h2a <- histones %>%
  filter(variant == "H2A") %>%
  dplyr::select(-variant)

histone_h2b <- histones %>%
  filter(variant == "H2B") %>%
  dplyr::select(-variant)

histone_h3 <- histones %>%
  filter(variant == "H3") %>%
  dplyr::select(-variant)

histone_h4 <- histones %>%
  filter(variant == "H4") %>%
  dplyr::select(-variant)

pdf("analysis/pictures/publication/v1/markers_heatmap.pdf", height = 12, width = 8)
## MARKERS HEATMAP
mmh <- .0558
mmw <- 4
h1 <- Heatmap(markers_mat,
              row_title = "Marker genes",
              row_title_rot = 0,
              km = 2,
              name = "z-score",
              height = unit(mmh, "mm") * nrow(markers_mat),
              width = unit(mmw, "mm") * ncol(markers_mat),
              show_row_names = FALSE,
              cluster_columns = as.dendrogram(hc),
              left_annotation = left_annot,
              top_annotation = top_annot,
              use_raster = TRUE)



## HELPERS TO BUILD HISTONES/TRANSITION PROTEINS/PROTAMINES HEATMAPS
as_matrix <- function(df) {
  if ("description" %in% colnames(df)){
    rn <- paste(df %>% pull(Name),
                df %>% pull(description),
                sep =  " - ")
    mat <- df %>% dplyr::select(-c(Name, description)) %>% as.matrix
  }else{
    rn <- df %>% pull(Name)
    mat <- df %>% dplyr::select(-Name) %>% as.matrix
  }
  rownames(mat) <- rn
  return(mat)
}

get_matrix <- function(m, h) m %>%
  as_tibble(rownames = "rn") %>%
  inner_join(h, by = c(rn = "gene_id")) %>%
  dplyr::select(-c(rn, mgi_type, type)) %>%
  as_matrix


mmh <- 1

## HISTONE H2A HEATMAP
mh2a <- get_matrix(m, histone_h2a)
h1h <- Heatmap(mh2a,
              height = unit(mmh, "mm")*nrow(mh2a),
              width = unit(mmw, "mm")*ncol(mh2a),
              name = "z-score",
              cluster_columns = F,
              show_row_names = F,
              row_title = "Histone H2A",
              row_title_rot = 0,
              # top_annotation = ha,
              row_names_max_width = max_text_width(
                rownames(mh2a),
                gp = gpar(fontsize = 12)
              ),
              use_raster = TRUE)

## HISTONE H2B HEATMAP
mh2b <- get_matrix(m, histone_h2b)
h2 <- Heatmap(mh2b,
              height = unit(mmh, "mm")*nrow(mh2a),
              width = unit(mmw, "mm")*ncol(mh2a),
              name = "z-score",
              cluster_columns = F,
              show_row_names = F,
              row_title = "Histone H2B",
              row_title_rot = 0,
              row_names_max_width = max_text_width(
                rownames(mh2b),
                gp = gpar(fontsize = 12)
              ),
              use_raster = TRUE)

## HISTONE H3 HEATMAP
mh3 <- get_matrix(m, histone_h3)
h3 <- Heatmap(mh3,
              height = unit(mmh, "mm")*nrow(mh3),
              width = unit(mmw, "mm")*ncol(mh2a),
              name = "z-score",
              cluster_columns = F,
              show_row_names = F,
              row_title = "Histone H3",
              row_title_rot = 0,
              row_names_max_width = max_text_width(
                rownames(mh3),
                gp = gpar(fontsize = 12)
              ),
              use_raster = TRUE)

## HISTONE H4 HEATMAP
mh4 <- get_matrix(m, histone_h4)
h4 <- Heatmap(mh4,
              height = unit(mmh, "mm")*nrow(mh4),
              width = unit(mmw, "mm")*ncol(mh2a),
              name = "z-score",
              cluster_columns = F,
              show_row_names = F,
              row_title = "Histone H4",
              row_title_rot = 0,
              row_names_max_width = max_text_width(
                rownames(mh4),
                gp = gpar(fontsize = 12)
              ),
              use_raster = TRUE)

## TRANSITION PROTEINS HEATMAP
mtp <- get_matrix(m, transition_proteins)
h5 <- Heatmap(mtp,
              height = unit(mmh, "mm")*nrow(mtp),
              width = unit(mmw, "mm")*ncol(mh2a),
              name = "z-score",
              cluster_columns = F,
              row_title = "Transition proteins",
              row_title_rot = 0,
              show_row_names = F,
              row_names_max_width = max_text_width(
                rownames(mtp),
                gp = gpar(fontsize = 12)
              ),
              use_raster = TRUE)

## PROTAMINES HEATMAP
mp <- get_matrix(m, protamines)
h6 <- Heatmap(mp,
              height = unit(mmh, "mm")*nrow(mp),
              width = unit(mmw, "mm")*ncol(mh2a),
              name = "z-score",
              cluster_columns = F,
              row_title = "Protamines",
              row_title_rot = 0,
              show_row_names = F,
              row_names_max_width = max_text_width(
                rownames(mp),
                gp = gpar(fontsize = 12)
              ),
              use_raster = TRUE)

# PLOT
hl <- h1 %v% h1h %v% h2 %v% h3 %v% h4 %v% h5 %v% h6
draw(hl, padding =unit(c(14, 0, 14, 0), "mm") )
dev.off()
