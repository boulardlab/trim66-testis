##########
## This script imports DESeq2 objects and output MA plots. Customized to work
## with a Snakemake pipeline output and compose a plot for figure in the
## publication.
##
## Author: Francesco Tabaro
## Email: francesco.tabaro@embl.it
##
##########

library(tidyverse)
library(ggrastr)
library(ggrepel)
library(DESeq2)

##########
##
## Some helper functions
##
##########

ggtheme <- function(legend.position = "bottom") {
  theme_classic(base_size = 14) %+replace%
    theme(aspect.ratio = 1,
          legend.position = legend.position,
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 15),
          legend.title = element_text(size = 15),
          plot.title = element_text(size = 16),
          strip.text = element_text(size = 14))
}

find_shape <- function(mean, lfc, qt) {
  shape <- "regular"
  if (mean > qt) {
    shape <- "right"
    if (lfc < -2) {
      # bottom
      shape <- "bottom-right"
    } else if (lfc > 2) {
      # top
      shape <- "top-right"
    }
  } else {
    # left side
    if (lfc < -2) {
      # bottom
      shape <- "bottom"
    } else if (lfc > 2) {
      # top
      shape <- "top"
    }
  }
  return(shape)
}

get_shape <- function(mean,
                      lfc,
                      qt_thr = 0.98,
                      find_shape_fun = find_shape) {

  qt <- quantile(mean, qt_thr)

  shape <- rep("regular", length(mean))

  for (i in seq_along(mean)) {
    m <- mean[i]
    l <- lfc[i]
    if (!any(is.na(c(m,l)))) shape[i] <- find_shape_fun(m, l, qt)
  }

  return(shape)

}

##########
##
## Paths to DESeq2 objects
##
##########

dds_paths <- c(
  "analysis/rdata/deseq2/TRIM66_Dec2020/dds.rds",
  "analysis/rdata/deseq2/TRIM66_Apr2021/dds.rds",
  "analysis/rdata/deseq2/TRIM66_elongatedSpermatids_Ago2021_noeGFP2_noeGFP3_noeGFP6/dds.rds"
)

##########
##
## Mapping for titles
##
##########

mains <- list(
  TRIM66_Apr2021 = "TRIM66-GFP\nRS total RNA-seq",
  TRIM66_Dec2020 = "TRIM66-GFP\nRS polyA RNA-seq",
  TRIM66_elongatedSpermatids_Ago2021_noeGFP2_noeGFP3_noeGFP6 = "TRIM66-GFP\nES total RNA-seq"
)

##########
##
## Import DESeq2 objects
##
##########

all_dds <- lapply(dds_paths, readRDS) %>% set_names(basename(dirname(dds_paths)))

##########
##
## Compute log2 fold change
##
##########

all_lfcShrink <- lapply(all_dds, lfcShrink,
                        coef = "genotype_KO_vs_WT",
                        type = "apeglm") %>%
  set_names(names(all_dds))

##########
##
## Extract DE genes
##
##########

get_deg <- function(dataset){
  dds <- all_dds[[dataset]]
  df <- all_lfcShrink[[dataset]]

  rr <- rowRanges(dds)
  rr <- as.data.frame(mcols(rr)) %>%
    dplyr::select(gene_id, Name, description, mgi_type)

  sbst <- subset(df, padj < 0.05)

  return(as_tibble(sbst, rownames = "gene_id") %>%
           left_join(rr) %>%
           dplyr::rename(ID = "gene_id"))
}

all_deg_tables <- lapply(names(all_dds), get_deg) %>% set_names(names(all_dds))

##########
##
## MA plot
##
##########

select_genes <- function (gene_ids, isDE) {
  ret <- rep("none", length(gene_ids))
  ret[isDE] <- "DE"

  selected_genes <- c(
    Setd1b  = "MGI:2652820",
    Kmt2e   = "MGI:1924825",
    Kmt2d   = "MGI:2682319",
    Brd4    = "MGI:1888520",
    Ep400   = "MGI:1276124",
    Chd2    = "MGI:2448567",
    Chd3    = "MGI:1344395",
    Arid4a  = "MGI:2444354",
    Dnmt3l  = "MGI:1859287",
    Chd2    = "MGI:2448567",
    Chd3    = "MGI:1344395",
    Smarca4 = "MGI:88192"

    )
  trim66 <- c(Trim66 = "MGI:2152406")

  for (i in seq_along(selected_genes)) {
    ii <- which(gene_ids == selected_genes[i])
    ret[ii] <- names(selected_genes)[i]
  }

  ii <- which(gene_ids == trim66)
  ret[ii] <- "Trim66"

  return(ret)
}

p <- lapply(all_lfcShrink, DESeq2::plotMA, returnData = TRUE, alpha = 0.05) %>%
  map(as_tibble, rownames = "gene_id") %>%
  map2(.y = names(all_lfcShrink), ~mutate(.x, dataset = .y)) %>%
  purrr::reduce(bind_rows) %>%
  mutate(shape = get_shape(mean, lfc),
         mean = if_else(mean > quantile(mean, 0.98), quantile(mean, 0.98,), mean),
         lfc = if_else(lfc < -2, -2, if_else(lfc > 2, 2, lfc)),
         color =  select_genes(gene_id, isDE),
         label = if_else(!color %in% c("DE","none"), color, NA_character_)) %>%
  ggplot(aes(x = mean, y = lfc, color = color, shape = shape, label = label)) +
  rasterize(geom_point(data = . %>% dplyr::filter(color == "none"), size = 2.8), dpi = 500) +
  rasterize(geom_point(data = . %>% dplyr::filter(color == "DE"), size = 2.8), dpi = 500) +
  rasterize(geom_point(data = . %>% dplyr::filter(!color %in% c("DE", "none")), size = 2.8),
            dpi = 500) +
  geom_text_repel(color = "black", force = 10) +
  geom_hline(yintercept = 0, color = "black") +
  scale_color_manual(values = c(DE = "red",
                                Trim66 = "blue",
                                Setd1b  = "black",
                                Kmt2e   = "black",
                                Kmt2d   = "black",
                                Brd4    = "black",
                                Ep400   = "black",
                                Chd2    = "black",
                                Chd3    = "black",
                                Arid4a  = "black",
                                Dnmt3l  = "black",
                                Chd2    = "black",
                                Chd3    = "black",
                                Smarca4 = "black",
                                none = rgb(203/255, 203/255, 203/255, 0.5)),
                     name = "") +
  scale_shape_manual(values=c(regular = "\U25CF",
                              right = "\U25BA",
                              bottom = "\U25BC",
                              top = "\U25B2",
                              `bottom-right` = "\U25E2",
                              `top-right` = "\U25E5"
  ), guide = "none") +
  facet_wrap(~ dataset, ncol = 3, labeller = labeller(dataset = function(labels) as.character(mains[labels]))) +
  xlab("mean of normalized counts") +
  ylab("log fold change") +
  ggtheme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(filename = "analysis/pictures/publication/v2/single_copy_genes_all_datasets.pdf", plot = p, width = 10, height = 7)

##########
##
## Gene Ontology Enrichment tests
##
##########

library(clusterProfiler)
library(org.Mm.eg.db)

go_enrichment <- function(dataset) {

  do_enrichment <- function(gene, universe) {
    require(clusterProfiler)
    if (length(gene) >= 5) {
      ego <- enrichGO(gene          = gene,
                      universe      = universe,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = "MGI",
                      ont           = "ALL",
                      pAdjustMethod = "fdr",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      minGSSize     = 5,
                      maxGSSize     = 500,
                      pool          = TRUE,
                      readable      = TRUE)

      df <- as.data.frame(ego)
      if (nrow(df) == 0) ego <- NULL
    } else {
      ego <- NULL
    }
    return(ego)
  }

  dds <- all_dds[[dataset]]
  universe <- paste0("MGI:", rownames(dds))

  deg_table <- all_deg_tables[[dataset]]

  gene_up <- deg_table %>%
    filter(log2FoldChange > 0) %>%
    mutate(ID = paste0("MGI:", ID)) %>%
    pull(ID)
  ego_up <- do_enrichment(gene_up, universe)

  gene_down <- deg_table %>%
    filter(log2FoldChange < 0) %>%
    mutate(ID = paste0("MGI:", ID)) %>%
    pull(ID)
  ego_down <- do_enrichment(gene_down, universe)

  return(list(upregulated = ego_up,
              downregulated = ego_down))
}

all_goenrichment <- lapply(names(all_dds), go_enrichment) %>% set_names(names(all_dds))

for (i in seq_along(all_goenrichment)){
  results <- all_goenrichment[[i]]
  df <- all_deg_tables[[i]]
  for (direction in c("upregulated", "downregulated")) {
    ego <- results[[direction]]
    tryCatch({
      p <- cnetplot(ego, foldChange = df %>% pull(log2FoldChange, name = Name))
      ggsave(filename = sprintf("analysis/pictures/publication/v2/single_copy_genes_cnet_%s_%s.pdf", names(all_goenrichment)[i], direction),
             plot = p, height = 10, width = 10)
    },
    error = function(er) {})
  }
}



########
##
## Volcano plots
##
######


get_up_down <- function(deseq2_results, pvalue_thr, lfc_thr) {
  x <- deseq2_results$log2FoldChange

  if ("svalue" %in% colnames(deseq2_results)){
    y <- -log10(deseq2_results$svalue)
  } else {
    y <- -log10(deseq2_results$padj)
  }

  k <- y > -log10(pvalue_thr)
  up <- k & x > lfc_thr
  down <- k & x < -lfc_thr
  tot <- up | down
  return(data.frame(row.names = rownames(deseq2_results), total=tot, up=up, down=down))
}

find_shape_volcano <- function(enrichment_score, qt) {
  shape <- "regular"
  if(enrichment_score > qt) {
    shape <- "up"
  }
  return(shape)
}

get_shape_volcano <- function (enrichment_score,
                               qt_thr = 0.98,
                               find_shape_fun = find_shape) {
  qt <- quantile(enrichment_score, qt_thr, na.rm = TRUE)
  shape <- rep("regular", length(enrichment_score))
  for (i in seq_along(enrichment_score)) {
    m <- enrichment_score[i]
    if(!is.na(m)) shape[i] <- find_shape_fun(m, qt)
  }
  return(shape)
}

select_genes_volcano <- function(gene_id, up, down) {
  ret <- rep("gray", length(gene_id))
  ret[up] <- "red"
  ret[down] <- "blue"

  selected_genes <- c(
    Setd1b = "MGI:2652820",
    Kmt2e  = "MGI:1924825",
    Kmt2d  = "MGI:2682319",
    # Brd4   = "MGI:1888520",
    Ep400  = "MGI:1276124"
    # Chd2   = "MGI:2448567",
    # Chd3   = "MGI:1344395",
    # Arid4a = "MGI:2444354",
    # Dnmt3l = "MGI:1859287"
  )
  trim66 <- c(Trim66 = "MGI:2152406")

  for (i in seq_along(selected_genes)) {
    ii <- which(gene_id == selected_genes[i])
    ret[ii] <- names(selected_genes)[i]
  }

  ii <- which(gene_id == trim66)
  ret[ii] <- "Trim66"

  return(ret)
}

get_volcano_data <- function(dds, deseq2_results, pvalue_thr = 0.05, lfc_thr = 1, qt = 1) {

  directions <- get_up_down(deseq2_results, pvalue_thr, lfc_thr)
  up <- directions$up
  down <- directions$down

  x <- deseq2_results$log2FoldChange
  if ("svalue" %in% colnames(deseq2_results)) {
    y <- -log10(deseq2_results$svalue)
  } else {
    y <- -log10(deseq2_results$padj)
  }

  df <- data.frame(gene_id = rownames(deseq2_results), x, y) %>%
    mutate(
      shape = get_shape_volcano(y, qt_thr = qt, find_shape_fun = find_shape_volcano),
      y = if_else(y > quantile(y, qt, na.rm = TRUE), quantile(y, qt, na.rm = TRUE), y),
      color = select_genes_volcano(gene_id, up, down),
      label = sub("gray|red|blue", NA_character_, color)
    )

  return(df)
}

volcano_plot <- function(dds, deseq2_results, title,
                         pvalue_thr = 0.05, lfc_thr = 1,
                         quantile_thr = 1) {

  df <- get_volcano_data(dds, deseq2_results, pvalue_thr, lfc_thr, quantile_thr)

  p <- ggplot(df, aes(x, y, color=color, shape = shape, label = label)) +
    geom_hline(yintercept = -log10(pvalue_thr)) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr)) +
    rasterize(geom_point(data = df %>% filter(!color %in% c("Setd1b", "Kmt2e", "Kmt2d", "Ep400", "Trim66")), size = 3), dpi = 500) +
    rasterize(geom_point(data = df %>% filter(color %in% c("Setd1b", "Kmt2e", "Kmt2d", "Ep400", "Trim66")), size = 4), dpi = 500) +
    geom_label_repel(color = "black") +
    scale_color_manual(values=c(red="firebrick1",
                                blue="cornflowerblue",
                                gray="gray",
                                Trim66="blue",
                                Setd1b = "springgreen2",
                                Kmt2e = "springgreen2",
                                Kmt2d = "springgreen2",
                                Ep400 = "springgreen2"),
                       labels = c(red = "Upregulated",
                                  blue = "Downregulated",
                                  gray = "Gene"),
                       name = "") +
    scale_shape_manual(values=c(regular = "\U25CF",
                                up = "\U25B2"),
                       guide = "none") +
    xlab("logFC") +
    ylab("-log10(p-value)") +
    ggtitle(title) +
    theme_bw()  +
    ggtheme(legend.position = "right")

  return(p)
}

all_qt <-c(TRIM66_Dec2020 = 1,
           TRIM66_Apr2021 = .99994,
           TRIM66_elongatedSpermatids_Ago2021_noeGFP2_noeGFP3_noeGFP6 = 1)

for (dataset in names(all_dds)) {
  dds <- all_dds[[dataset]]
  res <- all_lfcShrink[[dataset]]
  qt <- all_qt[[dataset]]
  p <- volcano_plot(dds, res, dataset, lfc_thr = 0, quantile_thr = qt)
  ggsave(filename = sprintf("analysis/pictures/publication/v2/single_copy_genes_volcano_%s.pdf", dataset),
         plot = p, width = 8, height = 8)
}


# dataset <- names(all_dds)[3]
# dds <- all_dds[[dataset]]
# res <- all_lfcShrink[[dataset]]
# qt <- all_qt[[dataset]]
# volcano_plot(dds, res, dataset, lfc_thr = 0, quantile_thr = 0.99996)



###########
##
## Boxplots of selected H3K4 methyl-transferases
##
###########

library(DESeq2)
library(rstatix)

selected_genes <- c(
  Setd1b  = "MGI:2652820",
  Kmt2e   = "MGI:1924825",
  Kmt2d   = "MGI:2682319",
  Trim66  = "MGI:2152406",
  Brd4    = "MGI:1888520",
  Ep400   = "MGI:1276124",
  Chd2    = "MGI:2448567",
  Chd3    = "MGI:1344395",
  Arid4a  = "MGI:2444354",
  Dnmt3l  = "MGI:1859287",
  Smarca4 = "MGI:88192"
)

dataset <- names(all_dds)[2]
dds1 <- all_dds[[dataset]]
sele <- dds1[selected_genes,]

tpm <- function(counts, len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

mat <- counts(dds1, normalized = FALSE)
w <- width(rowRanges(dds1))
tpm_mat <- tpm(mat, w)
meta <- colData(dds1) %>%
  as_tibble(rownames = "sample") %>%
  dplyr::select(sample, genotype)

gene_names <- tibble(gene_name = names(selected_genes),
  mgi_id = selected_genes)

dat <- tpm_mat[selected_genes,] %>%
  as_tibble(rownames = "mgi_id") %>%
  pivot_longer(cols = -mgi_id, names_to = "sample", values_to = "tpm") %>%
  inner_join(meta) %>%
  inner_join(gene_names)

ypositions <- dat %>%
  group_by(mgi_id) %>%
  summarise(y.position = max(log2(tpm+1)))

pvalues <- all_lfcShrink[[dataset]] %>%
  as_tibble(rownames = "mgi_id") %>%
  filter(mgi_id %in% selected_genes) %>%
  dplyr::select(mgi_id, log2FoldChange, padj) %>%
  inner_join(gene_names) %>%
  mutate(.y. = "tpm",
         group1 = "KO",
         group2 = "WT") %>%
  add_significance("padj") %>%
  inner_join(ypositions)

p <- ggplot(dat, aes(x = gene_name, y = log2(tpm+1))) +
  # geom_violin(mapping = aes(fill = genotype),
  #             trim = FALSE,
  #             position = position_dodge(width = 0.8),
  #             key_glyph = draw_key_rect) +
  geom_boxplot(mapping = aes(color = genotype),
               fill ="#E8E8E8",
               width=0.6,
               position = position_dodge(width = 0.8),
               outlier.color = NA) +
  geom_point(mapping = aes(color = genotype),
             position = position_jitterdodge(jitter.width = 0.4)) +
  geom_text(data = pvalues,
            mapping = aes(x = gene_name,
                          y = max(y.position) + 0.5,
                          label = paste0("P = ", format(padj,digits = 2))),
            size = 5,
            position = position_dodge(width = 1)) +
  geom_text(data = pvalues,
            mapping = aes(x = gene_name,
                          y = max(y.position) + 0.25,
                          label = padj.signif),
            size = 6,
            position = position_dodge(width = 1)) +
  xlab("Gene") +
  ylab("Log2 TPM") +
  scale_colour_discrete(name = "Genotype",
                      type = c(WT = "#6F6A87",
                               KO = "#2FAC66"),
                      labels = c(WT = "Wild type",
                                 KO = bquote("Trim66"^"gfp/gfp"))) +
  ggtheme() +
  theme(aspect.ratio = 0.5,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)); p
ggsave(plot = p, filename = "analysis/pictures/publication/v2/single_copy_genes_selected_histone_mt.pdf",
       height = 8, width = 17)



