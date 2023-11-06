starTE <- c("data/rna-seq/alignments/starTE//TRIM66_Apr2021/featureCount/random.txt",
            "data/rna-seq/alignments/starTE//TRIM66_Dec2020/featureCount/random.txt",
            "data/rna-seq/alignments/starTE//TRIM66_elongatedSpermatids_Ago2021_noeGFP2_noeGFP3_noeGFP6/featureCount/random.txt")

sample_sheets <- c(TRIM66_Apr2021 = file.path("data/rna-seq", "sample_sheet_trim66-gfp_totalRNA-seq.csv"),
                   TRIM66_Dec2020 = file.path("data/rna-seq", "sample_sheet_trim66-gfp.csv"),
                   TRIM66_elongatedSpermatids_Ago2021_noeGFP2=file.path("data/rna-seq", "sample_sheet_trim66-elongatedSpermatids_Ago2021_noeGFP2.csv"),
                   TRIM66_elongatedSpermatids_Ago2021_noeGFP2_noeGFP3_noeGFP6=file.path("data/rna-seq", "sample_sheet_trim66-elongatedSpermatids_Ago2021_noeGFP2_noeGFP3_noeGFP6.csv"),
                   TRIM66_elongatedSpermatids_Ago2021=file.path("data/rna-seq", "sample_sheet_trim66-elongatedSpermatids_Ago2021.csv"),
                   TRIM66_Oct2020=file.path("data/rna-seq", "sample_sheet.csv"))

sample_name_paterns <- c("tMBAO", "MBAO", "eGFP", "mbaj")

rmsk <- import("data/references/rmsk.gtf.gz")

datasets_starTE <- basename(dirname(dirname(starTE)))

mains_starTE <- mains[match(datasets_starTE, names(mains))]

to_matrix <- function(df, var = "Geneid") {
  rn <- df %>% pull({{ var }})
  mat <- df %>% dplyr::select(-{{ var }}) %>% as.matrix
  rownames(mat) <- rn
  return(mat)
}

find_count_mat <- function(p, get_fn, ...) {
  t <- read_tsv(p, comment = "#", ...)
  for (pat in sample_name_paterns) {
    patt <- paste0("[^t]", pat)
    if(any(grepl(patt, basename(colnames(t))))) {
      break
    }
  }
  return(t %>% get_fn(pat))
}

get_count_mat_starTE <- function(t, sample_name_pattern) t %>%
  rename_with(function(x, pat) {
    bn <- basename(x)
    bn <- sub(pat, "\\1", bn)
    return(bn)
  }, pat = sprintf(".*(%s[0-9]+).*", sample_name_pattern)) %>%
  dplyr::select(Geneid, starts_with(sample_name_pattern)) %>%
  mutate(across(starts_with(sample_name_pattern), as.integer)) %>%
  to_matrix(var = "Geneid")

starTE_count_mats <- starTE %>%
  map(find_count_mat, get_fn = get_count_mat_starTE) %>%
  set_names(basename(dirname(dirname(starTE))))

get_colData <- function(sample_sheet_path){
  colData <- read_csv(sample_sheet_path, col_names = TRUE) %>%
    dplyr::select(name, genotype) %>%
    DataFrame
  rownames(colData) <- colData$name
  colData$genotype <- factor(colData$genotype, levels = c("WT", "KO"))
  return(colData)
}

all_sample_sheets_starTE <- datasets_starTE %>%
  map(function(p) {
    sp <- sample_sheets[grepl(paste0(p,"$"), names(sample_sheets))]
    return(get_colData(sp))
  }) %>%
  set_names(datasets_starTE)

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

volcano_plot <- function(dds, deseq2_results, title, pvalue_thr = 0.05, lfc_thr = 1) {

  directions <- get_up_down(deseq2_results, pvalue_thr, lfc_thr)
  up <- directions$up
  down <- directions$down

  x <- deseq2_results$log2FoldChange
  if ("svalue" %in% colnames(deseq2_results)) {
    y <- -log10(deseq2_results$svalue)
  } else {
    y <- -log10(deseq2_results$padj)
  }

  color <- rep("gray", length(x))
  color[up] <- "red"
  color[down] <- "blue"

  df <- data.frame(x, y, color)

  p <- ggplot(df, aes(x, y, color=color)) +
    geom_hline(yintercept = -log10(pvalue_thr)) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr)) +
    geom_point(size=0.5) +
    scale_color_manual(values=c(red="red", blue="blue", gray="gray"), ) +
    xlab("logFC") + ylab("-log10(p-value)") + ggtitle(title) +
    theme_bw()  +
    ggtheme(legend.position = "none")

  return(p)
}

do_deseq <- function(count_mat, colData, main){
  require(DESeq2)
  dds <- DESeqDataSet(se = SummarizedExperiment(assays = SimpleList(counts = count_mat),
                                                colData = DataFrame(colData)),
                      design = ~ genotype)
  dds <- DESeq(dds)
  # MA <- plotMA(dds, main = main, returnData = TRUE)

  res <- lfcShrink(dds, coef = "genotype_KO_vs_WT")
  volcano <- volcano_plot(dds, res, main, lfc_thr = 0)

  return(list(
    dds=dds,
    res=res,
    volcano=volcano))
}

all_dds <- pmap(list(starTE_count_mats, all_sample_sheets_starTE, mains_starTE), do_deseq)

rmsk_agg <- rmsk %>%
  as_tibble %>%
  dplyr::select(gene_id, family_id, class_id) %>%
  group_by(gene_id, family_id, class_id) %>%
  summarize(n = n()) %>%
  arrange(-n)


baseMeans <- all_dds %>%
  map(`[[`, "res") %>%
  map(as_tibble, rownames = "gene_id") %>%
  map2(.y = names(all_dds), ~mutate(.x, dataset = .y))  %>%
  map(dplyr::select, gene_id, baseMean, dataset) %>%
  map(inner_join, rmsk_agg)

q <- seq(0, 1, by = .1)

quantile_values <- baseMeans %>%
  purrr::reduce(bind_rows) %>%
  dplyr::filter(baseMean > 0) %>%
  group_by(dataset) %>%
  summarize(qtls = paste(quantile(baseMean, probs = q), collapse = ",")) %>%   # get quantile values as a string
  separate(qtls, paste0(100*q,  "%"), sep = ",", convert = T) %>%
  group_split(dataset)

baseMeans %>%
  purrr::reduce(bind_rows) %>%
  mutate(log2BaseMean = log2(baseMean)) %>%
  ggplot(aes(dataset, log2BaseMean)) +
  geom_violin(draw_quantiles = q, trim = FALSE) +
  geom_jitter(width = 0.15, color = rgb(0,0,0,0.2)) +
  # scale_x_discrete(labels = function(x) mains[x]) +
  # facet_wrap(~ class_id) +
  ggtheme()

qt_thr <- 0.98

plot_dat <- map(names(all_dds), function(dataset) all_dds[[dataset]][["res"]]) %>%
  map(as_tibble, rownames = "gene_id") %>%
  set_names(names(all_dds)) %>%
  # map(dplyr::select, baseMean, log2FoldChange, padj) %>%
  map2(.y = names(all_dds),
       ~mutate(.x, dataset = .y)
  ) %>%
  map2(.y = quantile_values,
       ~mutate(.x, baseMeanOverQuantile = baseMean > .y$`40%`)
  ) %>%
  purrr::reduce(bind_rows) %>%
  mutate(
    shape = get_shape(baseMean, log2FoldChange),
    color = if_else(!baseMeanOverQuantile, "low expression", if_else(padj < 1e-2 & log2FoldChange > 1, "DE", "none")),
    baseMean = if_else(baseMean > quantile(baseMean, qt_thr), quantile(baseMean, qt_thr), baseMean),
    log2FoldChange = if_else(log2FoldChange < -2, -2, if_else(log2FoldChange > 2, 2, log2FoldChange))
  )

plot_dat %>%
  ggplot() +
  rasterize(geom_point(data = . %>% dplyr::filter(color != "DE"),
             aes(baseMean, log2FoldChange, color = color , shape = shape), size = 3), dpi = 500) +
  rasterize(geom_point(data = . %>% dplyr::filter(color == "DE"),
             aes(baseMean, log2FoldChange, color = color , shape = shape), size = 3), dpi = 500) +
  geom_hline(yintercept = 0, color = "black") +
  scale_shape_manual(values=c(regular = "\U25CF",
                              right = "\U25BA",
                              bottom = "\U25BC",
                              top = "\U25B2",
                              `bottom-right` = "\U25E2",
                              `top-right` = "\U25E5"),
                     guide = NULL) +
  scale_color_manual(values = c(DE = "red", none = rgb(203/255, 203/255, 203/255), `low expression` = "gray90"),
                     name = "Differential expression") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  facet_wrap(~ dataset, ncol = 3,
             labeller = labeller(dataset = function(labels) as.character(mains[labels]))) +
  xlab("mean of normalized counts") +
  ylab("log fold change") +
  ylim(-2,2) +
  ggtheme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("analysis/pictures/publication/v1/te_all_datasets.pdf", width = 10, height = 7)


# FPKM
rr <- as.data.frame(rmsk) %>%
  group_by(gene_id) %>%
  nest() %>%
  mutate(data = map(data, makeGRangesFromDataFrame, keep.extra.columns = TRUE)) %>%
  pull(data, name = "gene_id") %>%
  GRangesList

all_fpkm <- all_dds %>%
  map(function(x) x[["dds"]]) %>%
  map(function(x) { rowRanges(x) <- rr; x }) %>%
  map(fpkm) %>%
  set_names(names(all_dds))

fpkm_wide <- all_fpkm %>%
  map(as_tibble, rownames = "gene_id") %>%
  map(pivot_longer, cols = -gene_id) %>%
  map2(.y = names(all_dds), ~mutate(.x, dataset = .y)) %>%
  map2(.y = all_sample_sheets_starTE, inner_join, copy = T, by = "name") %>%
  purrr::reduce(bind_rows) %>%
  dplyr::filter(dataset != "TRIM66_Dec2020") %>%
  group_by(gene_id, genotype, dataset) %>%
  summarize(median_fpkm = log2(median(value) + 1)) %>%
  pivot_wider(id_cols = c(gene_id, dataset), names_from = genotype, values_from = median_fpkm)

cors <- fpkm_wide %>%
  group_by(dataset) %>%
  summarize(pearson = cor(WT, KO)) %>%
  mutate(label = sprintf("Pearson = %.3f", pearson))

targets <- c("L1Md_A", "L1Md_T", "ERVB4_2-LTR_MM", "MERVL-int")
colors <- c("L1Md_A" = "red",
            "L1Md_T" = "green",
            "ERVB4_2-LTR_MM" = "lightsteelblue",
            "MERVL-int" = "orange")

library(ggrepel)
library(ggnewscale)
fpkm_wide <- fpkm_wide %>%
  mutate(te_label = if_else(gene_id %in% targets, gene_id, NA_character_),
         col = if_else(gene_id %in% targets, colors[gene_id], "black"))

ggplot() +
  rasterise(geom_point(aes(x = WT, y = KO), data = fpkm_wide, size = 0.8),
            dpi = 500) +
  geom_abline(slope = 1, intercept = 0, color = "dodgerblue") +

  new_scale_color() +
  geom_point(aes(x = WT, y = KO, color = col),
             data = fpkm_wide %>% filter(!is.na(te_label)),
             size = 2) +
  scale_color_identity(guide = "legend",
                       name = "",
                       labels = c(setNames(names(colors), colors),
                                  black = "TE gene")) +

  # geom_label_repel(aes(x = WT, y = KO, label = te_label),
  #                  data = fpkm_wide, fill = "white", force = 2) +

  geom_text(data = cors, aes(label = label), x = 10, y = 1, size = 5) +

  facet_wrap(~ dataset, labeller = labeller(dataset = function(x) mains[x])) +

  xlab("log2 median FPKM WT") + ylab("log2 median FPKM KO") +
  ggtheme(legend.position = "bottom") +
  theme(panel.grid = element_blank())

ggsave("analysis/pictures/publication/v2/te_scatter_median_fpkm.pdf",
       width = 10, height = 7)
