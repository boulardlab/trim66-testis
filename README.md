# Trim66’s paternal deficiency causes intrauterine overgrowth
This repository holds the analysis code for the paper "Trim66’s paternal deficiency causes intrauterine overgrowth".

The content is organized as follows:

- `figures/` holds scripts to reproduce most of the figures in the paper.
- `notebooks/` directory holds both the rmarkdown and compiled version of a notebook comparing the RNA-seq datasets.
- `snakemake-chip-seq/` holds a Snakemake pipeline used to process the ChIP-seq data presented in the paper.
- `snakemake-rna-seq/` holds a Snakemake pipeline used to process the RNA-seq data presetend in the paper.

# Running the Snakemake pipelines
To run the pipelines users need to [install Snakemake as per the official docs](https://snakemake.readthedocs.io/en/v7.28.3/getting_started/installation.html#installation-via-conda-mamba). The pipelines were developed and run using Snakemake 7.28.3. Later version should also work.

Users will need to download our raw data from the ENA archive following the link given in the pubblication.

We strongly encourage users to familiarize themselves with Snakemake execution profiles, the Portable Encapsulted Projects standard and Singularity before running the pipelines.

# Reproducing figures
Once the pipelines executed correctly, users can build the notebook in `notebooks/compare_TRIM66_spermatid_libs.Rmd` to generate Figure 4b and 4c, among the others.

In addition, users can run the scripts in the `figures` directory to reproduce other figures:

- `figures/heatmap_markers.R` produces figure S5;
- `figures/plot_sashimi.R` produces figure 2;
- `figures/single_copy_genes_v2.R`: produces figure S6, among other figures;
- `figures/te_v1.R` produces figure 4a;
- `figures/scrna-seq-testis` contains scripts to download and prepare the data for a basic Seurat analysis to generate figure 1B and S1b;
- Figures 4d and 4e are generated from the Snakemake pipeline in `snakemake-chip-seq`;
- `figures/dotplot-npeaks.R` generates figure 4f.

# Reporting issues
If users encouter problems, please use the Issue tracker to open issues, we will get back to you.
