# Trim66’s paternal deficiency causes intrauterine overgrowth
This repository holds the analysis code for the paper [**Trim66’s paternal deficiency causes intrauterine overgrowth**](https://www.biorxiv.org/content/10.1101/2024.02.12.579976v1).

The content is organized as follows:

- `figures/` holds scripts to reproduce most of the figures in the paper.
- `pages/` directory holds files for the [GitHub Pages website](https://boulardlab.github.io/trim66-testis/).
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
- `figures/plot_sashimi_Trim66-PHD.R` produces figure XX. Bam files were generated with the [3t-seq pipeline](https://github.com/boulardlab/3t-seq);
- `figures/single_copy_genes_v2.R`: produces figure S6, among other figures;
- `figures/te_v1.R` produces figure 4a;
- `figures/scrna-seq-testis` contains scripts to download and prepare the data for a basic Seurat analysis to generate figure 1B and S1b;
- Figures 4d and 4e are generated from the Snakemake pipeline in `snakemake-chip-seq`;
- `figures/dotplot-npeaks.R` generates figure 4f.
- `figures/trim66-gfp-fpkm` contains scripts to download a version of the mm10 genome, insert GFP in the correct genomic location, run the STAR aligner on it, count reads mapping to the GFP interval and produce a boxplot comparing the number of reads detected in WT and mutant samples. The aim is quantify reads mapping the GFP sequence, shown in figure S3.

# Reporting issues
If users encouter problems, please use the Issue tracker to open issues, we will get back to you.

# Citing us

```bibtex
@article {Mielnicka2024.02.12.579976,
	author = {Monika Mielnicka and Francesco Tabaro and Rahul Sureka and Basilia Acurzio and Renata Paoletti and Ferdinando Scavizzi and Marcello Raspa and Alvaro H. Crevenna and Karine Lapouge and Kim Remans and Matthieu Boulard},
	title = {Trim66{\textquoteright}s paternal deficiency causes intrauterine overgrowth},
	elocation-id = {2024.02.12.579976},
	year = {2024},
	doi = {10.1101/2024.02.12.579976},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {The tripartite motif-containing protein 66 (TRIM66, also known as TIF1-delta) is a PHD-Bromo containing protein primarily expressed in post-meiotic male germ cells known as spermatids. Biophysical assays showed that TRIM66 PHD-Bromo domain binds to H3 N-terminus only when lysine 4 is unmethylated. We addressed TRIM66{\textquoteright}s role in reproduction by loss-of-function genetics in the mouse. Males homozygous for Trim66-null mutations produced functional spermatozoa. Round spermatids lacking TRIM66 upregulated a network of genes involved in histone acetylation and H3K4 methylation. Profiling of H3K4me3 patterns in the sperm produced by Trim66-null mutant showed minor alterations below statistical significance. Unexpectedly, Trim66-null males, but not females, sired pups overweight at birth, hence revealing that Trim66 mutations cause a paternal effect phenotype.Competing Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2024/02/14/2024.02.12.579976},
	eprint = {https://www.biorxiv.org/content/early/2024/02/14/2024.02.12.579976.full.pdf},
	journal = {bioRxiv}
}
```
