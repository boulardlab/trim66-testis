
def bw(wildcards):
    """
    This function will return the bigwigs associated with a specific 
    normalization method for all samples.
    """
    samples = get_samples(wildcards, as_list = False)
    return expand(
            config["bigwig_dir"].joinpath(f"{wildcards.normalization_method}", "{sample}.bw"),
            sample=samples,
        )

rule compute_matrix:
    input:
        # Please note that the -R and -S options are defined via input files
        bed=config["genome_gtf"],
        bigwig=bw,
    output:
        # Please note that --outFileName, --outFileNameMatrix and --outFileSortedRegions are exclusively defined via output files.
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/computematrix.html.
        matrix_gz=config["tables_dir"].joinpath(
            "deeptools", "{macs2_mode}", "scale-regions", "{normalization_method}.gz"
        ),
        # optional output files
        matrix_tab=config["tables_dir"].joinpath(
            "deeptools", "{macs2_mode}", "scale-regions", "{normalization_method}.tab"
        ),
        matrix_bed=config["tables_dir"].joinpath(
            "deeptools", "{macs2_mode}", "scale-regions", "{normalization_method}.bed"
        ),
    log:
        config["log_dir"].joinpath(
            "deeptools/compute_matrix", "{macs2_mode}", "scale-regions", "{normalization_method}.log"
        ),
    params:
        # required argument, choose "scale-regions" or "reference-point"
        command="scale-regions",
        # optional parameters
        extra=f"--regionBodyLength 500 -b 100 -a 100 --samplesLabel {' '.join(pep.sample_table['sample_name'].tolist())} --verbose",
    threads: 
        8
    resources:
        mem_mb_per_cpu=2000,
        runtime=360
    wrapper:
        "v1.20.0/bio/deeptools/computematrix"


use rule compute_matrix as compute_matrix_summit with:
    input:
        # Please note that the -R and -S options are defined via input files
        bed=config["peaks_dir"].joinpath("{macs2_mode}", "roi.bed"),
        bigwig=bw,
    output:
        # Please note that --outFileName, --outFileNameMatrix and --outFileSortedRegions are exclusively defined via output files.
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/computematrix.html.
        matrix_gz=config["tables_dir"].joinpath(
            "deeptools", "{macs2_mode}", "reference-point", "{normalization_method}.gz"
        ),
        # optional output files
        matrix_tab=config["tables_dir"].joinpath(
            "deeptools", "{macs2_mode}", "reference-point", "{normalization_method}.tab"
        ),
        matrix_bed=config["tables_dir"].joinpath(
            "deeptools", "{macs2_mode}", "reference-point", "{normalization_method}.bed"
        ),
    log:
        config["log_dir"].joinpath(
            "deeptools/compute_matrix", "{macs2_mode}", "reference-point", "{normalization_method}.log"
        ),
    params:
        # required argument, choose "scale-regions" or "reference-point"
        command="reference-point",
        # optional parameters
        extra=f"--referencePoint center -b 1000 -a 1000 --samplesLabel {' '.join(pep.sample_table['sample_name'].tolist())} --verbose",
    retries: 5
    resources:
        mem_mb_per_cpu=200,
        runtime=lambda wildcards, attempt: 60 * attempt


rule plot_heatmap:
    input:
         # matrix file from deepTools computeMatrix tool
        config["tables_dir"].joinpath(
            "deeptools", "{macs2_mode}", "{computeMatrix_mode}", "{normalization_method}.gz"
        ),
    output:
        # Please note that --outFileSortedRegions and --outFileNameMatrix are exclusively defined via output files.
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotheatmap.html.
        heatmap_img=report(
            config["pictures_dir"].joinpath("heatmap", "{macs2_mode}", "{computeMatrix_mode}", "{normalization_method}.pdf"),  # required
            category="Heatmap"
        ),  # required
        # optional output files
        regions=config["pictures_dir"].joinpath("heatmap", "{macs2_mode}", "{computeMatrix_mode}", "{normalization_method}.bed"),
        heatmap_matrix=config["pictures_dir"].joinpath("heatmap", "{macs2_mode}", "{computeMatrix_mode}", "{normalization_method}.tab")
    log:
        config["log_dir"].joinpath("deeptools", "plot_heatmap", "{macs2_mode}", "{computeMatrix_mode}", "{normalization_method}.log"),
    params:
        # optional parameters
        f"--plotType=fill --colorMap Purples"
    retries: 5
    resources:
        runtime=30,
        mem_mb=lambda wildcards, input, attempt: 2000 * attempt + input.size_mb
    wrapper:
        "v1.21.1/bio/deeptools/plotheatmap"