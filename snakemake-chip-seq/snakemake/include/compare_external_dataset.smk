def get_external_dataset(wildcards):
    selected_dataset = list(
        filter(
            lambda x: x["name"] == wildcards.external_dataset,
            list(config["external_datasets"]),
        )
    )[0]
    folder = selected_dataset["bigwig"]
    if folder.endswith("/"):
        folder = folder[:-1]
    ls = glob.glob(f"{folder}/*.bw") + glob.glob(f"{folder}/*.bigwig")
    return ls


use rule screenshot_config as screenshot_config_bw_compare with:
    input:
        bigwig=get_bw_by_norm,
        bed=get_sloped_peaks_bed,
        calls=config["peaks_dir"].joinpath(
            "{macs2_mode}", "intervene", "sets", "{set}.bed"
        ),
        genome=config["genome_fasta_filter_index"],
        external_dataset=get_external_dataset,
    output:
        batch=config["peaks_dir"].joinpath(
            "{external_dataset}",
            "{macs2_mode}",
            "{normalization_method}",
            "{set}.igv.batch",
        ),
    params:
        genome_label=config["genome"]["label"],
        snapshot_directory=lambda wildcards: config["pictures_dir"].joinpath(
            "peaks",
            wildcards.external_dataset,
            wildcards.macs2_mode,
            wildcards.normalization_method,
        ),
        sample_table=pep.sample_table.to_dict(),
    log:
        config["log_dir"].joinpath(
            "intervene",
            "screenshot_config",
            "{macs2_mode}",
            "{normalization_method}",
            "{external_dataset}",
            "{set}.log",
        ),


use rule igv_screenshot as igv_screenshot_external_dataset with:
    wildcard_constraints:
        normalization_method="|".join(normalization_methods),
        external_dataset="|".join(external_datasets),
    input:
        batchfile=config["peaks_dir"].joinpath(
            "{external_dataset}",
            "{macs2_mode}",
            "{normalization_method}",
            "{set}.igv.batch",
        ),
    output:
        temp(
            touch(
                config["pictures_dir"].joinpath(
                    "peaks",
                    "{external_dataset}",
                    "{macs2_mode}",
                    "{normalization_method}",
                    "{set}",
                    "igv.done",
                )
            )
        ),
    log:
        config["log_dir"].joinpath(
            "igv/{macs2_mode}/{normalization_method}/{external_dataset}/{set}.log"
        ),


def get_all_external_igv_flags(wildcards):
    intervene_venn_output = checkpoints.intervene_venn.get(**wildcards).output[2]
    return expand(
        config["pictures_dir"].joinpath(
            "peaks",
            wildcards.external_dataset,
            wildcards.macs2_mode,
            wildcards.normalization_method,
            "{intervene_set}",
            "igv.done",
        ),
        intervene_set=glob_wildcards(
            os.path.join(intervene_venn_output, "{intervene_set}.bed")
        ).intervene_set,
    )


use rule flag_igv as flag_igv_external with:
    wildcard_constraints:
        normalization_method="|".join(normalization_methods),
        external_dataset="|".join(external_datasets),
    input:
        get_all_external_igv_flags,
    output:
        touch(
            config["pictures_dir"].joinpath(
                "peaks",
                "{external_dataset}",
                "{macs2_mode}",
                "{normalization_method}",
                "igv.done",
            )
        ),


def get_path_to_external_dataset(wildcards):
    selected_dataset = list(
        filter(
            lambda x: x["name"] == wildcards.external_dataset,
            list(config["external_datasets"]),
        )
    )[0]
    folder = selected_dataset["bigwig"]
    ret = []
    if os.path.exists(os.path.join(folder, f"{wildcards.sample}.bw")):
        ret = os.path.join(folder, f"{wildcards.sample}.bw")
    elif os.path.exists(os.path.join(folder, f"{wildcards.sample}.bigwig")):
        ret = os.path.join(folder, f"{wildcards.sample}.bigwig")
    else:
        raise ValueError(f"{folder} contains unrecognized files.")
    return ret


rule bigwig2bedgraph:
    input:
        get_path_to_external_dataset,
    output:
        temp(
            config["results_dir"].joinpath(
                "external-comparisons", "{external_dataset}", "bdg", "{sample}.bdg"
            )
        ),
    conda:
        "../env/ucsc-bigwigtobedgraph.yml"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "bigWigToBedGraph {input} {output}"


use rule bigwig2bedgraph as bigwig2bedgraph_internal with:
    wildcard_constraints:
        normalization_method="|".join(normalization_methods),
    input:
        config["bigwig_dir"].joinpath("{normalization_method}", "{sample}.bw"),
    output:
        temp(config["bedgraph_dir"].joinpath("{normalization_method}", "{sample}.bdg")),


rule bedtools_sort_bdg:
    wildcard_constraints:
        normalization_method="|".join(normalization_methods),
    input:
        bdg=config["bedgraph_dir"].joinpath(
            "{normalization_method}", "{internal_sample}.bdg"
        ),
        genome=config["genome_fasta_filter_index"],
    output:
        temp(
            config["bedgraph_dir"].joinpath(
                "sorted", "{normalization_method}", "{internal_sample}.bdg"
            )
        ),
    conda:
        "../env/bedtools.yml"
    threads: 1
    retries: 3
    resources:
        mem_mb=lambda wildcards, input, attempt: 4000 + input.size_mb
        if attempt == 1
        else (4**attempt) * 1000 + input.size_mb,
        runtime=lambda wildcardsds, attempt: 90
        if attempt == 1
        else 90 + (3**attempt) * 10,
    log:
        config["log_dir"].joinpath(
            "external-comparison",
            "bedtools-sort",
            "{normalization_method}",
            "{internal_sample}.log",
        ),
    script:
        "../scripts/bedtools-sort.sh"


use rule bedtools_sort_bdg as bedtools_sort_bdg_external with:
    input:
        bdg=config["results_dir"].joinpath(
            "external-comparisons",
            "{external_dataset}",
            "bdg",
            "filter",
            "{external_sample}.bdg",
        ),
        genome=config["genome_fasta_filter_index"],
    output:
        temp(
            config["results_dir"].joinpath(
                "external-comparisons",
                "{external_dataset}",
                "bdg",
                "sorted",
                "{external_sample}.bdg",
            )
        ),
    log:
        config["log_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "bedtools-sort",
            "{external_sample}.log",
        ),


rule filter_external_bdg_by_shared_chromosomes:
    input:
        genome=config["genome_fasta_filter_index"],
        b=config["results_dir"].joinpath(
            "external-comparisons",
            "{external_dataset}",
            "bdg",
            "{external_sample}.bdg",
        ),
    output:
        config["results_dir"].joinpath(
            "external-comparisons",
            "{external_dataset}",
            "bdg",
            "filter",
            "{external_sample}.bdg",
        ),
    conda:
        "../env/bedtools.yml"
    threads: 1
    retries: 3
    resources:
        mem_mb=lambda wildcards, input, attempt: 8000
        if attempt == 1
        else (4**attempt) * 1000,
        runtime=120,
    log:
        config["log_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "filter-by-chromosome",
            "{external_sample}.log",
        ),
    script:
        "../scripts/match-chromosome-names.sh"


def get_tile_size(wildcards):
    tile_size = 1000 # default to 1kbp
    if len(config["external_datasets"]) > 0:
        external_dataset = wildcards.external_dataset
        for dataset in config["external_datasets"]:
            if dataset["name"] == external_dataset:
                tile_size = dataset["tile-size"]
    return tile_size

rule make_tiles:
    input:
        config["genome_fasta_filter_index"],
    output:
        config["results_dir"].joinpath(
            "external-comparisons",
            "{external_dataset}",
            "correlation-analysis",
            "tiles.bed",
        ),
    log:
        config["log_dir"].joinpath(
            "external-comparison", "{external_dataset}", "make-tiles.log"
        ),
    params:
        window_size=lambda wildcards: get_tile_size(wildcards),
    conda:
        "../env/bedtools.yml"
    threads: 1
    resources:
        runtime=30,
        mem_mb=4000,
    script:
        "../scripts/bedtools-makewindows.sh"


rule extract_bigwig_signals:
    wildcard_constraints:
        normalization_method="|".join(normalization_methods),
    input:
        sample=config["bedgraph_dir"].joinpath(
            "sorted", "{normalization_method}", "{internal_sample}.bdg"
        ),
        tiles=config["results_dir"].joinpath(
            "external-comparisons",
            "{external_dataset}",
            "correlation-analysis",
            "tiles.bed",
        ),
        genome=config["genome_fasta_filter_index"],
    output:
        config["results_dir"].joinpath(
            "external-comparisons",
            "{external_dataset}",
            "correlation-analysis",
            "{normalization_method}",
            "{internal_sample}.bed",
        ),
    params:
        fun="max"
    log:
        config["log_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "bedtools-map",
            "{normalization_method}",
            "{internal_sample}.log",
        ),
    conda:
        "../env/bedtools.yml"
    threads: 1
    resources:
        runtime=30,
        mem_mb=8000,
    script:
        "../scripts/bedtools-map.sh"


use rule extract_bigwig_signals as extract_bigwig_signals_external with:
    input:
        sample=config["results_dir"].joinpath(
            "external-comparisons",
            "{external_dataset}",
            "bdg",
            "filter",
            "{external_sample}.bdg",
        ),
        tiles=config["results_dir"].joinpath(
            "external-comparisons",
            "{external_dataset}",
            "correlation-analysis",
            "tiles.bed",
        ),
        genome=config["genome_fasta_filter_index"],
    output:
        config["results_dir"].joinpath(
            "external-comparisons",
            "{external_dataset}",
            "correlation-analysis",
            "{external_dataset}",
            "{external_sample}.bed",
        ),
    log:
        config["log_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "bedtools-map",
            "{external_dataset}",
            "{external_sample}.log",
        ),


rule join_bedgraph:
    wildcard_constraints:
        normalization_method="|".join(normalization_methods),
        external_dataset="|".join(get_external_datasets()),
    input:
        a=config["results_dir"].joinpath(
            "external-comparisons",
            "{external_dataset}",
            "correlation-analysis",
            "{normalization_method}",
            "{internal_sample}.bed",
        ),
        b=config["results_dir"].joinpath(
            "external-comparisons",
            "{external_dataset}",
            "correlation-analysis",
            "{external_dataset}",
            "{external_sample}.bed",
        )
    output:
        config["results_dir"].joinpath(
            "external-comparisons",
            "{external_dataset}",
            "correlation-analysis",
            "{normalization_method}",
            "{internal_sample}",
            "{external_sample}.txt",
        ),
    log:
        config["log_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "correlation-analysis",
            "{normalization_method}",
            "join-bedgraph",
            "{internal_sample}",
            "{external_sample}.log",
        ),
    threads: 1
    retries: 2
    resources:
        runtime=lambda wildcards, attempt: 90 * attempt,
        mem_mb=lambda wildcards, input, attempt: input.size_mb * 2
        + (2**attempt) * 1000,
    conda:
        "../env/tidyverse.yml"
    script:
        "../scripts/join-bedgraph.sh"


rule compare_signals:
    wildcard_constraints:
        normalization_method="|".join(normalization_methods),
        external_dataset="|".join(get_external_datasets()),
    input:
        config["results_dir"].joinpath(
            "external-comparisons",
            "{external_dataset}",
            "correlation-analysis",
            "{normalization_method}",
            "{internal_sample}",
            "{external_sample}.txt"
        ),
    output:
        corr=config["results_dir"].joinpath(
            "external-comparisons",
            "{external_dataset}",
            "correlation-analysis",
            "{normalization_method}",
            "{internal_sample}",
            "{external_sample}.corr",
        ),
        plot=config["results_dir"].joinpath(
            "external-comparisons",
            "{external_dataset}",
            "correlation-analysis",
            "{normalization_method}",
            "{internal_sample}",
            "{external_sample}.pdf",
        ),
    log:
        config["log_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "correlation-analysis",
            "{normalization_method}",
            "{internal_sample}",
            "{external_sample}.log",
        ),
    threads: 1
    retries: 1
    resources:
        runtime=lambda wildcards, attempt: 180 * attempt,
        mem_mb=lambda wildcards, input, attempt: input.size_mb * 2
        + 3 * attempt * 1000,
        disk_mb=lambda wildcards, input, attempt: input.size_mb * 2
        + attempt * 1000,
    conda:
        "../env/tidyverse.yml"
    script:
        "../scripts/compare-signals.R"

def get_all_external_sample_names(wildcards):
    # get external dataset name
    external_dataset = wildcards.external_dataset

    # get external samples
    external_samples = []
    for dataset in config["external_datasets"]:
        if dataset["name"] == external_dataset:
            folder = dataset["bigwig"]

            if folder.endswith("/"):
                folder = folder[:-1]

            for ext in [".bw", ".bigwig"]:
                for filename in glob.glob(f"{folder}/*{ext}"):
                    filename = os.path.basename(filename)
                    if filename.endswith(ext):
                        external_samples.append(filename.replace(ext, ""))
    
    return external_samples

def get_all_corr_files(wildcards):
    # get all samples names
    internal_samples = get_samples(wildcards)

    # get current normalization method
    normalization_method = wildcards.normalization_method

    external_samples = get_all_external_sample_names(wildcards)
    # build paths to the txt files from compare_signals
    ret = []
    for internal_sample in internal_samples:
        for external_sample in external_samples:
            ret.append(
                config["results_dir"].joinpath(
                    "external-comparisons",
                    wildcards.external_dataset,
                    "correlation-analysis",
                    normalization_method,
                    internal_sample,
                    f"{external_sample}.corr",
                )
            )
    return ret


rule plot_correlation_heatmap:
    input:
        get_all_corr_files,
    output:        
        config["results_dir"].joinpath(
            "external-comparisons",
            "{external_dataset}",
            "correlation-analysis",
            "{normalization_method}",
            "heatmap.pdf",
        )
    log:
        config["log_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "correlation-analysis",
            "{normalization_method}",
            "heatmap.log",
        ),
    threads: 1
    retries: 3
    resources:
        runtime=lambda wildcards, attempt: 40 * attempt,
        mem_mb=lambda wildcards, attempt: (4**attempt) * 1000,
    conda:
        "../env/pheatmap.yml"
    script:
        "../scripts/plot-correlation-heatmap.R"



def get_all_external_bigwig(wildcards):
    # get external dataset name
    external_dataset = wildcards.external_dataset

    # get external samples
    external_samples = []
    for dataset in config["external_datasets"]:
        if dataset["name"] == external_dataset:
            folder = dataset["bigwig"]

            if folder.endswith("/"):
                folder = folder[:-1]

            for ext in [".bw", ".bigwig"]:
                external_samples += glob.glob(f"{folder}/*{ext}")
                
    return external_samples


rule multiBigwigSummary:
    input:
        bigwig=get_all_external_bigwig,
        blacklist=config["encode_blacklist"]
    output:
        config["results_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "summary.npz",
        )
    params:
        labels=lambda wildcards: " ".join(get_all_external_sample_names(wildcards)),
        binsize=lambda wildcards: get_tile_size(wildcards),
    log:
        config["log_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "multibigwigsummary.log",
        ),
    threads: 8
    retries: 3
    resources:
        runtime=lambda wildcards, attempt: 40 * attempt,
        mem_mb=lambda wildcards, attempt: (4**attempt) * 1000,
    conda:
        "../env/deeptools.yml"
    shell:
        """
        multiBigwigSummary bins \
        -b {input.bigwig} \
        -o {output} \
        -l {params.labels} \
        -bs {params.binsize} \
        -bl {input.blacklist} \
        -p {threads}
        """
    
use rule plot_correlation_heatmap_with_summary as plot_correlation_heatmap_with_summary_external_dataset with:
    input:
        config["log_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "summary.npz",
        )
    output:
        config["results_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "correlation-heatmap.pdf",
        ),
        config["results_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "correlation-heatmap.txt",
        )
    log:
        config["log_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "plotHeatmap.log",
        ),
    params:
        config["plotCorrelation"]["params"],
        labels=lambda wildcards: " ".join(get_all_external_sample_names(wildcards)),


use rule plot_pca as plot_pca_external_dataset with:
    input:
        config["results_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "summary.npz",
        )
    output:
        config["results_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "correlation-pca.pdf",
        ),
        config["results_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "correlation-pca.txt",
        )
    log:
        config["log_dir"].joinpath(
            "external-comparison",
            "{external_dataset}",
            "plotPCA.log",
        ),
    params:
        labels=lambda wildcards: " ".join(get_all_external_sample_names(wildcards)),

