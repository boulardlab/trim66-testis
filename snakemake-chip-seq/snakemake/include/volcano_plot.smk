def quantify_roi_input(wildcards):
    d = get_filtered_bams(wildcards)
    d["roi"] = config["peaks_dir"].joinpath(wildcards.macs2_mode, "roi.bed")
    return d


rule quantify_roi:
    input:
        unpack(quantify_roi_input),
        # roi=config["peaks_dir"].joinpath("{macs2_mode}", "roi.bed"),
        # bam=config["alignments_dir"].joinpath("filtered", "sorted", "{sample}.bam"),
        # bai=config["alignments_dir"].joinpath("filtered", "sorted", "{sample}.bam.bai"),
        # genome=config["genome_fasta_filter_index"],
    output:
        config["coverage_dir"].joinpath("{macs2_mode}", "coverage.txt"),
    conda:
        "../env/bedtools.yml"
    params:
        mapq=20,
    log:
        config["log_dir"].joinpath("volcano_plot/{macs2_mode}/bedtools-multicov.log"),
    resources:
        runtime=10,
        mem_mb=lambda wildcards, input: 1500 + input.size_mb,
    script:
        "../scripts/bedtools-coverage.sh"


# def get_all_coverage(wildcards):
#     samples = get_samples(wildcards, as_list=True)
#     return expand(
#         config["coverage_dir"].joinpath("{{macs2_mode}}", "{sample}.bed"),
#         sample=samples,
#     )


# rule build_coverage_matrix:
#     input:
#         # get_all_coverage,
#     output:
#         config["coverage_dir"].joinpath("{macs2_mode}", "coverage.txt"),
#     params:
#         sample_sheet=pep.sample_table.to_dict(orient="list"),
#     conda:
#         "../env/tidyverse.yml"
#     resources:
#         runtime=30,
#         mem_mb=lambda wildcards, input: 1500 + input.size_mb,
#     log:
#         config["log_dir"].joinpath(
#             "volcano_plot/{macs2_mode}/build_coverage_matrix.log"
#         ),
#     script:
#         "../scripts/build_coverage_matrix.R"


rule get_num_mappd_reads:
    input:
        config["log_dir"].joinpath("bowtie2/{sample}.log"),
    output:
        config["alignments_dir"].joinpath("raw", "{sample}.total-aligned.txt"),
    conda:
        "../env/bash.yml"
    resources:
        runtime=10,
        mem_mb=2000,
    shell:
        """
        cat {input} | \
        sed -n '4,5'p | \
        awk '{{ print $1 }}' | \
        tr "\\n" "+" | \
        rev | \
        cut -c 2- | \
        rev | \
        bc > {output}
        """


def get_all_fragsizes(wildcards):
    samples = get_samples(wildcards, as_list=True)
    return expand(
        config["alignments_dir"].joinpath(
            "filtered", "sorted", "{sample}.fragsize.txt"
        ),
        sample=samples,
    )


def get_all_numbers_mapped_reads(wildcards):
    samples = get_samples(wildcards, as_list=True)
    return expand(
        config["alignments_dir"].joinpath("raw", "{sample}.total-aligned.txt"),
        sample=samples,
    )


def get_all_effective_genome_sizes():
    read_length = pep.sample_table["read_length"].tolist()
    genome_label = config["genome"]["label"]
    ret = [effective_genome_size[genome_label][rl] for rl in read_length]
    return ret


rule normalize_matrix_se:
    input:
        matrix=config["coverage_dir"].joinpath("{macs2_mode}", "coverage.txt"),
        fragment_size=get_all_fragsizes,
        number_of_mapped_reads=get_all_numbers_mapped_reads,
    output:
        config["coverage_dir"].joinpath("{macs2_mode}", "coverage-norm.txt"),
    conda:
        "../env/tidyverse.yml"
    params:
        effective_genome_size=get_all_effective_genome_sizes(),
        sample_names=pep.sample_table["sample_name"].tolist(),
    resources:
        runtime=30,
        mem_mb=lambda wildcards, input: 1500 + input.size_mb,
    log:
        config["log_dir"].joinpath(
            "volcano_plot/{macs2_mode}/normalize_coverage_matrix.log"
        ),
    script:
        "../scripts/normalize_coverage_matrix.R"


rule calculate_logfold_and_test:
    input:
        config["coverage_dir"].joinpath("{macs2_mode}", "coverage-norm.txt"),
    output:
        config["coverage_dir"].joinpath("{macs2_mode}", "coverage-lfc-pvalues.txt"),
    conda:
        "../env/tidyverse.yml"
    params:
        sample_sheet=pep.sample_table.to_dict(orient="list"),
        contrast_column="genotype",
        contrast_levels=["KO", "WT"],
    resources:
        runtime=30,
        mem_mb=lambda wildcards, input: 3000 + input.size_mb,
    log:
        config["log_dir"].joinpath(
            "volcano_plot/{macs2_mode}/calculate_log2_fold_change.log"
        ),
    script:
        "../scripts/calculate_log2_fold_change.R"


rule plot_volcano:
    input:
        config["coverage_dir"].joinpath("{macs2_mode}", "coverage-lfc-pvalues.txt"),
    output:
        report(
            config["pictures_dir"].joinpath(
                "differential-binding", "{macs2_mode}", "volcano.pdf"
            ),
            category="Differential binding",
        ),
    conda:
        "../env/tidyverse.yml"
    resources:
        runtime=30,
        mem_mb=lambda wildcards, input: 3000 + input.size_mb,
    log:
        config["log_dir"].joinpath("differential-binding/{macs2_mode}/plot_volcano.log"),
    script:
        "../scripts/plot_volcano.R"


rule plot_scatterplot:
    input:
        config["coverage_dir"].joinpath("{macs2_mode}", "coverage-norm.txt"),
    output:
        report(
            config["pictures_dir"].joinpath(
                "differential-binding", "{macs2_mode}", "scatterplot.pdf"
            ),
            category="Differential binding",
        ),
    conda:
        "../env/tidyverse.yml"
    params:
        sample_sheet=pep.sample_table.to_dict(orient="list"),
        contrast_column="genotype",
        contrast_levels=["KO", "WT"],
    resources:
        runtime=10,
        mem_mb=lambda wildcards, input: 3000 + input.size_mb,
    log:
        config["log_dir"].joinpath(
            "differential-binding/{macs2_mode}/plot_differential_peaks_scatterplot.log"
        ),
    script:
        "../scripts/plot_differential_peaks_scatterplot.R"


use rule compute_matrix as compute_matrix_roi with:
    input:
        bed=config["peaks_dir"].joinpath("{macs2_mode}", "roi.bed"),
        bigwig=get_bw_by_norm,
    output:
        # Please note that --outFileName, --outFileNameMatrix and --outFileSortedRegions are exclusively defined via output files.
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/computematrix.html.
        matrix_gz=config["tables_dir"].joinpath(
            "differential-binding",
            "{macs2_mode}",
            "scale-regions",
            "{normalization_method}.gz",
        ),
        # optional output files
        matrix_tab=config["tables_dir"].joinpath(
            "differential-binding",
            "{macs2_mode}",
            "scale-regions",
            "{normalization_method}.tab",
        ),
        matrix_bed=config["tables_dir"].joinpath(
            "differential-binding",
            "{macs2_mode}",
            "scale-regions",
            "{normalization_method}.bed",
        ),
    log:
        config["log_dir"].joinpath(
            "differential-binding/compute_matrix",
            "{macs2_mode}",
            "scale-regions",
            "{normalization_method}.log",
        ),
    params:
        # required argument, choose "scale-regions" or "reference-point"
        command="scale-regions",
        # optional parameters
        extra=f"-bs 1 -a 0 -b 0 --samplesLabel {' '.join(pep.sample_table['sample_name'].tolist())} --verbose",
    retries: 5
    resources:
        mem_mb_per_cpu=200,
        runtime=lambda wildcards, attempt: 60 * attempt,


use rule plot_heatmap as plot_heatmap_roi with:
    input:
        config["tables_dir"].joinpath(
            "differential-binding",
            "{macs2_mode}",
            "scale-regions",
            "{normalization_method}.gz",
        ),
    output:
        heatmap_img=config["pictures_dir"].joinpath(
            "heatmap-roi",
            "{macs2_mode}",
            "{computeMatrix_mode}",
            "{normalization_method}.pdf",
        ),
        # required
        regions=config["pictures_dir"].joinpath(
            "heatmap",
            "{macs2_mode}",
            "{computeMatrix_mode}",
            "{normalization_method}.bed",
        ),
        heatmap_matrix=config["pictures_dir"].joinpath(
            "heatmap",
            "{macs2_mode}",
            "{computeMatrix_mode}",
            "{normalization_method}.tab",
        ),
    log:
        config["log_dir"].joinpath(
            "differential-binding",
            "plot_heatmap",
            "{macs2_mode}",
            "{computeMatrix_mode}",
            "{normalization_method}.log",
        ),
    params:
        # optional parameters
        f"--plotType=fill --colorMap Purples",
