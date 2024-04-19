def get_all_homer(wildcards):
    samples = get_samples(wildcards, as_list=False)
    return expand(config["peaks_dir"].joinpath(
            wildcards.macs2_mode, 
            "{sample}", 
            "{sample}_homer.txt"),
            sample=samples)


rule get_promoter_peaks_bed:
    input: 
        annotations=get_all_homer
    output:
        config["peaks_dir"].joinpath(
            "{macs2_mode}",
            "promoter_peaks.bed")
    conda: 
        "../env/bash.yml"
    resources:
        runtime=5,
        mem_mb=500
    threads: 1
    script:
        "../scripts/get_promoter_peaks_bed.sh"

rule get_promoter_peaks_gtf:
    input:
        annotations=get_all_homer,
        gtf=str(config["genome_gtf"]).replace(".gz", "")
    output: 
         config["peaks_dir"].joinpath(
            "{macs2_mode}",
            "promoter_peaks.gtf")
    conda:
        "../env/bash.yml"
    resources: 
        runtime=10
    script:
        "../scripts/get_promoter_peaks_gtf.sh"

rule compute_matrix_promoters:
    input:
        # Please note that the -R and -S options are defined via input files
        bed=config["peaks_dir"].joinpath(
            "{macs2_mode}",
            "promoter_peaks.bed"),
        bigwig=bw,
    output:
        # Please note that --outFileName, --outFileNameMatrix and --outFileSortedRegions are exclusively defined via output files.
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/computematrix.html.
        matrix_gz=config["tables_dir"].joinpath(
            "deeptools-promoters", "{macs2_mode}", "reference-point", "{normalization_method}.gz"
        ),
        # optional output files
        matrix_tab=config["tables_dir"].joinpath(
            "deeptools-promoters", "{macs2_mode}", "reference-point", "{normalization_method}.tab"
        ),
        matrix_bed=config["tables_dir"].joinpath(
            "deeptools-promoters", "{macs2_mode}", "reference-point", "{normalization_method}.bed"
        ),
    log:
        config["log_dir"].joinpath(
            "deeptools-promoters", "compute_matrix", "{macs2_mode}", "reference-point", "{normalization_method}.log"
        ),
    params:
        # required argument, choose "scale-regions" or "reference-point"
        command="reference-point",
        # optional parameters
        extra=f"--referencePoint center -b 500 -a 500 --samplesLabel {' '.join(pep.sample_table['sample_name'].tolist())} --verbose",
    retries: 5
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
        runtime=lambda wildcards, attempt: 180 * attempt
    conda:
        "../env/deeptools.yml"
    shell:
        """
        if [ $(file -b {input.bed}) = "empty" ]
        then
            touch {output.matrix_gz}
            touch {output.matrix_tab}
            touch {output.matrix_bed}
            echo "No promoter peaks detected" | tee {log} > $(dirname {output.matrix_gz})/README
        else
            computeMatrix \
            {params.command} \
            {params.extra} \
            --numberOfProcessors {threads} \
            -R {input.bed} \
            -S {input.bigwig} \
            -o {output.matrix_gz} \
            --outFileNameMatrix {output.matrix_tab} \
            --outFileSortedRegions {output.matrix_bed} |& tee {log}
        fi
        """


use rule compute_matrix_promoters as compute_matrix_promoters_gene_bodies with:
    input:
        # Please note that the -R and -S options are defined via input files
        bed=config["peaks_dir"].joinpath(
            "{macs2_mode}",
            "promoter_peaks.gtf"),
        bigwig=bw,
    output:
        # Please note that --outFileName, --outFileNameMatrix and --outFileSortedRegions are exclusively defined via output files.
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/computematrix.html.
        matrix_gz=config["tables_dir"].joinpath(
            "deeptools-promoters", "{macs2_mode}", "scale-regions", "{normalization_method}.gz"
        ),
        # optional output files
        matrix_tab=config["tables_dir"].joinpath(
            "deeptools-promoters", "{macs2_mode}", "scale-regions", "{normalization_method}.tab"
        ),
        matrix_bed=config["tables_dir"].joinpath(
            "deeptools-promoters", "{macs2_mode}", "scale-regions", "{normalization_method}.bed"
        ),
    log:
        config["log_dir"].joinpath(
            "deeptools-promoters", "compute_matrix", "{macs2_mode}", "scale-regions", "{normalization_method}.log"
        ),
    params:
        # required argument, choose "scale-regions" or "reference-point"
        command="scale-regions",
        # optional parameters
        extra=f"--regionBodyLength 500 -b 250 -a 250 --samplesLabel {' '.join(pep.sample_table['sample_name'].tolist())} --verbose",
    resources:
        mem_mb_per_cpu=lambda wildcards, attempt: 150 + 100 * (attempt - 1),
        runtime=lambda wildcards, attempt: 180 * attempt


rule plot_promoters:
    input:
        matrix_gz=config["tables_dir"].joinpath(
            "deeptools-promoters", "{macs2_mode}", "{computeMatrix_mode}", "{normalization_method}.gz"
        )
    output:
        # Please note that --outFileSortedRegions and --outFileNameMatrix are exclusively defined via output files.
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotheatmap.html.
        heatmap_img=report(
            config["pictures_dir"].joinpath("heatmap-promoters", "{macs2_mode}", "{computeMatrix_mode}", "{normalization_method}.pdf"),  # required
            category="Promoters heatmap"
        ),  # required
        # optional output files
        regions=config["pictures_dir"].joinpath("heatmap-promoters", "{macs2_mode}", "{computeMatrix_mode}", "{normalization_method}.bed"),
        heatmap_matrix=config["pictures_dir"].joinpath("heatmap-promoters", "{macs2_mode}", "{computeMatrix_mode}", "{normalization_method}.tab")
    log:
        config["log_dir"].joinpath("deeptools", "plot_promoters", "{macs2_mode}", "{computeMatrix_mode}", "{normalization_method}.log")
    params:
        # optional parameters
        f"--plotType=fill --colorMap Purples --samplesLabel {' '.join(pep.sample_table['sample_name'].tolist())}"
    retries: 5
    resources:
        runtime=30,
        mem_mb=lambda wildcards, input, attempt: 2000 * attempt + input.size_mb
    conda:
        "../env/deeptools.yml"
    shell:
        """
        if [ $(file -b {input.matrix_gz}) = "empty" ]
        then
            touch {output.heatmap_img}
            touch {output.heatmap_matrix}
            touch {output.regions}
            echo "No promoter peaks detected" |& tee {log} > $(dirname {output.heatmap_img})/README
        else
            plotHeatmap \
            -m {input.matrix_gz} \
            -o {output.heatmap_img} \
            --outFileSortedRegions {output.regions} \
            --outFileNameMatrix {output.heatmap_matrix} \
            {params} |& tee {log}
        fi 
        """