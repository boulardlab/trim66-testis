
rule samples_correlation:
    input:
        unpack(get_filtered_bams),
    output:
        config["results_dir"].joinpath("multiBamSummary.npz"),
    params:
        config["multiBamSummary"]["params"],
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: (200 * attempt)
        + (input.size // 1000000),
        runtime=30,
    log:
        config["log_dir"].joinpath("deeptools/multiBamSummary.log"),
    conda:
        "../env/deeptools.yml"
    shell:
        """
        multiBamSummary bins {params} \
        -p {threads} \
        --bamfiles {input.bam} \
        -o {output}
        """


rule plot_correlation_heatmap_with_summary:
    input:
        config["results_dir"].joinpath("multiBamSummary.npz"),
    output:
        report(
            config["pictures_dir"].joinpath("multiBamSummary", "heatmap.pdf"),
            category="MultiBam Summary",
            labels={"figure": "heatmap"},
        ),
        config["pictures_dir"].joinpath("multiBamSummary", "heatmap.txt"),
    params:
        config["plotCorrelation"]["params"],
        labels=" ".join(pep.sample_table["sample_name"].tolist()),
    log:
        config["log_dir"].joinpath("deeptools/plot_correlation.log"),
    conda:
        "../env/deeptools.yml"
    threads: 1
    retries: 5
    resources:
        mem_mb=lambda wildcards, input, attempt: (20000 * attempt)
        + (input.size // 1000000),
        runtime=lambda wildcards, attempt: 60 * attempt,
    shell:
        """
        plotCorrelation -in {input} \
        --whatToPlot heatmap \
        --plotFileFormat pdf \
        --outFileCorMatrix {output[1]} \
        --labels {params.labels} \
        {params[0]} \
        -o {output[0]}
        """


rule plot_pca:
    input:
        config["results_dir"].joinpath("multiBamSummary.npz"),
    output:
        report(
            config["pictures_dir"].joinpath("multiBamSummary", "pca.pdf"),
            category="MultiBam Summary",
            labels={"figure": "PCA"},
        ),
        config["pictures_dir"].joinpath("multiBamSummary", "pca.txt"),
    params:
        labels=" ".join(pep.sample_table["sample_name"].tolist()),
    log:
        config["log_dir"].joinpath("deeptools/plot_pca.log"),
    conda:
        "../env/deeptools.yml"
    threads: 1
    retries: 5
    resources:
        mem_mb=lambda wildcards, attempt: 20000 * attempt,
        runtime=lambda wildcards, attempt: 60 * attempt,
    shell:
        """
        plotPCA --corData {input} \
        --plotFile {output[0]} \
        --labels {params.labels} \
        --outFileNameData {output[1]} \
        """
