rule fragment_size:
    input:
        bamFile=config["alignments_dir"].joinpath("filtered", "sorted", "{sample}.bam"),
        baiindex=config["alignments_dir"].joinpath(
            "filtered", "sorted", "{sample}.bam.bai"
        ),
    output:
        report=config["alignments_dir"].joinpath(
            "filtered", "sorted", "{sample}.fragsize.txt"
        ),
        figure=report(
            config["alignments_dir"].joinpath(
                "filtered", "sorted", "{sample}.fragsize.png"
            ),
            category="Fragment size estimation",
            subcategory="{sample}",
        ),
    params:
        max_delay=500,
    log:
        config["log_dir"].joinpath("estimate_fragment_size/{sample}.log"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, threads, attempt: (200 * attempt)
        + (input.size // 1000000),
        runtime=10,
    conda:
        "../env/csaw.yml"
    script:
        "../scripts/get_fragment_size.R"
