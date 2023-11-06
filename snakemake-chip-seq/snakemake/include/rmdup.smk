
rule mark_duplicates:
    # group:
    #     "rmdup"
    input:
        bams=config["alignments_dir"].joinpath("sorted", "{sample}.bam"),
        bai=config["alignments_dir"].joinpath("sorted", "{sample}.bam.bai"),
    output:
        bam=config["alignments_dir"].joinpath("rmdup", "{sample}.bam"),
        metrics=config["alignments_dir"].joinpath("rmdup", "{sample}.metrics.txt"),
    log:
        config["log_dir"].joinpath("picard/dedup/{sample}.log"),
    params:
        extra=config["picard-markduplicates"]["params"],
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    retries: 5
    resources:
        mem_mb=lambda wildcards, input: 1024 + (input.size // 1000000),
        runtime=lambda wildcards, attempt: 30 * attempt,
    threads: 1
    wrapper:
        "v1.20.0/bio/picard/markduplicates"


rule samtools_index_rmdup:
    # group:
    #     "rmdup"
    input:
        config["alignments_dir"].joinpath("rmdup", "{sample}.bam"),
    output:
        config["alignments_dir"].joinpath("rmdup", "{sample}.bam.bai"),
    log:
        config["log_dir"].joinpath("samtools_index/{sample}.log"),
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    resources:
        mem_mb=lambda wildcards, input: 300 + (input.size // 1000000),
        runtime=10,
    wrapper:
        "v1.20.0/bio/samtools/index"
