rule samtools_idxstats:
    input:
        bam=config["alignments_dir"].joinpath("filtered", "sorted", "{sample}.bam"),
        idx=config["alignments_dir"].joinpath("filtered", "sorted", "{sample}.bam.bai"),
    output:
        config["alignments_dir"].joinpath("filtered", "sorted", "{sample}.idxstats.txt"),
    log:
        config["log_dir"].joinpath("samtools_idxstats/{sample}.log"),
    params:
        extra="",  # optional params string
    threads: 1
    retries: 5
    resources:
        mem_mb=lambda wildcards, input, threads, attempt: (200 * attempt)
        + (input.size // 1000000),
        runtime=10,
    wrapper:
        "v1.20.0/bio/samtools/idxstats"
