rule bowtie2_build:
    input:
        ref=str(config["genome_fasta"]).replace(".fa", ".filter.fa"),
    output:
        protected(
            multiext(
                str(config["bowtie2_index_root"]),
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            )
        ),
    log:
        config["log_dir"].joinpath("bowtie2_build/build.log"),
    params:
        extra="",  # optional parameters
    retries: 5
    threads: 8
    resources:
        mem_mb_per_cpu=lambda wildcards, input, threads, attempt: 1500 * attempt,
        runtime=lambda wildcards, input, threads, attempt: 120 * attempt,
    wrapper:
        "v1.20.0/bio/bowtie2/build"


def build_bowtie_sample(wildcards):
    protocol = get_protocol(wildcards, wildcards.sample)
    if protocol == "single-end":
        return [config["trim_dir"].joinpath("{sample}.fastq.gz")]
    elif protocol == "paired-end":
        return [
            config["trim_dir"].joinpath("{sample}_1.fastq.gz"),
            config["trim_dir"].joinpath("{sample}_2.fastq.gz"),
        ]


rule bowtie2:
    # group:
    #     "alignment"
    input:
        idx=multiext(
            str(config["bowtie2_index_root"]),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
        sample=build_bowtie_sample,
    output:
        temp(config["alignments_dir"].joinpath("raw", "{sample}.bam")),
    log:
        config["log_dir"].joinpath("bowtie2/{sample}.log"),
    params:
        extra="--local --very-sensitive --no-mixed --no-discordant --dovetail",  # optional parameters
    retries: 5
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: (500 * attempt)
        + (input.size // 1000000),
        runtime=lambda wildcards, attempt: 720 * attempt,
    wrapper:
        "v1.20.0/bio/bowtie2/align"


def get_samtools_sort_mem(wildcards, input, attempt):
    mem = (20000 * attempt) + input.size_mb
    return mem


# def get_samtools_sort_mem_formatted(wildcards, input, resources):
#     num = resources["mem_mb"]
#     for unit in ["", "K", "M", "G"]:
#         if abs(num) < 1000.0:
#             return f"-m {num:d}{unit}"
#         num /= 1000.0
#     raise ValueError(
#         "Check Samtools sort memory requirement. "
#         + "Are your really requesting more than 1000G of memory?"
#     )


rule samtools_sort:
    input:
        config["alignments_dir"].joinpath("raw", "{sample}.bam"),
    output:
        temp(config["alignments_dir"].joinpath("sorted", "{sample}.bam")),
    log:
        config["log_dir"].joinpath("samtools_sort/{sample}.log"),
    threads: 8
    retries: 5
    params:
        # extra=lambda wildcards, input, resources: get_samtools_sort_mem_formatted(
        #     wildcards, input, resources
        # ),
        extra="",
    resources:
        mem_mb=lambda wildcards, input, attempt: get_samtools_sort_mem(
            wildcards, input, attempt
        ),
        runtime=lambda wildcards, attempt: 30 * attempt,
    wrapper:
        "v1.20.0/bio/samtools/sort"


rule samtools_index_bowtie2:
    # group:
    #     "alignment"
    input:
        config["alignments_dir"].joinpath("sorted", "{sample}.bam"),
    output:
        temp(config["alignments_dir"].joinpath("sorted", "{sample}.bam.bai")),
    log:
        config["log_dir"].joinpath("samtools_index/{sample}.log"),
    params:
        extra="",  # optional params string
    threads: 8  # This value - 1 will be sent to -@
    resources:
        mem_mb=lambda wildcards, input: 300 + (input.size // 1000000),
        runtime=10,
    wrapper:
        "v1.20.0/bio/samtools/index"
