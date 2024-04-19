ruleorder: trimmomatic_pe > trimmomatic_se


rule trimmomatic_se:
    input:
        config["reads_dir"].joinpath("{sample}.fastq.gz"),
    output:
        config["trim_dir"].joinpath("{sample}.fastq.gz"),
        summary=config["trim_dir"].joinpath("{sample}.summary.txt"),
    log:
        config["log_dir"].joinpath("trimmomatic/{sample}.log"),
    params:
        # list of trimmers (see manual)
        trimmer=config["trimmomatic"]["trimmer"]["se"].split(" "),
        # optional parameters
        extra=f"-trimlog {config['trim_dir'].joinpath('{sample}.summary.txt')}",
        # optional compression levels from -0 to -9 and -11
        compression_level="-9",
    params:
        # list of trimmers (see manual)
        trimmer=config["trimmomatic"]["trimmer"]["se"].split(" "),
        # optional parameters
        extra=config["trimmomatic"]["params"]["se"]
        + f"-trimlog {config['trim_dir'].joinpath('{sample}.summary.txt')}",
        compression_level="-9",
    threads: 32
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    retries: 5
    resources:
        mem_mb=lambda widlcards, attempt: 1024 * attempt,
        runtime=lambda widlcards, attempt: 120 * attempt,
    wrapper:
        "v1.20.0/bio/trimmomatic/se"


rule trimmomatic_pe:
    input:
        r1=config["reads_dir"].joinpath("{sample}_1.fastq.gz"),
        r2=config["reads_dir"].joinpath("{sample}_2.fastq.gz"),
    output:
        r1=config["trim_dir"].joinpath("{sample}_1.fastq.gz"),
        r2=config["trim_dir"].joinpath("{sample}_2.fastq.gz"),
        # reads where trimming entirely removed the mate
        r1_unpaired=config["trim_dir"].joinpath("{sample}_1.unpaired.fastq.gz"),
        r2_unpaired=config["trim_dir"].joinpath("{sample}_2.unpaired.fastq.gz"),
        summary=config["trim_dir"].joinpath("{sample}.summary.txt"),
    log:
        "logs/trimmomatic/{sample}.log",
    params:
        # list of trimmers (see manual)
        trimmer=config["trimmomatic"]["trimmer"]["pe"].split(" "),
        # optional parameters
        extra=config["trimmomatic"]["params"]["pe"]
        + f"-trimlog {config['trim_dir'].joinpath('{sample}.summary.txt')}",
        compression_level="-9",
    threads: 32
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    retries: 5
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt,
        runtime=lambda wildcards, attempt: 120 * attempt,
    wrapper:
        "v1.20.0/bio/trimmomatic/pe"


use rule fastqc as fastqc_trim with:
    input:
        config["trim_dir"].joinpath("{sample}.fastq.gz"),
    output:
        html=report(
            config["results_dir"].joinpath("qc/fastqc-trim/{sample}.html"),
            category="FastQC trimmed reads",
            subcategory="{sample}",
        ),
        zip=config["results_dir"].joinpath("qc/fastqc-trim/{sample}_fastqc.zip"),
    log:
        config["log_dir"].joinpath("fastqc-trim/{sample}.log"),


def get_fastqc_trimmed(wildcards):
    s = get_samples(wildcards)
    return expand(
        config["results_dir"].joinpath("qc", "fastqc-trim", "{sample}_fastqc.zip"),
        sample=s,
    )


def get_trimmomatic_report(wildcards):
    s = get_samples(wildcards)
    return expand(config["trim_dir"].joinpath("{sample}.summary.txt"), sample=s)


use rule multiqc_fastqc_raw as multiqc_fastqc_trim with:
    input:
        get_fastqc,
        get_trimmomatic_report,
    output:
        report(
            config["results_dir"].joinpath("qc", "multiqc", "fastqc-trimmed.html"),
            category="MultiQC",
            subcategory="FastQC trimmed reads",
        ),
    log:
        config["log_dir"].joinpath("multiqc/fastqc-trim.log"),
