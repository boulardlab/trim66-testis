## This will work with both se and pe samples
rule fastqc:
    input:
        config["reads_dir"].joinpath("{sample}.fastq.gz"),
    output:
        html=report(
            config["results_dir"].joinpath("qc/fastqc/{sample}.html"),
            category="FastQC",
            subcategory="{sample}",
        ),
        zip=config["results_dir"].joinpath("qc/fastqc/{sample}_fastqc.zip"),
    params:
        "--quiet " + config["fastqc"]["params"],
    log:
        config["log_dir"].joinpath("fastqc/{sample}.log"),
    retries: 5
    threads: 8
    resources:
        mem_mb_per_cpu=lambda wildcards, input, threads, attempt: 200 * attempt,
        runtime=lambda wildcards, input, threads, attempt: 20 * attempt,
    wrapper:
        "v1.20.0/bio/fastqc"


def get_fastqc(wildcards):
    s = get_samples(wildcards)
    return expand(
        config["results_dir"].joinpath("qc", "fastqc", "{sample}_fastqc.zip"), sample=s
    )


rule multiqc_fastqc_raw:
    input:
        get_fastqc,
    output:
        report(
            config["results_dir"].joinpath("qc", "multiqc", "fastqc-raw.html"),
            category="MultiQC",
            subcategory="FastQC raw reads",
        ),
    params:
        extra="",  # Optional: extra parameters for multiqc.
        use_input_files_only=True,
    log:
        config["log_dir"].joinpath("multiqc/fastqc-raw.log"),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10,
    wrapper:
        "v1.20.0/bio/multiqc"
