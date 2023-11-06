def get_blacklist_url(config_label):
    label = "mm10" if config_label == "mm9" else config_label
    return encode_blacklists[label]


use rule download_fasta as get_blacklist with:
    output:
        protected(config["encode_blacklist"]),
    params:
        url=get_blacklist_url(config["genome"]["label"]),
    log:
        config["log_dir"].joinpath("download_blacklist.log"),


rule complement_bed:
    input:
        config["encode_blacklist"],
        config["genome_fasta_filter_index"],
    output:
        str(config["encode_blacklist"]).replace(".bed.gz", ".complement.bed.gz"),
    log:
        config["log_dir"].joinpath("complement_bed.log"),
    conda:
        "../env/bedtools.yml"
    resources:
        runtime=30,
        mem_mb=2000,
    shell:
        """
        gzip -c -d {input[0]} | \
        bedtools sort -i - -g {input[1]} | \
        bedtools complement -i - -g {input[1]} | \
        gzip -c > {output}
        """


rule samtools_view_blacklist:
    input:
        bam=config["alignments_dir"].joinpath("rmdup", "{sample}.bam"),
        bai=config["alignments_dir"].joinpath("rmdup", "{sample}.bam.bai"),
        blacklist=str(config["encode_blacklist"]).replace(
            ".bed.gz", ".complement.bed.gz"
        ),
    output:
        bam=temp(
            config["alignments_dir"].joinpath(
                "filtered", "samtools-view", "{sample}.bam"
            )
        ),
        # idx=config["filter_dir"].joinpath("{sample}.bam.bai")
    log:
        config["log_dir"].joinpath(
            "filter_alignements", "samtools", "view", "{sample}.log"
        ),
    params:
        extra=f"-F 0x004 -q 1 -L {str(config['encode_blacklist']).replace('.bed.gz', '.complement.bed.gz')}",  # optional params string
        region="",  # optional region string
    threads: 7
    resources:
        runtime=20,
        mem_mb=lambda wildcards, input: 1000 + input.size_mb,
    wrapper:
        "v1.21.2/bio/samtools/view"


def get_bamtools_script(wildcards):
    sample = wildcards.sample
    protocol = get_protocol(wildcards, sample)
    if protocol == "single-end":
        protocol = "se"
    elif protocol == "paired-end":
        protocol = "pe"
    return config["assets_dir"].joinpath(Path(f"bamtools_filter_{protocol}.json"))


rule bamtools_filter:
    input:
        config["alignments_dir"].joinpath("filtered", "samtools-view", "{sample}.bam"),
    output:
        temp(config["alignments_dir"].joinpath("filtered", "bamtools", "{sample}.bam")),
    params:
        json=lambda wildcards: get_bamtools_script(wildcards),
        region="",  # optional parameter for defining a specific region, e.g. "chr1:500..chr3:750"
    log:
        config["log_dir"].joinpath(
            "filter_alignments", "bamtools", "filter", "{sample}.log"
        ),
    resources:
        runtime=60,
        mem_mb=2000,
    wrapper:
        "v1.21.2/bio/bamtools/filter_json"


use rule samtools_sort as samtools_sort_filtered with:
    input:
        config["alignments_dir"].joinpath("filtered", "bamtools", "{sample}.bam"),
    output:
        config["alignments_dir"].joinpath("filtered", "sorted", "{sample}.bam"),


use rule samtools_index_rmdup as samtools_index_filtered with:
    input:
        config["alignments_dir"].joinpath("filtered", "sorted", "{sample}.bam"),
    output:
        config["alignments_dir"].joinpath("filtered", "sorted", "{sample}.bam.bai"),
