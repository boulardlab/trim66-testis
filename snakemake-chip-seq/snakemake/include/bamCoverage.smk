def get_bamCoverage_input(wildcards):
    samples = get_samples(wildcards, as_list=False)
    df = pep.sample_table[samples == wildcards.sample]
    protocol = df["sequencing_protocol"].iloc[0]
    if protocol == "single-end":
        ret = {
            "bam": config["alignments_dir"].joinpath(
                "filtered", "sorted", "{sample}.bam"
            ),
            "bai": config["alignments_dir"].joinpath(
                "filtered", "sorted", "{sample}.bam.bai"
            ),
            "fragsize": config["alignments_dir"].joinpath(
                "filtered", "sorted", "{sample}.fragsize.txt"
            ),
        }
    else:
        ret = {
            "bam": config["alignments_dir"].joinpath(
                "filtered", "sorted", "{sample}.bam"
            ),
            "bai": config["alignments_dir"].joinpath(
                "filtered", "sorted", "{sample}.bam.bai"
            ),
        }
    return ret


def get_bamCoverage_params(wildcards, input):
    ret = config["bamCoverage"]["params"]
    if (
        wildcards.normalization_method.lower() == "raw"
        or wildcards.normalization_method.lower() == "none"
    ):
        ret += " --normalizeUsing None"
    else:
        ret += f" --normalizeUsing {wildcards.normalization_method.upper()}"

    if hasattr(input, "fragsize"):
        with open(input.fragsize) as handler:
            fragsize = handler.readline().strip()
            ret += f" --extendReads {fragsize}"

    return ret


rule bamCoverage:
    wildcard_constraints:
        normalization_method="|".join(normalization_methods),
    input:
        unpack(get_bamCoverage_input),
    output:
        config["bigwig_dir"].joinpath("{normalization_method}", "{sample}.bw"),
    log:
        config["log_dir"].joinpath(
            "bamCoverage", "{normalization_method}", "{sample}.log"
        ),
    conda:
        "../env/deeptools.yml"
    params:
        effective_genome_size=get_effective_genome_size,
        minMappingQuality=config["bamCoverage"]["minMappingQuality"],
        others=lambda wildcards, input: get_bamCoverage_params(wildcards, input),
    threads: 16
    resources:
        mem_mb_per_cpu=lambda wildcards, attempt: 150 + 100 * (attempt - 1),
        runtime=lambda wildcards, attempt: 60 + 30 * (attempt - 1),
    shell:
        """        
        bamCoverage \
        -b {input.bam} \
        -o {output} \
        -of bigwig \
        -p {threads} \
        --effectiveGenomeSize {params.effective_genome_size} \
        --minMappingQuality {params.minMappingQuality} \
        {params.others} |& \
        tee {log}
        """
