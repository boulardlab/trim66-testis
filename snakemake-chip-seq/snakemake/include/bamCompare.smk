def build_sample_wildcard(wildcards, samples, wildcard_name, index):
    # if wildcards.serie in samples["pe"]:
    #     case_control = samples["pe"][wildcards.serie]["case-control"]
    # else:
    #     case_control = samples["se"][wildcards.serie]["case-control"]

    s = [x[index] for x in case_control]
    s = list(set(s))
    w = "|".join(s)
    w = "{" + wildcard_name + "," + w + "}"
    return w


def get_treatments(wildcards):
    return {
        "case": str(config["dedup_dir"].joinpath("{{sample}}" + ".bam")),
        "case_index": str(config["dedup_dir"].joinpath("{{sample}}" + ".bam.bai")),
    }


def get_controls(wildcards):
    return {
        "control": str(config["dedup_dir"].joinpath("{{sample}}" + ".bam")),
        "control_index": str(config["dedup_dir"].joinpath("{{sample}}" + ".bam.bai")),
    }


rule normalize:
    input:
        unpack(get_dedup_bams),
    output:
        config["bigwig_dir"].joinpath("{treatment}_vs_{control}.bw"),
    params:
        config["bamCompare"]["params"],
    threads: 2
    log:
        config["log_dir"].joinpath("bamCompare-{serie}-{treatment}_vs_{control}.log"),
    conda:
        "../env/deeptools.yml"
    shell:
        """
        bamCompare -b1 {input.case} -b2 {input.control} \
        -o {output} -of bigwig \
        -p {threads} \
        {params} |& \
        tee {log}
        """
