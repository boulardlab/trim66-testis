gsize = {"mm10": "mm", "mm9": "mm", "hg19": "hg", "hg38": "hg"}

macs_mode_params = {
    "model": "",
    "nomodel_broad": "--nomodel --broad --broad-cutoff 0.05",
    "nomodel_nobroad": "--nomodel",
}


def get_macs2_call_input_files(wildcards):
    sample = wildcards.sample
    protocol = get_protocol(wildcards, sample)
    ret = {
        "bam": config["alignments_dir"].joinpath("filtered", "sorted", "{sample}.bam")
    }
    if wildcards.macs2_mode != "model" and protocol == "single-end":
        ret["fragsize"] = config["alignments_dir"].joinpath(
            "filtered", "sorted", "{sample}.fragsize.txt"
        )
    return ret


def get_macs2_mode_params(wildcards, snakemake_inputs):
    sample = wildcards.sample
    protocol = get_protocol(wildcards, sample)
    ret = config["macs2-callpeak"]["params"]
    ret += f" {macs_mode_params[wildcards.macs2_mode]}"
    if (
        hasattr(snakemake_inputs, "fragsize")
        and wildcards.macs2_mode != "model"
        and protocol == "single-end"
    ):
        with open(snakemake_inputs.fragsize) as handler:
            fragsize = handler.readline().strip()
            ret += f" --extsize {fragsize}"
    return ret


rule macs2_call_nomodel:
    wildcard_constraints:
        macs2_mode="nomodel.*",
    input:
        unpack(get_macs2_call_input_files),
    output:
        config["peaks_dir"].joinpath("{macs2_mode}", "{sample}", "{sample}_peaks.xls"),
        config["peaks_dir"].joinpath("{macs2_mode}", "{sample}", "{sample}_summits.bed"),
        # config["peaks_dir"].joinpath("{macs2_mode}", "{sample}", "{sample}_model.r"),
    params:
        gsize=gsize[config["genome"]["label"]],
        outdir=lambda wildcards: config["peaks_dir"].joinpath(wildcards.macs2_mode),
        others=lambda wildcards, input: get_macs2_mode_params(wildcards, input),
        protocol=lambda wildcards: get_protocol(wildcards, wildcards.sample),
    conda:
        "../env/macs2.yml"
    log:
        config["log_dir"].joinpath("macs2_callpeak/{macs2_mode}/{sample}.log"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, threads, attempt: (200 * attempt)
        + (input.size // 1000000),
        runtime=30,
    shell:
        """
        macs2 callpeak \
        --treatment {input.bam} \
        --name {wildcards.sample} \
        --format AUTO \
        --gsize {params.gsize} \
        --outdir {params.outdir}/{wildcards.sample} \
        {params.others} |& tee {log}
        
        cd {params.outdir}/{wildcards.sample}
        if [ -f {wildcards.sample}_model.r ]; then 
          Rscript {wildcards.sample}_model.r |& tee -a {log}
        fi
        """


use rule macs2_call_nomodel as macs2_call_model with:
    wildcard_constraints:
        macs2_mode="(?<!no)model.*",


rule convert_macs2_to_bed:
    input:
        peaks=config["peaks_dir"].joinpath(
            "{macs2_mode}", "{sample}", "{sample}_peaks.xls"
        ),
        genome=config["genome_fasta_filter_index"],
    output:
        config["peaks_dir"].joinpath("{macs2_mode}", "{sample}", "{sample}_peaks.bed"),
    conda:
        "../env/bedtools.yml"
    log:
        config["log_dir"].joinpath(
            "convert_macs2_to_bed/{macs2_mode}", "{sample}", "{sample}_peaks.log"
        ),
    threads: 1
    resources:
        runtime=5,
        mem_mb=2000,
    shell:
        """
        cat {input.peaks} |\
        grep -v '#' |\
        grep -v '^$' |\
        tail -n+2 |\
        awk -v OFS="\\t" '{{print $1,$2,$3,$10,".","*"}}' | \
        bedtools sort -g {input.genome} -i - >{output} 
        """


# rule macs2_plot_model:
#     input:
#         config["peaks_dir"].joinpath("{macs2_mode}", "{sample}", "{sample}_model.r"),
#     output:
#         report(
#             config["peaks_dir"].joinpath(
#                 "{macs2_mode}", "{sample}", "{sample}_model.pdf"
#             ),
#             category="MACS2 model",
#             subcategory="{sample}",
#         ),
#     conda:
#         "../env/tidyverse.yml"
#     log:
#         config["log_dir"].joinpath(
#             "macs2_callpeak/{macs2_mode}/plot_model/{sample}.log"
#         ),
#     threads: 1
#     resources:
#         mem_mb=lambda wildcards, input, threads, attempt: (500 * attempt)
#         + (input.size // 1000000),
#         runtime=10,
#     shell:
#         """
#         BN=$(dirname {input})
#         cd $BN
#         Rscript {input}
#         """


rule macs2_refine:
    input:
        peaks=config["peaks_dir"].joinpath(
            "{macs2_mode}", "{sample}", "{sample}_peaks.bed"
        ),
        bam=config["alignments_dir"].joinpath("filtered", "sorted", "{sample}.bam"),
        bai=config["alignments_dir"].joinpath("filtered", "sorted", "{sample}.bam.bai"),
    output:
        config["peaks_dir"].joinpath(
            "{macs2_mode}", "{sample}", "refined", "{sample}_peaks.bed"
        ),
    conda:
        "../env/macs2.yml"
    log:
        config["log_dir"].joinpath("macs2_callpeak/{macs2_mode}/{sample}.log"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, threads, attempt: (200 * attempt)
        + (input.size // 1000000),
        runtime=30,
    shell:
        """
        macs2 refinepeak \
        -b {input.peaks} \
        -i {input.bam} \
        -f AUTO \
        -o {output} 2>&1 >{log}
        """


rule count_peaks:
    input:
        config["peaks_dir"].joinpath("{macs2_mode}", "{sample}", "{sample}_peaks.xls"),
    output:
        config["peaks_dir"].joinpath("{macs2_mode}", "{sample}", "npeaks.txt"),
    conda:
        "../env/bash.yml"
    log:
        config["log_dir"].joinpath("count_peaks/{macs2_mode}/{sample}.log"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, threads, attempt: (200 * attempt)
        + (input.size // 1000000),
        runtime=10,
    shell:
        """
        grep -v '#' {input} | \
        grep -v '^$' | \
        tail -n +2 | \
        wc -l |& \
        tee {log} > {output}
        """


rule plot_peak_counts:
    input:
        get_npeaks_files,
    output:
        report(
            config["pictures_dir"].joinpath("peaks", "{macs2_mode}", "npeaks.pdf"),
            category="Number of peaks",
        ),
    params:
        sample_sheet=pep.sample_table.to_dict(orient="list"),
        contrast_column="genotype",
        contrast_levels=["KO", "WT"],
    conda:
        "../env/tidyverse.yml"
    log:
        config["log_dir"].joinpath("plot_peak_counts/{macs2_mode}.log"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, threads, attempt: (2000 * attempt)
        + (input.size // 1000000),
        runtime=10,
    script:
        "../scripts/plot_npeaks.R"
