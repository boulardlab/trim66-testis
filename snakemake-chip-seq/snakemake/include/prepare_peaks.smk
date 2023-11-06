import random
import os
import glob


def get_summits(wildcards):
    samples = get_samples(wildcards, as_list=False)
    return expand(
        config["peaks_dir"].joinpath(
            f"{wildcards.macs2_mode}", "{sample}", "refined", "{sample}_peaks.bed"
        ),
        sample=samples,
    )


rule prepare_peaks:
    wildcard_constraints:
        macs2_mode="model|nomodel_nobroad",
    input:
        summits=get_summits,
        genome=config["genome_fasta_filter_index"],
    output:
        merged_peaks=config["peaks_dir"].joinpath("{macs2_mode}", "roi.bed"),
    params:
        extension=100,  # extend both directions
    conda:
        "../env/bedtools.yml"
    log:
        config["log_dir"].joinpath("prepare_peaks/{macs2_mode}.log"),
    threads: 1
    resources:
        mem_mb_per_cpu=500,
        runtime=30,
    shell:
        """
        cat {input.summits} | \
        grep -v '#' | \
        sed '/^$/d' | \
        awk -F "\\t" -v OFS="\\t" '$1!="chr"{{print $1,$2,$3}}' | \
        sort -k1,1 -k2,2n | \
        bedtools slop -i - -g {input.genome} -b {params.extension} | \
        bedtools merge -i - | \
        bedtools sort -g {input.genome} -i - > {output}
        """


rule slop_summits:
    input:
        summits=config["peaks_dir"].joinpath(
            "{macs2_mode}", "{sample}", "refined", "{sample}_peaks.bed"
        ),
        genome=config["genome_fasta_filter_index"],
    output:
        config["peaks_dir"].joinpath(
            "{macs2_mode}", "{sample}", "slopped", "{sample}_peaks.bed"
        ),
    params:
        extension=200,
    conda:
        "../env/bedtools.yml"
    log:
        config["log_dir"].joinpath("slop_peaks/{macs2_mode}/{sample}.log"),
    threads: 1
    resources:
        mem_mb_per_cpu=500,
        runtime=20,
    shell:
        """
        bedtools sort -g {input.genome} -i {input.summits} | \
        bedtools slop -i - -g {input.genome} -b {params.extension} > {output}
        """


checkpoint intervene_venn:
    input:
        get_sloped_peaks_bed,
    output:
        plot=config["peaks_dir"].joinpath(
            "{macs2_mode}", "intervene", "Intervene_upset.pdf"
        ),
        folder=directory(config["peaks_dir"].joinpath("{macs2_mode}", "intervene")),
        sets=directory(
            config["peaks_dir"].joinpath("{macs2_mode}", "intervene", "sets")
        ),
    params:
        sample_names=",".join(pep.sample_table["sample_name"].tolist()),
    log:
        config["log_dir"].joinpath("intervene/{macs2_mode}.log"),
    conda:
        "../env/intervene.yml"
    threads: 1
    resources:
        mem_mb_per_cpu=500,
        runtime=10,
    shell:
        "intervene upset -i {input} --names={params.sample_names} -o {output.folder} --figtype pdf --save-overlaps"


rule screenshot_config:
    input:
        bigwig=get_bw_by_norm,
        bed=get_sloped_peaks_bed,
        calls=config["peaks_dir"].joinpath(
            "{macs2_mode}", "intervene", "sets", "{set}.bed"
        ),
        genome=config["genome_fasta_filter_index"],
    output:
        batch=config["peaks_dir"].joinpath(
            "{macs2_mode}",
            "intervene",
            "sets",
            "{normalization_method}",
            "{set}.igv.batch",
        ),
    params:
        genome_label=config["genome"]["label"],
        snapshot_directory=lambda wildcards: config["pictures_dir"].joinpath(
            "peaks", wildcards.macs2_mode, wildcards.normalization_method
        ),
        sample_table=pep.sample_table.to_dict(),
    conda:
        "../env/deeptools.yml"
    threads: 1
    resources:
        runtime=20,
        mem_mb=2000,
    log:
        config["log_dir"].joinpath(
            "intervene",
            "screenshot_config",
            "{macs2_mode}",
            "{normalization_method}",
            "{set}.log",
        ),
    script:
        "../scripts/generate_igv_script.py"


rule bedtools_multicov:
    input:
        unpack(get_filtered_bams),
        bed=config["peaks_dir"].joinpath(
            "{macs2_mode}", "intervene", "sets", "{set}.bed"
        ),
    output:
        config["peaks_dir"].joinpath(
            "{macs2_mode}", "intervene", "sets", "{set}.multicov.txt"
        ),
    log:
        config["log_dir"].joinpath("intervene", "multicov", "{macs2_mode}", "{set}.log"),
    conda:
        "../env/bedtools.yml"
    threads: 1
    resources:
        runtime=10,
        mem_mb=2000,
    shell:
        """
        cat {input.bed} | \
        awk -v OFS="\t" '{{print $1,$2,$3,$4}}' | \
        bedtools multicov \
        -bams {input.bam} \
        -bed - \
        > {output} |& \
        tee {log} 
        """


use rule screenshot_config as screenshot_config_bam with:
    input:
        unpack(get_filtered_bams),
        bed=get_peaks_bed,
        calls=config["peaks_dir"].joinpath(
            "{macs2_mode}", "intervene", "sets", "{set}.bed"
        ),
        coverage=config["peaks_dir"].joinpath(
            "{macs2_mode}", "intervene", "sets", "{set}.multicov.txt"
        ),
        genome=config["genome_fasta_filter_index"],
    output:
        batch=config["peaks_dir"].joinpath(
            "{macs2_mode}",
            "intervene",
            "sets",
            "bam",
            "{set}.igv.batch",
        ),
    params:
        genome_label=config["genome"]["label"],
        snapshot_directory=lambda wildcards: config["pictures_dir"].joinpath(
            "peaks", wildcards.macs2_mode, "intervene", "bam", wildcards.set
        ),
        sample_table=pep.sample_table.to_dict(),
    log:
        config["log_dir"].joinpath(
            "intervene", "screenshot_config", "{macs2_mode}", "bam", "{set}.log"
        ),


rule igv_screenshot:
    wildcard_constraints:
        normalization_method="|".join(normalization_methods + ["bam"]),
    input:
        batchfile=config["peaks_dir"].joinpath(
            "{macs2_mode}",
            "intervene",
            "sets",
            "{normalization_method}",
            "{set}.igv.batch",
        ),
    output:
        temp(
            touch(
                config["pictures_dir"].joinpath(
                    "peaks",
                    "{macs2_mode}",
                    "intervene",
                    "{normalization_method}",
                    "{set}",
                    "igv.done",
                )
            )
        ),
    conda:
        "../env/igv.yml"
    log:
        config["log_dir"].joinpath("igv/{macs2_mode}/{normalization_method}/{set}.log"),
    params:
        session_port=60151 + random.randint(0, 1000),
    retries: 3
    threads: 2
    resources:
        runtime=lambda wildcards, attempt: 20 + 10 * attempt,
        mem_mb_per_cpu=lambda wildcards, attempt: 1000 + 1500 * attempt
        if attempt > 1
        else 2000,
    envmodules:
        "Xvfb/21.1.6-GCCcore-12.2.0",
    shell:
        """
        xvfb-run \
        --auto-servernum \
        --server-args="-screen 0 3200x2400x24" \
        igv \
        -p {params.session_port} \
        -b {input.batchfile} |& \
        tee {log}
        """


def get_all_igv_flags(wildcards):
    intervene_venn_output = checkpoints.intervene_venn.get(**wildcards).output[2]
    return expand(
        config["pictures_dir"].joinpath(
            "peaks",
            wildcards.macs2_mode,
            "intervene",
            wildcards.normalization_method,
            "{intervene_set}",
            "igv.done",
        ),
        intervene_set=glob_wildcards(
            os.path.join(intervene_venn_output, "{intervene_set}.bed")
        ).intervene_set,
    )


rule flag_igv:
    wildcard_constraints:
        normalization_method="|".join(normalization_methods + ["bam"]),
    input:
        get_all_igv_flags,
    output:
        touch(
            config["pictures_dir"].joinpath(
                "peaks",
                "{macs2_mode}",
                "intervene",
                "{normalization_method}",
                "igv.done",
            )
        ),
