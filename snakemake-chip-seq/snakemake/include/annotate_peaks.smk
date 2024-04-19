
rule mkgenome_for_homer:
    input:
        config["genome_fasta_filter"],
    output:
        temp(str(config["genome_fasta_filter"]).replace(".gz", "")),
    threads: 1
    resources:
        runtime=10,
        mem_mb=2000,
    conda:
        "../env/bash.yml"
    shell:
        """
        pigz -d -k {input}
        """


use rule mkgenome_for_homer as mkgtf_for_homer with:
    input:
        config["genome_gtf"],
    output:
        temp(str(config["genome_gtf"]).replace(".gz", "")),


def get_genome_for_homer(wildcards):
    p = str(config["genome_fasta_filter"])
    if p.endswith(".gz"):
        return p.replace(".gz", "")
    else:
        return p


def get_gtf_for_homer(wildcards):
    p = str(config["genome_gtf"])
    if p.endswith(".gz"):
        return p.replace(".gz", "")
    else:
        return p


rule homer_annotatepeaks:
    input:
        peaks=config["peaks_dir"].joinpath(
            "{macs2_mode}", "{sample}", "{sample}_peaks.bed"
        ),
        genome=get_genome_for_homer,
        # optional input files
        gtf=get_gtf_for_homer,  # implicitly sets the -gtf flag
        # gene="", # implicitly sets the -gene flag for gene data file to add gene expression or other data types
        #motif_files="peaks_refs/motives.txt", # implicitly sets the -m flag
        # filter_motiv="", # implicitly sets the -fm flag
        # center="",  # implicitly sets the -center flag
        #nearest_peak="peaks_refs/b.peaks", # implicitly sets the -p flag
        # tag="",  # implicitly sets the -d flag for tagDirectories
        # vcf="", # implicitly sets the -vcf flag
        # bed_graph="", # implicitly sets the -bedGraph flag
        # wig="", # implicitly sets the -wig flag
        # map="", # implicitly sets the -map flag
        # cmp_genome="", # implicitly sets the -cmpGenome flag
        # cmp_Liftover="", # implicitly sets the -cmpLiftover flag
        # advanced_annotation=""  # optional, implicitly sets the -ann flag, see http://homer.ucsd.edu/homer/ngs/advancedAnnotation.html
    output:
        annotations=config["peaks_dir"].joinpath(
            "{macs2_mode}", "{sample}", "{sample}_homer.txt"
        ),
        # optional output, implicitly sets the -matrix flag, requires motif_files as input
        # matrix=multiext("{sample}",
        #                 ".count.matrix.txt",
        #                 ".ratio.matrix.txt",
        #                 ".logPvalue.matrix.txt",
        #                 ".stats.txt"
        #                 ),
        # optional output, implicitly sets the -mfasta flag, requires motif_files as input
        # mfasta="{sample}_motif.fasta",
        # # optional output, implicitly sets the -mbed flag, requires motif_files as input
        # mbed="{sample}_motif.bed",
        # # optional output, implicitly sets the -mlogic flag, requires motif_files as input
        # mlogic="{sample}_motif.logic"
    threads: 2
    params:
        mode="tss " + config["genome"]["label"],  # add tss, tts or rna mode and options here, i.e. "tss mm8"
        extra=lambda wildcards, output: f"-gid -annStats {str(output).replace('_homer.txt', '_homer.annStats.txt')}",
        # optional params, see http://homer.ucsd.edu/homer/ngs/annotation.html
    resources:
        runtime=20,
        mem_mb=8000,
    log:
        config["log_dir"].joinpath("annotatePeaks", "{macs2_mode}", "{sample}.log"),
    wrapper:
        "v1.21.1/bio/homer/annotatePeaks"
