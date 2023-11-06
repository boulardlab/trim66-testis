def get_fingerprint_inputs(wildcards):
    d = get_filtered_bams(wildcards)
    return {"bam_files": d["bam"], "bam_idx": d["bam_idx"]}


rule plot_fingerprint:
    input:
        unpack(get_fingerprint_inputs),
    output:
        # Please note that --plotFile and --outRawCounts are exclusively defined via output files.
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotfingerprint.html.
        fingerprint=report(
            config["pictures_dir"].joinpath("fingerprint/fingerprint.png"),  # required
            category="Fingerprint"
        ),
        # optional output
        counts=config["pictures_dir"].joinpath("fingerprint/raw_counts.tab"),
        qc_metrics=config["pictures_dir"].joinpath("fingerprint/qc_metrics.txt"),
    log:
        config["log_dir"].joinpath("deeptools/plot_fingerprint.log"),
    params:
        # optional parameters
        f"--numberOfSamples 250000 --labels {' '.join(pep.sample_table['sample_name'].tolist())}",
    threads: 
        8
    resources:
        mem_mb=lambda wildcards, input, attempt: (200 * attempt) + (input.size // 1000000),
        runtime=30
    wrapper:
        "v1.20.0/bio/deeptools/plotfingerprint"
