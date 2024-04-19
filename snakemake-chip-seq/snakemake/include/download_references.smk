rule download_fasta:
    output:
        temp(config["genome_fasta"]),
    params:
        url=config["genome"]["fasta_link"],
    log:
        config["log_dir"].joinpath("download_fasta.log"),
    conda:
        "../env/wget.yml"
    threads: 1
    resources:
        mem_mb_per_cpu=lambda wildcards, input, threads, attempt: 500 * attempt,
        runtime=60,
    retries: 5
    shell:
        """
        wget -O {output} -o {log} {params.url}
        """


use rule download_fasta as dowload_gtf with:
    output:
        protected(config["genome_gtf"]),
    params:
        url=config["genome"]["gtf_link"],
    log:
        config["log_dir"].joinpath("dowload_gtf.log"),


rule filter_index_genome:
    input:
        config["genome_fasta"],
    output:
        filter_fasta=protected(config["genome_fasta_filter"]),
        filter_fasta_index=protected(config["genome_fasta_filter_index"]),
    params:
        chromosomes=build_chromosome_names(),
    log:
        config["log_dir"].joinpath("filter_genome.log"),
    conda:
        "../env/bowtie2-samtools.yml"
    threads: 8
    resources:
        mem_mb_per_cpu=200,
        runtime=30,
        tmpdir="/tmp",
    shell:
        """
        if [[ "{input}" == *.gz ]]; then
            export TMPDIR="/tmp"
            I=$(mktemp -t genome.XXX.fa)
            pigz -p {threads} -d -c "{input}" > $I            
        else 
            I={input}
        fi

        samtools faidx $I {params.chromosomes} >"{output.filter_fasta}" 2>{log}
        samtools faidx "{output.filter_fasta}" 2>>{log}
        """
