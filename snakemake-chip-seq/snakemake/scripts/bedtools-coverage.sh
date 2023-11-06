#!/usr/bin/env bash

# bedtools coverage \
# -a "${snakemake_input[roi]}" \
# -b "${snakemake_input[bam]}" \
# -g "${snakemake_input[genome]}" \
# -sorted \
# -counts > "${snakemake_output}" 

bedtools multicov -q {snakemake_params[mapq]} \
-bams ${snakemake_input[bam]} \
-bed ${snakemake_input[roi]} > ${snakemake_output}
