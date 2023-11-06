#!/usr/bin/env bash

set -x

bedtools \
map \
-a "${snakemake_input[tiles]}" \
-b "${snakemake_input[sample]}" \
-c 4 \
-o "${snakemake_params[fun]}" \
-g "${snakemake_input[genome]}" \
-sorted \
>"${snakemake_output}" \
2>"${snakemake_log}"