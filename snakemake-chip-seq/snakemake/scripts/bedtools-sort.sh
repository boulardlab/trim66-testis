#!/usr/bin/env bash

set -x 

bedtools sort \
-faidx "${snakemake_input[genome]}" \
-i "${snakemake_input[bdg]}" \
>"${snakemake_output}" \
2>"${snakemake_log}"