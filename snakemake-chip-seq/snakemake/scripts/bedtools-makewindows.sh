#!/usr/bin/env bash

set -x 

bedtools \
makewindows \
-g "${snakemake_input}" \
-w "${snakemake_params[window_size]}" \
>"${snakemake_output}" \
2>"${snakemake_log}"