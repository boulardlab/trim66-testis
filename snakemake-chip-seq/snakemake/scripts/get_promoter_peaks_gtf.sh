#!/usr/bin/env bash

set -x
set -e

GENES=$(cat ${snakemake_input[annotations]} | \
grep -v PeakID | \
grep -i promoter | \
cut -f 12 | \
sort -u)

if [ ${#GENES} -gt 0 ]; then
    PAT=$(for GENE in $GENES; do echo -n "-e $GENE "; done)
    grep $PAT ${snakemake_input[gtf]} > ${snakemake_output[0]}
else
    touch ${snakemake_output[0]}
fi