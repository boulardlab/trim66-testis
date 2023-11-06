#!/usr/bin/env bash

set -x
set -e

cat ${snakemake_input[0]} | \
grep -v PeakID | \
grep -i promoter | \
awk -v OFS="\t" '{print $2,$3,$4,$1,$6,$5}' | \
sort -k2,2 -k3,3n > ${snakemake_output[0]}