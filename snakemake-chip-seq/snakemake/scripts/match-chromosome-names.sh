#!/usr/bin/env bash

set -x

CHR1=$(cut -f1 "${snakemake_input[genome]}" | sort | uniq 2>"${snakemake_log}")
CHR2=$(cut -f1 "${snakemake_input[b]}" | sort | uniq 2>>"${snakemake_log}")

SHARED=$(comm -12 <(echo "$CHR1") <(echo "$CHR2") | tr "\n" " " 2>>"${snakemake_log}")

awk \
-v x="$SHARED" \
'BEGIN{FS=OFS="\t";split(x,myarr," ")}{for(k in myarr){if($1==myarr[k]){print; break}}}' \
"${snakemake_input[b]}" \
>"${snakemake_output}" \
2>>"${snakemake_log}"
