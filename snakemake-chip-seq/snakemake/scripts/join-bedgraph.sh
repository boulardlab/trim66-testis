#!/usr/bin/env bash

join -t $'\t' -j1 -o1.2,1.3,1.4,1.5,2.5 \
<(cat ${snakemake_input[a]} | awk -v OFS="\t" '{print $1"-"$2"-"$3,$0}' | sort -k1,1) \
<(cat ${snakemake_input[b]} | awk -v OFS="\t" '{print $1"-"$2"-"$3,$0}' | sort -k1,1) \
>${snakemake_output} 2>{snakemake_log}