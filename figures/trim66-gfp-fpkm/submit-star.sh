#!/usr/bin/env bash

mkdir log
rm log/*

find reads/raw/ \
-type d \
\( -name TRIM66_Dec2020 \
-or -name TRIM66_Apr2021 \
-or -name TRIM66_elongatedSpermatids_Ago2021_noeGFP2_noeGFP3_noeGFP6 \) | \
xargs -I {} find {} -name '*_1_sequence*' | \
xargs -I {} sbatch star.sh {}
