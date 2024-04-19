#!/usr/bin/env bash

rm log/*

find alignments \
-type f \
-name '*.Aligned.sortedByCoord.out.bam' | \
xargs -I {} sbatch calculate-library-size.sh {}


