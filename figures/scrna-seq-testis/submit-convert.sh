#!/usr/bin/env bash

set -e

JOBS=""
for i in $(find -type f -name 'convert-*.R' | grep -v -i 'emtab6946'); do
  name=$(echo $i | sed -r 's/\.\/convert-(.*)\.R/\1/g')
  name=${name^^}
  JOBID=$(sbatch \
    -t 30:00 \
    -c 1 \
    --mem-per-cpu=64G \
    -A boulard \
    -J $name \
    -e log/${i%.R}-%a.log \
    -o log/${i%.R}-%a.log \
    -a 0,1 \
    --wrap="echo Task \$SLURM_ARRAY_TASK_ID; module load R-bundle-Bioconductor; Rscript $i \$SLURM_ARRAY_TASK_ID" | \
    sed 's/Submitted batch job //' )
  JOBS="$JOBS:$JOBID"
done

echo "Dependency list:$(echo $JOBS | tr ':' ' ')"
sbatch \
    -d afterok$JOBS \
    -t 2:00:00 \
    -c 1 \
    --mem=32G \
    -A boulard \
    -J integrate \
    -e log/integrate.log \
    -o log/integrate.log \
    --wrap="sleep 5; module purge; module load R-bundle-Bioconductor; Rscript integrate.R"
