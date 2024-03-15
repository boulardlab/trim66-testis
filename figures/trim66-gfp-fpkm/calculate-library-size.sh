#!/usr/bin/env bash

#SBATCH -A boulard 
#SBATCH -p htc-el8
#SBATCH -C bergamo
#SBATCH -n 2
#SBATCH -N 1 
#SBATCH --mem 2992
#SBATCH -t 30:00
#SBATCH --overcommit 
#SBATCH --nice 
#SBATCH --oversubscribe 
#SBATCH -J picard-CollectAlignmentSummaryMetrics
#SBATCH -e log/%x-%j.log
#SBATCH -o log/%x-%j.log

set -e

BAM=$1
GENOME="references/trim66-gfp/mm10-gfp.fa"

module load picard/3.1.0-Java-17
java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
      R=$GENOME \
      I=$BAM \
      O=${BAM/bam/alignment_metrics.txt}

grep "PF_BASES" ${BAM/bam/alignment_metrics.txt} > ${BAM/bam/library_size.txt}

