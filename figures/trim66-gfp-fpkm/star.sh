#!/usr/bin/env bash

#SBATCH -A boulard 
#SBATCH -p bigmem
#SBATCH -C rome
#SBATCH -n 16
#SBATCH -N 1 
#SBATCH --mem 128312
#SBATCH -t 12:00:00
#SBATCH --overcommit 
#SBATCH --nice 
#SBATCH --oversubscribe 
#SBATCH -J star
#SBATCH -e log/%x-%j.log
#SBATCH -o log/%x-%j.log

set -e

M1=$1
M2=${M1/_1_sequence/_2_sequence}
DIR=$(dirname $M1)
DIR=$(basename $DIR)
FN=$(basename $M1)
FN=${FN/_1_sequence.fq.gz/.}
PREFIX=alignments/$DIR/$FN
GENOME="references/trim66-gfp/STAR"


module load STAR/2.7.11a-GCC-12.3.0
STAR \
--runMode alignReads \
--outSAMtype BAM SortedByCoordinate \
--outTmpDir $TEMPDIR/STAR \
--runThreadN $SLURM_NTASKS \
--outBAMsortingThreadN $SLURM_NTASKS \
--genomeDir $GENOME \
--readFilesCommand zcat \
--outFileNamePrefix $PREFIX \
--limitBAMsortRAM $(( $SLURM_MEM_PER_NODE * 1000000 )) \
--genomeLoad NoSharedMemory \
--readFilesIn $M1 $M2 \
--bamRemoveDuplicatesType UniqueIdentical \
--outBAMcompression -1

module load SAMtools/1.18-GCC-12.3.0

samtools index \
-@ $SLURM_NTASKS \
"${PREFIX}Aligned.sortedByCoord.out.bam"

samtools coverage \
--region chr7:109484669-109485543 \
"${PREFIX}Aligned.sortedByCoord.out.bam" \
> "${PREFIX}coverage.txt"
