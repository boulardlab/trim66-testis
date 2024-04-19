#!/usr/bin/env bash

#SBATCH -A boulard 
#SBATCH -C bergamo 
#SBATCH -n 16 
#SBATCH -N 1 
#SBATCH --mem 23936 
#SBATCH -p htc-el8 
#SBATCH -t 1-00:00:00
#SBATCH --overcommit 
#SBATCH --nice 
#SBATCH --oversubscribe 

module load STAR/2.7.11a-GCC-12.3.0
STAR --runMode genomeGenerate \
--outTmpDir $TEMPDIR/STAR \
--runThreadN 16 \
--genomeDir references/trim66-gfp/STAR \
--genomeFastaFiles references/trim66-gfp/mm10-gfp.fa \
--genomeSAindexNbases 8 \
--genomeLoad NoSharedMemory
