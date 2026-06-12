#!/bin/bash

#SBATCH --qos=1day
#SBATCH --time=08:00:00
#SBATCH --mem=80g
#SBATCH --output=run.out
#SBATCH --error=run.error
#SBATCH --job-name=tree_root_sand
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=jan.waelchli@unibas.ch
#SBATCH --mail-type=ALL

#align
module load MAFFT/7.520-GCC-12.3.0-with-extensions
mafft-linsi SEQ_phytoplankton.fasta > SEQ_phytoplankton_aligned.fasta

#create tree
source /scicore/home/schlae0003/GROUP/software/miniforge3/bin/activate
conda activate Jan_iqtree_3.0.1
iqtree -s SEQ_phytoplankton_aligned.fasta -T 16