#!/bin/bash

#SBATCH --time=3-00:00:00
#SBATCH --qos=1week
#SBATCH --mem=40g
#SBATCH --output=log.out
#SBATCH --error=log.error
#SBATCH --job-name=blastn
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=jan.waelchli@unibas.ch
#SBATCH --mail-type=ALL

#load modules

#get the taxa with blast
module load Boost.MPI/1.83.0-gompi-2023b
module load BLAST+/2.14.1-gompi-2023b

#blast
blastn -db /scicore/data/managed/BLAST_FASTA/latest/nt -query ../../02_dadapipe/10_OUTPUT/bacteria_SEQ100.fasta \
-evalue 1e-5 \
-max_target_seqs 5 \
-num_threads $SLURM_CPUS_PER_TASK \
-outfmt "6 qseqid pident staxids" \
-out taxid.tab
