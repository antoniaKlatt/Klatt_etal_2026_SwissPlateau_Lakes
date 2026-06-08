#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --qos=6hours
#SBATCH --mem=20G
#SBATCH --error=cutadapt.log
#SBATCH --output=cutadapt.log
#SBATCH --job-name=cutadapt
#SBATCH --cpus-per-task=24
#SBATCH --mail-user=jan.waelchli@unibas.ch
#SBATCH --mail-type=ALL

## --------------------------------------------------------------------
## primer cutting
## --------------------------------------------------------------------

#load modules
module load cutadapt/4.2-GCCcore-11.3.0

#primers
primer_forward=GGACAGAAAGACCCTATGAA
primer_reverse=TCAGCCTGTTATCCCTAGAG

#create log folder
mkdir -p noprimer noprimer_log

#loop over each file
for fwd_file in ../00_data/*R1*; do

  #get run, name and reverse file
  run=$(basename ${fwd_file} | cut -f2 -d "_")
  name=$(basename ${fwd_file} | cut -f4 -d "_")
  rev_file=$(echo ${fwd_file} | sed 's/R1/R2/')

  #trim
  cutadapt -g ^${primer_forward} \
            -G ^${primer_reverse} \
            -o noprimer/${run}_${name}_F.fastq.gz \
            -p noprimer/${run}_${name}_R.fastq.gz \
            -e 0.1 \
            --cores 24 \
            --no-indels \
            ${fwd_file} ${rev_file} > noprimer_log/${name}.log 2>&1 #log file

done