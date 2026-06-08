#!/bin/bash

#SBATCH --time=1-00:00:00
#SBATCH --qos=1day
#SBATCH --mem=1g
#SBATCH --error=dadapipe.log
#SBATCH --job-name=dadapipe
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=jan.waelchli@unibas.ch
#SBATCH --mail-type=ALL

#create folder for log files
mkdir -p log/jobs 2> /dev/null

# #decide if only bacteria, only fungi or bacteria and fungi samples
# #which taxa is present
# taxa=$(cat 01_INPUT/02_design/design.csv | awk -F ',' '{print $2}' | sort | uniq | grep -v "taxa")
# bacteria=$(echo $taxa | grep -c "bacteria")
# fungi=$(echo $taxa | grep -c "fungi")
# #create config label
# if [[ $bacteria -gt 0 ]]; then label='"bacteria"'; fi
# if [[ $fungi -gt 0 ]]; then label='"fungi"'; fi
# if [[ $bacteria -gt 0 &&  $fungi -gt 0 ]]; then label='["bacteria", "fungi"]'; fi
# #add or overwrite label
# cat dada_config.yaml | grep -v "^taxa:" > temp.yaml #taxa: line removed from dada_config if there
# echo taxa: $label >> temp.yaml #concat new label to end
# cat temp.yaml > dada_config.yaml #overwrite
# rm temp.yaml


#activate mambaforge and snakemake
source /scicore/home/schlae0003/GROUP/software/miniforge3/bin/activate
conda activate Jan_snakemake_8.14.0

#unlock folder if needed
snakemake --unlock

#start snakemake
wd=$(pwd)
snakemake --workflow-profile ${wd}


#move log files
#mv dadapipe.out log
mv dadapipe.log log
mv SLURM.log log
