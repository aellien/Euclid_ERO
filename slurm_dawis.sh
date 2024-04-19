#!/usr/bin/env bash
sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=euclid
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --output /n03data/ellien/Euclid_ICL/logs/%x.%j.out 
#SBATCH --error  /n03data/ellien/Euclid_ICL/logs/%x.%j.err

source /home/ellien/.bashrc
conda activate dawis

python -W"ignore" /home/ellien/Euclid_ICL/Euclid_scripts/dawis_ERO.py $1 $2 $3

exit 0
EOT