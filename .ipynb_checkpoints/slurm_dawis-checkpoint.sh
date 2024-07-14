#!/usr/bin/env bash
sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=euclid_ERO
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output /n03data/ellien/Euclid_ERO/logs/%x.%j.out 
#SBATCH --error  /n03data/ellien/Euclid_ERO/logs/%x.%j.err

source /home/ellien/.bashrc
conda activate dawis

python -W"ignore" /home/ellien/Euclid_ERO/Euclid_ERO_scripts/dawis_ERO.py $1
exit 0
EOT
