#!/bin/bash
#SBATCH --job-name=euclid_icl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --output /n03data/ellien/Euclid_ICL/logs/synthesis.%x.%j.out 
#SBATCH --error  /n03data/ellien/Euclid_ICL/logs/synthesis.%x.%j.err

source /home/ellien/.bashrc
conda activate dawis

python /home/ellien/Euclid_ICL/Euclid_scripts/synthesis_simulations.py

exit 0