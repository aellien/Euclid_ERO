#!/bin/bash
#SBATCH --job-name=euclid_ERO
#SBATCH --nodelist=n07
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output /n03data/ellien/Euclid_ERO/logs/%x.%j.out 
#SBATCH --error  /n03data/ellien/Euclid_ERO/logs/%x.%j.err

source /home/ellien/.bashrc
conda activate dawis
ray start --head --port=6379 &
sleep 10
python -u -W"ignore" /home/ellien/Euclid_ERO/Euclid_ERO_scripts/dawis_ERO_ray_test.py
exit 0