#!/usr/bin/env bash
sbatch <<EOT
#!/bin/bash
#SBATCH --output /home/ellien/Euclid_ERO/Euclid_ERO_scripts/test_ray_slurm/worker_%x.%j.out 
#SBATCH --error  /home/ellien/Euclid_ERO/Euclid_ERO_scripts/test_ray_slurm/worker_%x.%j.err
#SBATCH --time=00:10:00
### This script works for any number of nodes, Ray will find and manage all resources
#SBATCH --nodes=1
### Give all resources on each node to a single Ray task, ray can manage the resources internally
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
# Load modules or your own conda environment here
conda init bash
source /home/ellien/.bashrc
conda activate dawis

echp $1
ray start --address $1
sleep 5

python -u dawis_ray_slurm_test.py

exit
EOT
