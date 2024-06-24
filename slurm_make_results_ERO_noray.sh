#!/bin/bash
#SBATCH --job-name=synth_Euclid_ERO
#SBATCH --nodelist=n09
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --output /n03data/ellien/Euclid_ERO/logs/%x.%j.out 
#SBATCH --error  /n03data/ellien/Euclid_ERO/logs/%x.%j.err

source /home/ellien/.bashrc
conda activate dawis

path_data='/n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3'
cd $path_data
echo $(pwd)
for nf in Euclid-NISP-?-ERO-Abell2390-LSB.DR3.crop.fits
do
echo "launch $nf."
python /home/ellien/Euclid_ERO/Euclid_ERO_scripts/make_results_ERO.py $nf 5 80 1000 6 10 9 10 &
done

exit 0