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
for nf in Euclid-NISP-?-ERO-Abell2390-LSB.DR3.crop.fits
do
echo "launch $nf."
    for lvl_sep in 5 6
    do
        for size_sep in 60 80 100 140
        do
            python /home/ellien/Euclid_ERO/Euclid_ERO_scripts/make_results_ERO_noray.py $nf $lvl_sep $size_sep 1000 6 10 9 10&
        done
    done
done
wait
exit 0