#!/bin/bash
for file in /n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3/vignets/*fits
do
  njobs=$((${njobs}+1))
  echo "Launch Dawis on file ${file}"
  n=$(basename "$file")
  ncl=${n:0:-5}
  indir=/n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3/vignets
  outdir=/n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3/wavelets/out14
  bash slurm_dawis_distributed.sh $indir $n $outdir
  
done
