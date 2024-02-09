#!/bin/bash
#
#

for file in Euclid-NISP-H-ERO-Abell2390-LSB.DR3.crop.fits Euclid-NISP-J-ERO-Abell2390-LSB.DR3.crop.fits Euclid-NISP-Y-ERO-Abell2390-LSB.DR3.crop.fits
      echo "Launch Dawis on file $file"
      qsub qsub_dawis_ERO.sh -v n=${file},ncl=A2390 #${file:0:-5}
done
