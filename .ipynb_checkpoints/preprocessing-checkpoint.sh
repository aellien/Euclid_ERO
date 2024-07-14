#!/bin/bash
# conda activate gnuastro

for f in /home/aellien/Euclid_ERO/data/Euclid-NISP-Stack-ERO-Abell2390.DR3/Euclid-NISP-?-ERO-Abell2390-LSB.DR3.crop.fits
do

warp=${f:0:-5}_warp.fits
astwarp $f -h0 --scale=1./2.0 --output=${warp}
astfits $warp --copy=1 --output=${warp:0:-5}_input.fits --primaryimghdu

done