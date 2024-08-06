#!/bin/bash
# Crop image (hand made with box region)

path_data=/home/aellien/Euclid_ERO/data/Euclid-NISP-Stack-ERO-Abell2390.DR3/Abell2390_NISP_starsubtraction
path_scripts=/home/aellien/Euclid_ERO/Euclid_ERO_scripts

for nfp in /home/aellien/Euclid_ERO/data/Euclid-NISP-Stack-ERO-Abell2390.DR3/Abell2390_NISP_starsubtraction/final/Euclid-NISP-H-ERO-Abell2390-LSB.DR3.mocks.crop.synth.cirrus_lvl5.bkg_sub.fits
do
    cd ${path_data}/vignets

    echo $nfp
    vignet_num=0
    fn="$(basename -- $nfp)"
    while read -r ra dec
    do 
        
        vignet_num=$(($vignet_num+1))
        out=${fn:0:-5}_vignet_${vignet_num}.fits
        echo $out
        astcrop --mode=wcs --center=$ra,$dec \
                --width=0.2,0.2 \
                -h2 \
                --output=$out \
                --primaryimghdu \
                $nfp
    done < /home/aellien/Euclid_ERO/analysis/callum_mock_clusters_circles_pos_wcs.txt
done

for nfp in /home/aellien/Euclid_ERO/data/Euclid-NISP-Stack-ERO-Abell2390.DR3/Abell2390_NISP_starsubtraction/Euclid-NISP-H-ERO-Abell2390-LSB.DR3.mock.revrot.fits
do
    cd ${path_data}/vignets

    echo $nfp
    vignet_num=0
    fn="$(basename -- $nfp)"
    while read -r ra dec
    do 
        
        vignet_num=$(($vignet_num+1))
        out=${fn:0:-5}_vignet_${vignet_num}.fits
        echo $out
        astcrop --mode=wcs --center=$ra,$dec \
                --width=0.2,0.2 \
                -h1 \
                --output=$out \
                --primaryimghdu \
                $nfp
    done < /home/aellien/Euclid_ERO/analysis/callum_mock_clusters_circles_pos_wcs.txt
done

cd $path_scripts


cd $path_scripts
