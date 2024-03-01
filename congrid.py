#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 15:30:09 2024

@author: aellien
"""

import dawis as d
import os
import glob
from astropy.io import fits

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/home/aellien/Euclid_ERO/data/Euclid-NISP-Stack-ERO-Abell2390.DR3'
    path_scripts = '/home/aellien/Euclid_ERO/Euclid_ERO_scripts'
    path_wavelets = '/home/aellien/Euclid_ERO/wavelets/out5/'
    path_plots = '/home/aellien/Euclid_ERO/plots'
    path_analysis = '/home/aellien/Euclid_ERO/analysis/'
    
    anx, any = 5000, 5000
    for nfp in glob.glob(os.path.join(path_data, '*crop.fits')):
        
        hdu = fits.open(nfp)
        im = hdu[0].data
        cim = d.congrid.congrid(im, (anx, any), method = 'spline')
        
        hduo = fits.PrimaryHDU(cim, header = hdu[0].header)
        hduo.writeto(nfp[:-5]+'_congrid.fits', overwrite = True)