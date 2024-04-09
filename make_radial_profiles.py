#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 18:44:06 2024

@author: aellien
"""
from astropy.io import fits
from astropy.visualization import *
from photutils.profiles import *
from photutils.centroids import centroid_quadratic
from astropy.convolution import convolve
from photutils.segmentation import make_2dgaussian_kernel
from scipy.ndimage import gaussian_filter
import numpy as np
import matplotlib.pyplot as plt
import glob as glob
import os as os


if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/home/aellien/Euclid_ERO/data/Euclid-NISP-Stack-ERO-Abell2390.DR3'
    path_scripts = '/home/aellien/Euclid_ERO/Euclid_ERO_scripts'
    path_wavelets = '/home/aellien/Euclid_ERO/wavelets/out4/'
    path_plots = '/home/aellien/Euclid_ERO/plots'
    path_analysis = '/home/aellien/Euclid_ERO/analysis/'
        
    for input_file in glob.glob(os.path.join(path_data, '*crop.fits')):
        
        hdu = fits.open(input_file)
        oim = hdu[0].data
        xs, ys = oim.shape
        print(input_file)
        
        nf = input_file.split('/')[-1][:-5]
        nfp = os.path.join(path_wavelets, nf + 'synth.icl.bcgwavsizesepmask_006_010.fits')
        icl = fits.getdata(nfp)
        nfp = os.path.join(path_wavelets, nf + 'synth.icl_err.bcgwavsizesepmask_006_010.fits')
        icl_err = fits.getdata(nfp)
        
        kernel = make_2dgaussian_kernel(20.0, size=101)  # FWHM = 3.0
        #icl_smooth = convolve(icl, kernel)
        icl_smooth = gaussian_filter(icl, 50)
        
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (15,5))
        ax1.imshow(oim, origin = 'lower', cmap = 'grey', norm = ImageNormalize(oim, interval = ZScaleInterval(), stretch = LinearStretch()))
        ax2.imshow(icl_smooth, origin = 'lower', cmap = 'grey', norm = ImageNormalize(icl, interval = MinMaxInterval(), stretch = LogStretch()))
        
        plt.suptitle('Bande %s'%nf[12])
        
        
        #edge_radii = np.logspace(0, 3, 50, base = 10)
        edge_radii = np.linspace(0, 1000, 50)
        xycen = centroid_quadratic(icl)
        rp = RadialProfile(icl, xycen, edge_radii, error=icl_err, mask=None)
        mu = -2.5 * np.log10(rp.profile / (0.3**2) ) + 30
        mu_err = abs(-2.5 * np.log10( rp.profile_error / (0.3**2) ) * rp.profile_error / rp.profile / (0.3**2))

        mu[mu>30] = np.nan
        rad_kpc = rp.radius * 0.3 * 3.68
        
        mask = ~np.isnan(mu)
        mu = mu[mask]
        mu_err = mu_err[mask]
        rad_kpc = rad_kpc[mask]
        ax3.invert_yaxis()
        ax3.errorbar(rad_kpc, mu, yerr = mu_err, label = '%s band'%nf[12], \
                     linestyle='None', fmt = 'o', markersize = 4, capsize=2, capthick=1, alpha = 0.7, barsabove = True)
        ax3.fill_between(rad_kpc, mu - mu_err, mu + mu_err, color='gray', alpha=0.2)

        
        ax3.set_xlabel('r[kpc]')
        ax3.set_ylabel(r'$\mu$ [mag/arcsec$^2$]')
        ax3.set_title('BCG+ICL radial profile')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(path_plots,'icl_rec_rp_%s_band.png'%nf[12]), format = 'png')        

    plt.show()
