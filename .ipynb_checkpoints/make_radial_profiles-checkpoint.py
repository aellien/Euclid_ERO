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
import pyregion as pyr
import os as os


if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/home/aellien/Euclid_ERO/data/Euclid-NISP-Stack-ERO-Abell2390.DR3'
    path_scripts = '/home/aellien/Euclid_ERO/Euclid_ERO_scripts'
    path_wavelets = '/home/aellien/Euclid_ERO/wavelets/out7/'
    path_plots = '/home/aellien/Euclid_ERO/plots'
    path_analysis = '/home/aellien/Euclid_ERO/analysis/'
    
    bcg_coo = [ 1763, 1759 ]
    outside_coo = [1360, 1016]
    rpl = []
    r = pyr.open('/home/aellien/Euclid_ERO/analysis/mask_top_half_cl.reg')
    
    
    for input_file in glob.glob(os.path.join(path_data, '*crop.fits')):
        
        hdu = fits.open(input_file)
        oim = hdu[0].data
        xs, ys = oim.shape
        print(input_file)
        
        nf = input_file.split('/')[-1][:-5]
        nfp = os.path.join(path_wavelets, nf + '.synth.bcgwavsizesepmask_005_100.fits')
        hdu = fits.open(nfp)
        
        icl = hdu[1].data
        icl_err = hdu[5].data
        gal = hdu[2].data
        gal_err = hdu[6].data
        
        m = r.get_mask(hdu = hdu[1])
        #icl[m] = 0
        #kernel = make_2dgaussian_kernel(20.0)  # FWHM = 3.0
        #icl = convolve(icl, kernel)
        #icl = gaussian_filter(icl, 20)
        
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (15,5))
        ax1.imshow(oim, origin = 'lower', cmap = 'grey', norm = ImageNormalize(oim, interval = ZScaleInterval(), stretch = LinearStretch()))
        ax2.imshow(icl, origin = 'lower', cmap = 'grey', norm = ImageNormalize(icl, interval = MinMaxInterval(), stretch = LogStretch()))
        
        plt.suptitle('Bande %s'%nf[12])
        
        
        #edge_radii = np.logspace(0, 3, 50, base = 10)
        edge_radii = np.linspace(0, 1000, 100)
        xycen = bcg_coo
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
        
        mag = -2.5 * np.log10(rp.profile)
        
        rpl.append((mag[mask],rad_kpc))
        
        

    # 0 - J, 1 - Y, 2 - H
    plt.figure()
    rad_kpc = rp.radius * 0.3 * 3.68
    idx = np.min([np.size(rpl[1][0]), np.size(rpl[0][0])])
    plt.plot(rad_kpc[:idx], rpl[0][0][:idx]-rpl[1][0][:idx], 'bo', label = 'Y-J')
    plt.title('color radial profile')

    idx = np.min([np.size(rpl[2][0]), np.size(rpl[1][0])])
    plt.plot(rad_kpc[:idx], rpl[2][0][:idx]-rpl[1][0][:idx], 'ro', label = 'H-J')

    idx = np.min([np.size(rpl[2][0]), np.size(rpl[0][0])])
    plt.plot(rad_kpc[:idx], rpl[2][0][:idx]-rpl[0][0][:idx], 'go', label = 'H-Y')
    
    plt.legend()
    plt.gca().invert_yaxis()
    plt.savefig(os.path.join(path_plots,'rcolor.png'), format = 'png')
    
    plt.figure()
    for input_file in glob.glob(os.path.join(path_data, '*crop.fits')):
        
        hdu = fits.open(input_file)
        oim = hdu[0].data
        xs, ys = oim.shape
        print(input_file)
        
        nf = input_file.split('/')[-1][:-5]
        nfp = os.path.join(path_wavelets, nf + '.synth.bcgwavsizesepmask_005_100.fits')
        hdu = fits.open(nfp)
        
        icl = hdu[1].data
        icl_err = hdu[5].data
        gal = hdu[2].data
        gal_err = hdu[6].data
        
        m = r.get_mask(hdu = hdu[1])
        #icl[m] = 0
                
        edge_radii = np.linspace(0, 1000, 100)
        xycen = bcg_coo
        rp = RadialProfile(icl, xycen, edge_radii, error=icl_err, mask=None)
        mu = -2.5 * np.log10(rp.profile / (0.3**2) ) + 30
        mu_err = abs(-2.5 * np.log10( rp.profile_error / (0.3**2) ) * rp.profile_error / rp.profile / (0.3**2))

        mu[mu>30] = np.nan
        rad_kpc = rp.radius * 0.3 * 3.68
        
        mask = ~np.isnan(mu)
        mu = mu[mask]
        mu_err = mu_err[mask]
        rad_kpc = rad_kpc[mask]
        plt.errorbar(rad_kpc, mu, yerr = mu_err, label = '%s band'%nf[12], \
                     linestyle='None', fmt = 'o', markersize = 4, capsize=2, capthick=1, alpha = 0.7, barsabove = True)
        plt.fill_between(rad_kpc, mu - mu_err, mu + mu_err, color='gray', alpha=0.2)

    ax = plt.gca()
    ax.invert_yaxis()

    ax.set_xlabel('r[kpc]')
    ax.set_ylabel(r'$\mu$ [mag/arcsec$^2$]')
    ax.set_title('BCG+ICL radial profile')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(path_plots,'icl_rec_rp_all_band.png'), format = 'png')
    plt.show()