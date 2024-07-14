#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 15:09:09 2024

@author: aellien
"""
from astropy.io import fits
from astropy.visualization import *
from photutils.profiles import *
from photutils.centroids import centroid_quadratic
from astropy.convolution import convolve
from photutils.segmentation import make_2dgaussian_kernel
from astropy.visualization import make_lupton_rgb
from scipy.ndimage import gaussian_filter
import numpy as np
import matplotlib.pyplot as plt
import glob as glob
import os as os
import matplotlib as mpl
from rebin import rebin
from vorbin.voronoi_2d_binning import *

def plot_RGB_fullfield(r, g, b, artist_list):
    plt.figure()
    sdc = np.copy([r, g, b])
    sdc[sdc < 0.4] = 0
    sdc[sdc > 0] = 1
    print(sdc.shape)
    master_sup = np.prod(sdc, axis = 0)
    b *= master_sup*255.0/b.max()
    g *= master_sup*255.0/g.max()
    r *= master_sup*255.0/r.max()
    image = make_lupton_rgb(r, g, b, Q = 10, stretch = 0.5)
    for a in artist_list:
        #a.set_marker('x')
        a.set_markeredgecolor('blue')
        a.set_markersize('2')
        plt.gca().add_artist(a)
        
    plt.imshow(image, origin= 'lower')
    plt.savefig('/home/aellien/Euclid_ERO/plots/RGB_raw.png', format = 'png', dpi = 1000)

def display_pixels(x, y, xmin, xmax, ymin, ymax, counts, pixelSize):
    """
    Display pixels at coordinates (x, y) coloured with "counts".
    This routine is fast but not fully general as it assumes the spaxels
    are on a regular grid. This needs not be the case for Voronoi binning.

    """
    #xmin, xmax = np.min(x), np.max(x)
    #ymin, ymax = np.min(y), np.max(y)
    nx = round((xmax - xmin)/pixelSize) + 1
    ny = round((ymax - ymin)/pixelSize) + 1
    img = np.full((nx, ny), np.nan)  # use nan for missing data
    j = np.round((x - xmin)/pixelSize).astype(int)
    k = np.round((y - ymin)/pixelSize).astype(int)
    img[j, k] = counts

    plt.imshow(img, interpolation='none', cmap='prism',
               extent=[xmin - pixelSize/2, xmax + pixelSize/2,
                       ymin - pixelSize/2, ymax + pixelSize/2], origin = 'lower')
    
    return img

if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/home/aellien/Euclid_ERO/data/Euclid-NISP-Stack-ERO-Abell2390.DR3'
    path_scripts = '/home/aellien/Euclid_ERO/Euclid_ERO_scripts'
    path_wavelets = '/home/aellien/Euclid_ERO/wavelets/out7/'
    path_plots = '/home/aellien/Euclid_ERO/plots'
    path_analysis = '/home/aellien/Euclid_ERO/analysis/'
    
    #plt.ion()

    
    for input_file in glob.glob(os.path.join(path_data, '*-J-*crop.fits')):
        
        hdu = fits.open(input_file)
        oim = hdu[0].data
        xs, ys = oim.shape
        print(input_file)
        
        nf = input_file.split('/')[-1][:-4]
        nfp = os.path.join(path_wavelets, nf + 'synth.bcgwavsizesepmask_005_080.fits')

        hdu = fits.open(nfp)
        icl = hdu[1].data
        
        
        icl[icl < 0.4] = 0
        
        bin_icl = rebin(icl, 10, 10, 'sum')
        xs, ys = np.shape(bin_icl)
        snr = np.sum(bin_icl) / 1000
        
        signal_idx = np.where(bin_icl != 0.)
        flat_icl = bin_icl[signal_idx]
        x = signal_idx[0]
        y = signal_idx[1]
        
        noise = np.ones(flat_icl.shape)        
        bin_number, x_gen, y_gen, x_bar, y_bar, sn, nPixels, scale = voronoi_2d_binning(x, y, flat_icl, noise, snr, cvt=False, pixelsize=1, plot=False, quiet=False, sn_func=None, wvt=True)
        print(bin_number.shape, nPixels.shape)
        plt.show()
        
        #display_pixels(x, y, 0, xs, 0, ys, flat_icl, 1)
        #plt.show()
    
    vorl = []
    colors = []
    rgb = []
    #2 = H, 1 = Y, 0 = J
    #lambda Y->J->H 
    for input_file in glob.glob(os.path.join(path_data, '*crop.fits')):
        
        hdu = fits.open(input_file)
        oim = hdu[0].data
        print(input_file)
        
        nf = input_file.split('/')[-1][:-4]
        nfp = os.path.join(path_wavelets, nf + 'synth.bcgwavsizesepmask_005_100.fits')

        hdu = fits.open(nfp)
        icl = hdu[1].data
        gal = hdu[2].data
        
        bin_icl = rebin(icl, 10, 10, 'sum')
        flat_icl = bin_icl[signal_idx]
        
        vor_icl = np.copy(flat_icl)
        for idx in np.unique(bin_number):
             i = np.where(bin_number == idx)
             flux = np.sum(flat_icl[i])
             vor_icl[i] = flux
        vor_icl = display_pixels(x, y, 0, xs, 0, ys, vor_icl, 1)
        vorl.append(vor_icl)
        
        vor_icl[vor_icl==0] = 1e-10
        mag_vor_icl = -2.5 * np.log10(vor_icl) + 30
        #mag_vor_icl[mag_vor_icl > 30] = 0
        colors.append( mag_vor_icl)
        rgb.append(gal)
        
    
    # RGB
    vrgb = make_lupton_rgb(vorl[2], vorl[0], vorl[1], Q=10, stretch=0.1)
    plt.imshow(vrgb, origin='lower')
    plt.savefig(os.path.join(path_plots, 'RGB_vor_YJH.png'))
    
    r, g, b = rgb[2], rgb[0], rgb[1]
    plot_RGB_fullfield(r, g, b, [])
    
    # HY
    plt.figure(1)
    cmap = mpl.colormaps.get_cmap('seismic_r')
    cmap.set_bad(color='white')
    plt.imshow(colors[2] - colors[1], cmap = cmap, origin='lower', vmax = 2, vmin = -2)
    plt.suptitle("H-Y color map")
    plt.colorbar()
    plt.savefig(os.path.join(path_plots, 'color_vor_HY.png'))
    
    # HJ
    plt.figure(2)
    plt.imshow(colors[2] - colors[0], cmap = cmap, origin='lower', vmax = 2, vmin = -2)
    plt.suptitle("H-J color map")
    plt.colorbar()
    plt.savefig(os.path.join(path_plots, 'color_vor_HJ.png'))
    
    # JY
    plt.figure(3)
    plt.imshow(colors[0] - colors[1], cmap = cmap, origin='lower', vmax = 2, vmin = -2)
    plt.suptitle("J-Y color map")
    plt.colorbar()
    plt.savefig(os.path.join(path_plots, 'color_vor_JY.png'))
    
    
    plt.show()
