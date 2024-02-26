#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 13:44:37 2024

@author: aellien
"""
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Modules
import sys
sys.path.append("/home/aellien/dawis")
import glob as glob
import os
import numpy as np
import pyregion as pyr
import random
import pandas as pd
import dawis as d
import pandas as pd
import photutils as phut
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import *
from astropy.convolution import convolve

def try_iptd_seg():
    bkg_estimator = phut.background.MedianBackground()
    bkg = phut.background.Background2D(data, (50, 50), filter_size=(3, 3), bkg_estimator=bkg_estimator)
    data -= bkg.background
    kernel = phut.segmentation.make_2dgaussian_kernel(5.0, size=25)  # FWHM = 3.0
    convolved_data = convolve(data, kernel)
    
    '''norm = ImageNormalize(data, interval = ZScaleInterval())
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
    ax1.imshow(data, origin='lower', cmap='Greys_r', norm=norm)
    ax1.set_title('Background-subtracted Data')
    ax2.imshow(segment_map, origin='lower', cmap=segment_map.cmap, interpolation='nearest')
    ax2.set_title('Segmentation Image')'''

    threshold = 2. * bkg.background_rms
    finder = phut.segmentation.SourceFinder(npixels = 10, progress_bar = True )
    #segment_map = finder(convolved_data, threshold)
    segment_map = phut.segmentation.detect_sources(convolved_data, threshold, npixels=10)
    
    bboxl = [ x.extent for x in segment_map.bbox ]
    centerl = [ x.center for x in segment_map.bbox ]
    extentl = []
    for i, (bbox, center, area) in enumerate( zip(bboxl, centerl, segment_map.areas)):
        tot_area = (bbox[1] - bbox[0]) * (bbox[3] - bbox[2])
        extent = area / tot_area
        extentl.append(extent)
    extentl = np.array(extentl)
    
    sorted_labels = segment_map.labels[ np.argsort( segment_map.areas )][-20:]
    sorted_extentl = extentl[ np.argsort( segment_map.areas )][-20:]
    
    mask = np.zeros(data.shape)
    for i in sorted_labels:
        if segment_map.data[1800, 1800] != i:
            mask[np.where(segment_map.data == i)] = 1

    ''' norm = ImageNormalize(data, interval = AsymmetricPercentileInterval(40, 99.5), stretch = LogStretch())
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 12.5))
    ax1.imshow(data, origin='lower', cmap='Greys_r', norm=norm)
    ax1.set_title('Background-subtracted Data')
    ax2.imshow(mask, origin='lower', cmap=segment_map.cmap,
               interpolation='nearest')
    ax2.set_title('Segmentation Image')
    plt.savefig('/home/aellien/Euclid_ERO/plots/seg.png', format = 'png', dpi = 1000)
    '''
    iptd_im = np.copy(im)
    noise = d.sample_noise(iptd_im, n_sigmas = 3)
    draws = np.random.normal(noise.mean(), noise.std(), iptd_im.shape)
    mask *= draws
    iptd_im[ mask > 1 ] = 0.
    iptd_im += mask
    
    norm = ImageNormalize(iptd_im, interval = AsymmetricPercentileInterval(40, 99.5), stretch = LogStretch())
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 12.5))
    ax1.imshow(iptd_im, origin='lower', cmap='Greys_r', norm=norm)
    ax1.set_title('Background-subtracted Data')
    ax2.imshow(mask, origin='lower', cmap=segment_map.cmap,
               interpolation='nearest')
    ax2.set_title('Segmentation Image')
    plt.savefig('/home/aellien/Euclid_ERO/plots/seg.png', format = 'png', dpi = 1000)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/home/aellien/Euclid_ERO/data/Euclid-NISP-Stack-ERO-Abell2390.DR3'
    path_scripts = '/home/aellien/Euclid_ERO/Euclid_ERO_scripts'
    path_wavelets = '/home/aellien/Euclid_ERO/wavelets/out5/'
    path_plots = '/home/aellien/Euclid_ERO/plots'
    path_analysis = '/home/aellien/Euclid_ERO/analysis/'
    
    hdu = fits.open(os.path.join(path_data, 'Euclid-NISP-H-ERO-Abell2390-LSB.DR3.crop.fits'))
    im = hdu[0].data
    head = hdu[0].header
    w = WCS(head)
    
    data = np.copy(im)
    peaks_tbl = phut.detection.find_peaks(data, threshold = 10E4)
    peaks_tbl['peak_value'].info.format = '%.8g'
    
    size = 25
    hsize = (size - 1) / 2
    x = peaks_tbl['x_peak']  
    y = peaks_tbl['y_peak']  

    mask = ((x > hsize) & (x < (data.shape[1] -1 - hsize)) &\
        (y > hsize) & (y < (data.shape[0] -1 - hsize))) 