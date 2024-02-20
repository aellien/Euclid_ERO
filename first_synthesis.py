#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 13:42:17 2024

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
import ray
import dawis as d
from astropy.io import fits
from astropy.visualization import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.morphology import binary_dilation
from scipy.stats import kurtosis

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def read_image_atoms( nfp, filter_it = None, verbose = False ):

    # Object lists
    if filter_it == None:
        opath = nfp + '*ol.it*.pkl'
        itpath = nfp + '*itl.it*.pkl'
    else:
        opath = nfp + '*ol.it' + filter_it  + '.pkl'
        itpath = nfp + '*itl.it' + filter_it + '.pkl'

    opathl = glob.glob(opath)
    opathl.sort()

    # Interscale tree lists

    itpathl = glob.glob(itpath)
    itpathl.sort()

    tol = []
    titl = []

    if verbose:
        print('Reading %s.'%(opath))
        print('Reading %s.'%(itpath))

    for i, ( op, itlp ) in enumerate( zip( opathl, itpathl )):

        if verbose :
            print('Iteration %d' %(i), end ='\r')

        ol = d.read_objects_from_pickle( op )
        itl = d.read_interscale_trees_from_pickle( itlp )

        for j, o in enumerate(ol):

            tol.append(o)
            titl.append(itl[j])

    return tol, titl

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def synthesis_fullfield( oim, nfp, xs, ys, n_levels, C, write_fits = True ):
    '''Synthesis of the full astronomical field (e.g. sum of all atoms)
    --- Args:
    oim         # Original astronomical field
    nfp         # root path of *.pkl
    gamma       # attenuation factor
    lvl_sep_big # wavelet scale at which gamma set to 1
    lvl_sep     # wavelet scale threshold for the separation
    xs, ys      # image size
    n_levels    # number of wavelet scales
    plot_vignet # plot pdf vignet of output
    --- Output:
    rec         # synthesis image with all atoms
    res         # residuals (original - rec)
    wei         # atom weight image
    '''
    # path, list & variables
    res = np.zeros( (xs, ys) )
    rec = np.zeros( (xs, ys) )
    wei = np.zeros( (xs, ys) )
    xc = xs / 2.
    yc = ys / 2.

    # Read atoms
    ol, itl = read_image_atoms( nfp, verbose = False )

    for j, o in enumerate(ol):

        lvlo = o.level
        x_min, y_min, x_max, y_max = o.bbox

        rec[ x_min : x_max, y_min : y_max ] += (o.image - C) #/!\

        # atom weight map
        o.image[o.image > 0.] = 1.
        wei[ x_min : x_max, y_min : y_max ] += ( o.image - C) #/!\

    res = oim - rec
    if write_fits == True:

        print('\nFULLFIELD -- write fits as %s*fits'%(nfp))

        hduo = fits.PrimaryHDU(rec)
        hduo.writeto( nfp + 'synth.restored.fits', overwrite = True )

        hduo = fits.PrimaryHDU(res)
        hduo.writeto( nfp + 'synth.residuals.fits', overwrite = True )

        hduo = fits.PrimaryHDU(wei)
        hduo.writeto( nfp + 'synth.weight.fits', overwrite = True )

    return rec, res, wei


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/home/aellien/Euclid_ERO/data/'
    path_scripts = '/home/aellien/Euclid_ERO/Euclid_ERO_scripts'
    path_wavelets = '/home/aellien/Euclid_ERO/wavelets/out4/'
    path_plots = '/home/aellien/Euclid_ERO/plots'
    path_analysis = '/home/aellien/Euclid_ERO/analysis/'
    
    n_levels = 11
    
    for input_file in glob.glob(os.path.join(path_data, 'Euclid-NISP-Stack-ERO-Abell2390.DR3/*crop*')):
        
        hdu = fits.open(input_file)
        oim = hdu[0].data
        xs, ys = oim.shape
        C = np.min(oim)
        oim -= C
        head = hdu[0].header
        
        nf = input_file.split('/')[-1][:-4]
        nfp = os.path.join(path_wavelets, nf)
        print(input_file)
        output = synthesis_fullfield( oim, nfp, xs, ys, n_levels, C = C, write_fits = True )
