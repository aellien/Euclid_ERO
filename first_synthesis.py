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
import dawis as d
from astropy.io import fits
from astropy.visualization import *
from scipy.stats import kurtosis

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def read_image_atoms( nfp, filter_it = None, verbose = True ):

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
def synthesis_fullfield( oim, nfp, xs, ys, n_levels, write_fits = True ):
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
    rec_dc = np.zeros((n_levels, xs, ys))
    xc = xs / 2.
    yc = ys / 2.

    # Read atoms
    #ol, itl = read_image_atoms( nfp, verbose = True ) MEMORY ISSUE
    
    ######################################## MEMORY
    opath = nfp + '*ol.it*.pkl'
    itpath = nfp + '*itl.it*.pkl'
    opathl = glob.glob(opath)
    opathl.sort()

    # Interscale tree lists

    itpathl = glob.glob(itpath)
    itpathl.sort()

    tol = []
    titl = []

    for i, ( op, itlp ) in enumerate( zip( opathl, itpathl )):
        
        ol = d.read_objects_from_pickle( op )
        itl = d.read_interscale_trees_from_pickle( itlp )
        ############################################################## MEMORY
        for j, o in enumerate(ol):
            
            lvlo = o.level
            x_min, y_min, x_max, y_max = o.bbox
            rec[ x_min : x_max, y_min : y_max ] += o.image
            rec_dc[ lvlo, x_min : x_max, y_min : y_max ] += o.image
            
            # atom weight map
            o.image[o.image > 0.] = 1.
            wei[ x_min : x_max, y_min : y_max ] += o.image
        
    res = oim - rec
    if write_fits == True:

        print('\nFULLFIELD -- write fits as %s*fits'%(nfp))

        hduo = fits.PrimaryHDU(rec)
        hduo.writeto( nfp + 'synth.restored.fits', overwrite = True )

        hduo = fits.PrimaryHDU(res)
        hduo.writeto( nfp + 'synth.residuals.fits', overwrite = True )

        hduo = fits.PrimaryHDU(wei)
        hduo.writeto( nfp + 'synth.weight.fits', overwrite = True )
        
        hduo = fits.PrimaryHDU(rec_dc)
        hduo.writeto( nfp + 'synth.restored_dc.fits', overwrite = True )

    return rec, res, wei

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def synthesis_bcgwavsizesep_with_masks( nfp, lvl_sep, lvl_sep_max, lvl_sep_bcg, size_sep_pix, xs, ys, rc_pix, n_levels, mscell, mscbcg, kurt_filt = True, write_fits = True ):
    '''Wavelet Separation + Spatial filtering.
    ICL --> Atoms with z > lvl_sep, with maximum coordinates within ellipse mask 'mscell' and with size > size_sep_pix.
    Galaxies --> Satellites + BCG, so a bit complicated:
        - Atoms not classified as ICL but with maximum coordinates within ellipse mask 'mscbcg'
        - Atoms within radius 'rc' of a member galaxy position
        - In case lvl_sep > 5 (arbitrary...), atoms with z > 5 and within 'mscell' are BCG
    Unclassified --> rest of atoms
    '''
    # path, list & variables
    icl = np.zeros( (xs, ys) )
    gal = np.zeros( (xs, ys) )
    im_art = np.zeros( (xs, ys) )
    im_unclass = np.zeros( (xs, ys) )

    icl_al = []
    gal_al = []
    noticl_al = []
    unclass_al = []

    #%
    at_test = []
    #%

    xc = xs / 2.
    yc = ys / 2.

    # Read atoms
    #ol, itl = read_image_atoms( nfp, verbose = True ) MEMORY ISSUE
    
    ######################################## MEMORY
    opath = nfp + '*ol.it*.pkl'
    itpath = nfp + '*itl.it*.pkl'
    opathl = glob.glob(opath)
    opathl.sort()

    # Interscale tree lists

    itpathl = glob.glob(itpath)
    itpathl.sort()

    tol = []
    titl = []

    for i, ( op, itlp ) in enumerate( zip( opathl, itpathl )):
        
        ol = d.read_objects_from_pickle( op )
        itl = d.read_interscale_trees_from_pickle( itlp )
        ############################################################## MEMORY

        # Kurtosis + ICL+BCG
        for j, o in enumerate(ol):
    
            x_min, y_min, x_max, y_max = o.bbox
            sx = x_max - x_min
            sy = y_max - y_min
            itm = itl[j].interscale_maximum
            xco = itm.x_max
            yco = itm.y_max
            lvlo = o.level
        
            if kurt_filt == True:
                k = kurtosis(o.image.flatten(), fisher=True)
                if k < 0:
                    im_art[ x_min : x_max, y_min : y_max ] += o.image
                    continue
    
            # Remove background
            if o.level >= lvl_sep_max:
                continue
    
            # ICL + BCG
            if mscell[xco, yco] == 1:
    
                # BCG
                xbcg, ybcg = xs, ys
                if mscbcg[xco, yco] == 1:
    
                    dr = np.sqrt( (xbcg - xco)**2 + (ybcg - yco)**2 )
                    if (o.level <= 3) & (dr < rc_pix):
    
                        icl[ x_min : x_max, y_min : y_max ] += o.image
                        icl_al.append([o, xco, yco])
    
                    elif o.level > 3 :
    
                        icl[ x_min : x_max, y_min : y_max ] += o.image
                        icl_al.append([o, xco, yco])
    
                # ICL
                else:
    
                    if (o.level >= lvl_sep) & (sx >= size_sep_pix) & (sy >= size_sep_pix):
                        
                        icl[ x_min : x_max, y_min : y_max ] += o.image
                        icl_al.append([o, xco, yco])
                        at_test.append([xco, yco])
                        
                    else:
                        
                        #gal[ x_min : x_max, y_min : y_max ] += o.image
                        noticl_al.append([o, xco, yco])
    
            else:
                noticl_al.append([ o, xco, yco ])
            
    if write_fits == True:
        print('\nWS + SF + SS -- ICL+BCG -- write fits as %s*'%(nfp))

        # write to fits
        hduo = fits.PrimaryHDU(icl)
        hduo.writeto( nfp + 'synth.icl.bcgwavsizesepmask_%03d_%03d.fits'%(lvl_sep, size_sep_pix), overwrite = True )

        hduo = fits.PrimaryHDU(gal)
        hduo.writeto( nfp + 'synth.gal.bcgwavsizesepmask_%03d_%03d.fits'%(lvl_sep, size_sep_pix), overwrite = True )

    return icl, gal

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3'
    path_scripts = '/home/ellien/Euclid_ERO/Euclid_ERO_scripts'
    path_wavelets = '/n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3/wavelets/out4/'
    path_plots = '/home/ellien/Euclid_ERO/plots'
    path_analysis = '/n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3/'
    
    n_levels = 11
    lvl_sep = 6
    lvl_sep_max = 12
    lvl_sep_bcg = 5
    size_sep_pix = 10
    rc_pix = 10
    mscell = fits.getdata(os.path.join(path_data,'mscell.fits'))
    mscbcg = fits.getdata(os.path.join(path_data,'mscbcg.fits'))

    for input_file in glob.glob(os.path.join(path_data, '*crop.fits')):
        
        hdu = fits.open(input_file)
        oim = hdu[0].data
        xs, ys = oim.shape

        head = hdu[0].header

        nf = input_file.split('/')[-1][:-5]
        nfp = os.path.join(path_wavelets, nf)
        print(input_file)
        print(nfp)
        output = synthesis_fullfield( oim, nfp, xs, ys, n_levels, write_fits = True )
        output = synthesis_bcgwavsizesep_with_masks( nfp, lvl_sep, lvl_sep_max, lvl_sep_bcg, size_sep_pix, xs, ys, rc_pix, n_levels, mscell, mscbcg, kurt_filt = False, write_fits = True)