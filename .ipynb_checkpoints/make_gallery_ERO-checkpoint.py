#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 18:16:13 2024

@author: aellien
"""

import sys
sys.path.append("/home/aellien/dawis")
import glob as glob
import os
import gc
import numpy as np
import pyregion as pyr
import random
import pandas as pd
import ray
import dawis as d
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import kurtosis
from photutils.background import Background2D, MedianBackground
from mpl_toolkits.axes_grid1 import make_axes_locatable


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def gallery_synthesis( nfp, path_plot, oim, lvl_sep, lvl_sep_max, lvl_sep_bcg, size_sep_pix, xs, ys, rc_pix, n_levels, mscell, mscbcg, mscstar, kurt_filt = True, write_fits = True ):
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
    res = np.copy(oim)
    rec = np.zeros( (xs, ys) )
    star = np.zeros( (xs, ys) )
    im_art = np.zeros( (xs, ys) )
    im_unclass = np.zeros( (xs, ys) )
    icl_err = np.zeros( (xs, ys) )

    icl_al = []
    gal_al = []
    noticl_al = []
    unclass_al = []

    #%
    at_test = []
    #%

    xc = xs / 2.
    yc = ys / 2.
    
    tn_std = d.table_noise_stds(11,'BSPL')

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
        
        gc.collect()
        if i!=118:
            print(i)
            ol = d.read_objects_from_pickle( op )
            itl = d.read_interscale_trees_from_pickle( itlp )
            ############################################################## MEMORY
    
            atom = np.zeros((xs, ys))
            
            if i!=0:
                # Kurtosis + ICL+BCG
                for j, o in enumerate(ol):
            
                    x_min, y_min, x_max, y_max = o.bbox
                    sx = x_max - x_min
                    sy = y_max - y_min
                    itm = itl[j].interscale_maximum
                    xco = itm.x_max
                    yco = itm.y_max
                    lvlo = o.level
                    
                    rec[ x_min : x_max, y_min : y_max ] += o.image
                    atom[ x_min : x_max, y_min : y_max ] += o.image


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
                                print(o.level, tn_std[o.level] * o.image.size, np.sum(o.image))
                                
                        
                                icl_err[ x_min : x_max, y_min : y_max ] += tn_std[o.level] * o.image.astype(bool)
            
                            elif o.level > 3 :
                                
                                print(o.level, tn_std[o.level] * o.image.size, np.sum(o.image))
                                icl[ x_min : x_max, y_min : y_max ] += o.image * o.image.astype(bool)
                                icl_err[ x_min : x_max, y_min : y_max ] += tn_std[o.level]
                                icl_al.append([o, xco, yco])
                                
                            else:
                                gal[ x_min : x_max, y_min : y_max ] += o.image
            
                        # ICL
                        else:
            
                            if (o.level >= lvl_sep) & (sx >= size_sep_pix) & (sy >= size_sep_pix):
                                
                                print(o.level, tn_std[o.level] * o.image.size, np.sum(o.image))
                                icl[ x_min : x_max, y_min : y_max ] += o.image * o.image.astype(bool)
                                icl_err[ x_min : x_max, y_min : y_max ] += tn_std[o.level]
                                icl_al.append([o, xco, yco])
                                at_test.append([xco, yco])
                        
                            elif mscstar[xco, yco] == 1:
                                star[ x_min : x_max, y_min : y_max ] += o.image
                    
                            else:
                                noticl_al.append([ o, xco, yco ])
                                gal[ x_min : x_max, y_min : y_max ] += o.image
            
                    elif mscstar[xco, yco] == 1:
                        star[ x_min : x_max, y_min : y_max ] += o.image
            
                    else:
                        noticl_al.append([ o, xco, yco ])
                        gal[ x_min : x_max, y_min : y_max ] += o.image
                
                res -= atom
            else:continue
            
            fig, axes = plt.subplots(2, 3)
            
            
            norm1 = ImageNormalize(res, interval = ZScaleInterval(), stretch = LinearStretch())
            cmap = 'binary'
            
            ax = axes[0][0]
            axim = ax.imshow(oim, norm = norm1, cmap = cmap, origin = 'lower')
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.0)
            caxre = fig.colorbar(axim, cax = cax, aspect = 10, pad = 0, \
                                        orientation = 'vertical')
            caxre.ax.tick_params(labelsize = 5)
            ax.set_title('Original Image', fontsize = 8)
            
            #
            ax = axes[0][1]
            axim = ax.imshow(rec, norm = norm1, cmap = cmap, origin = 'lower')
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.0)
            caxre = fig.colorbar(axim, cax = cax, aspect = 10, pad = 0, \
                                        orientation = 'vertical')
            caxre.ax.tick_params(labelsize = 5)
            ax.set_title('Reconstructed Image', fontsize = 8)
            
            #
            ax = axes[0][2]
            ax.imshow(res, norm = norm1, cmap = cmap, origin = 'lower')
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.0)
            caxre = fig.colorbar(axim, cax = cax, aspect = 10, pad = 0, \
                                        orientation = 'vertical')
            caxre.ax.tick_params(labelsize = 5)
            ax.set_title('Residuals', fontsize = 8)
            
            #
            ax = axes[1][0]
            norm = ImageNormalize(star, interval = ZScaleInterval(), stretch = LinearStretch())
            ax.imshow(star, norm = norm1, cmap = cmap, origin = 'lower')
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.0)
            caxre = fig.colorbar(axim, cax = cax, aspect = 10, pad = 0, \
                                        orientation = 'vertical')
            caxre.ax.tick_params(labelsize = 5)
            ax.set_title('Foreground MW stars', fontsize = 8)
            
            #
            ax = axes[1][1]
            norm = ImageNormalize(gal, interval = ZScaleInterval(), stretch = LinearStretch())
            ax.imshow(gal, norm = norm1, cmap = cmap, origin = 'lower')
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.0)
            caxre = fig.colorbar(axim, cax = cax, aspect = 10, pad = 0, \
                                        orientation = 'vertical')
            caxre.ax.tick_params(labelsize = 5)
            ax.set_title('Galaxies', fontsize = 8)
            
            #
            ax = axes[1][2]
            norm = ImageNormalize(icl, interval = ZScaleInterval(), stretch = LinearStretch())
            ax.imshow(icl, norm = norm1, cmap = cmap, origin = 'lower')
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.0)
            caxre = fig.colorbar(axim, cax = cax, aspect = 10, pad = 0, \
                                        orientation = 'vertical')
            caxre.ax.tick_params(labelsize = 5)
            ax.set_title('ICL model', fontsize = 8)
            
            #plt.tight_layout()
            plt.suptitle('Iteration %03d'%i)
            plt.subplots_adjust( left=0.05, bottom=0.05, right=0.93, top=0.98, wspace=0.18, hspace=0.01)
            plt.savefig(os.path.join(path_plots, 'ERO_H_gallery.it%03d.png'%i), dpi = 1000, format = 'png')
            plt.close()
        else:
            pass
        
    if write_fits == True:
        print('\nWS + SF + SS -- ICL+BCG -- write fits as %s*'%(nfp))

        # write to fits
        hduo = fits.PrimaryHDU(icl)
        hduo.writeto( nfp + 'synth.icl.bcgwavsizesepmask_%03d_%03d.fits'%(lvl_sep, size_sep_pix), overwrite = True )

        hduo = fits.PrimaryHDU(gal)
        hduo.writeto( nfp + 'synth.gal.bcgwavsizesepmask_%03d_%03d.fits'%(lvl_sep, size_sep_pix), overwrite = True )

        hduo = fits.PrimaryHDU(icl_err)
        hduo.writeto( nfp + 'synth.icl_err.bcgwavsizesepmask_%03d_%03d.fits'%(lvl_sep, size_sep_pix), overwrite = True )

    return icl, gal


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/home/aellien/Euclid_ERO/data/Euclid-NISP-Stack-ERO-Abell2390.DR3'
    path_scripts = '/home/aellien/Euclid_ERO/Euclid_ERO_scripts'
    path_wavelets = '/home/aellien/Euclid_ERO/wavelets/out4/'
    path_plots = '/home/aellien/Euclid_ERO/plots'
    path_analysis = '/home/aellien/Euclid_ERO/analysis/'
    
    n_levels = 11
    lvl_sep = 6
    lvl_sep_max = 12
    lvl_sep_bcg = 5
    size_sep_pix = 10
    rc_pix = 10

    mscell = fits.getdata(os.path.join(path_analysis,'mscell.fits'))
    mscbcg = fits.getdata(os.path.join(path_analysis,'mscbcg.fits'))
    mscstar = fits.getdata(os.path.join(path_analysis,'mscstar.fits'))

    for input_file in glob.glob(os.path.join(path_data, '*crop.fits')):
        
        hdu = fits.open(input_file)
        oim = hdu[0].data
        xs, ys = oim.shape
        nx, ny = int(xs / 25.), int(ys / 25.)
        head = hdu[0].header
        bkg_estimator = MedianBackground()
        bkg = Background2D(oim, (nx, ny), filter_size = (3, 3), bkg_estimator = bkg_estimator)
        oim -= bkg.background_median

        nf = input_file.split('/')[-1][:-5]
        nfp = os.path.join(path_wavelets, nf)
        print(input_file)
        print(nfp)
        output = gallery_synthesis( nfp, path_plots, oim, lvl_sep, lvl_sep_max, lvl_sep_bcg, size_sep_pix, xs, ys, rc_pix, n_levels, mscell, mscbcg, mscstar, kurt_filt = False, write_fits = True)