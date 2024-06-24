#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 14:31:33 2024

@author: aellien
"""

import sys
import dawis as d
import glob as glob
import os
import numpy as np
import pyregion as pyr
import random
import pandas as pd
import ray
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.morphology import binary_dilation
from scipy.stats import kurtosis
import gc
from power_ratio import *
from datetime import datetime
import h5py

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def selection_error(atom_in_list, atom_out_list, M, percent, xs, ys, flux_lim, mscsedl):
    '''Computation of classification error on flux.
    '''
    # Output array
    sed_sample = []

    # Sample
    size_sample = np.random.uniform(low = int( len(atom_in_list) * (1. - percent)), \
                               high = int( len(atom_in_list) + len(atom_in_list) * percent ), \
                               size = M).astype(int)
    replace_sample = []
    for s in size_sample:
        replace_sample.append(int(np.random.uniform(low = 0, high = int( s * percent ))))
    replace_sample = np.array(replace_sample)

    flux_sample = []
    flux_err_sample = []
    for i, (s, r) in enumerate(zip(size_sample, replace_sample)):
        
        #print(i, s, r, len(atom_in_list), len(atom_out_list))
        im_s = np.zeros((xs, ys))
        im_s_err = np.zeros((xs, ys))
        if s < len(atom_in_list):
            #print('cc1')
            flux = 0
            draw = random.sample(atom_in_list, s)

        if s >= len(atom_in_list):
            #print('cc2')

            flux = 0
            draw1 = random.sample(atom_in_list, len(atom_in_list) - r)
            draw2 = random.sample(atom_out_list, s - len(atom_in_list) + r)
            draw = draw1 + draw2
            
        #print(i, len(atom_in_list), len(atom_out_list), len(draw), len(draw[0]), draw[0])
        for (image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo) in draw:
            im_s[ x_min : x_max, y_min : y_max ] += image
            im_s_err[ x_min : x_max, y_min : y_max ] += det_err_image
            #flux += np.sum(o.image)

        flux = np.sum(im_s[im_s >= flux_lim])
        flux_err = np.sqrt(np.sum(im_s_err[im_s >= flux_lim]**2))
        flux_sample.append(flux)
        flux_err_sample.append(flux_err)

        line = []
        for mscsed in mscsedl:
            flux_sed = np.sum(im_s[mscsed.astype(bool)])
            line.append(flux_sed)

        sed_sample.append(line)

    sed_sample = np.array(sed_sample)
    mean_sed = np.median(sed_sample, axis = 0)
    up_err_sed = np.percentile(sed_sample, 95, axis = 0)
    low_err_sed = np.percentile(sed_sample, 5, axis = 0)
    out_sed = np.array([ mean_sed, low_err_sed, up_err_sed ] ).swapaxes(0, 1).flatten() # 1d size n_sed_region x 3

    flux_sample = np.array(flux_sample)
    flux_err_sample = np.array(flux_err_sample)
    mean_flux = np.median(flux_sample)
    up_flux = np.percentile(flux_sample, 95)
    low_flux = np.percentile(flux_sample, 5)
    mean_flux_err = np.median(flux_sample)
    up_flux_err = np.percentile(flux_sample, 95)
    low_flux_err = np.percentile(flux_sample, 5)

    #plt.figure()
    #plt.hist(flux_sample, bins = 10)
    #plt.show()

    return mean_flux, low_flux, up_flux, mean_flux_err, up_flux_err, low_flux_err, out_sed

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def PR_with_selection_error(atom_in_list, atom_out_list, M, percent, R, xs, ys):
    '''Computation of classification error on PR.
    '''
    size_sample = np.random.uniform(low = int( len(atom_in_list) * (1. - percent)), \
                               high = int( len(atom_in_list) + len(atom_in_list) * percent ), \
                               size = M).astype(int)
    replace_sample = []
    for s in size_sample:
        replace_sample.append(int(np.random.uniform(low = 0, high = int( s * percent ))))
    replace_sample = np.array(replace_sample)

    PR_sample = []
    for s, r in zip(size_sample, replace_sample):
        im = np.zeros((xs, ys))

        if s < len(atom_in_list):

            draw = random.sample(atom_in_list, s)
            for (image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo) in draw:
                im[ x_min : x_max, y_min : y_max ] += image

            orderl = []
            for order in range(1, 5):
                PR = power_ratio( image = im, order = order, radius = R )
                orderl.append(PR)
            PR_sample.append(orderl)

        if s >= len(atom_in_list):

            flux = 0
            draw1 = random.sample(atom_in_list, len(atom_in_list) - r)
            draw2 = random.sample(atom_out_list, s - len(atom_in_list) + r)
            draw = draw1 + draw2
            for (image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo) in draw:
                im[ x_min : x_max, y_min : y_max ] += image

            orderl = []

            for order in range(1, 5):
                PR = power_ratio( image = im, order = order, radius = R )
                orderl.append(PR)
            PR_sample.append(orderl)

        #%---
        #interval = AsymmetricPercentileInterval(5, 99.5) # meilleur rendu que MinMax or ZScale pour images reconstruites
        #fig, ax = plt.subplots(1)
        #poim = ax.imshow(im, norm = ImageNormalize( im, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')
        #%---

    PR_sample = np.array(PR_sample)
    PR_results = []
    for i in range(1, 5):
        mean_PR = np.median(PR_sample[:, i - 1])
        up_err = np.percentile(PR_sample[:, i - 1], 95)
        low_err = np.percentile(PR_sample[:, i - 1], 5)
        PR_results.append([mean_PR, up_err, low_err])
    #%---
    #plt.figure()
    #for i in range(1,5):
    #    plt.hist(PR_sample[:, i-1], bins = 10, alpha = 0.5)
    #plt.show()
    #%---

    return PR_results

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def synthesis_bcgwavsizesep_with_masks( nfp, lvl_sep, lvl_sep_max, lvl_sep_bcg, size_sep, size_sep_pix, xs, ys, n_levels, mscstar, mscell, mscbcg, mscsedl, msat, R, rc_pix, N_err, per_err, flux_lim, kurt_filt = True, plot_vignet = False, write_fits = True, measure_PR = False ):
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
    icl_dei = np.zeros( (xs, ys) )
    gal = np.zeros( (xs, ys) )
    gal_dei = np.zeros( (xs, ys) )
    im_art = np.zeros( (xs, ys) )
    im_unclass = np.zeros( (xs, ys) )
    im_unclass_dei = np.zeros( (xs, ys) )

    tot_icl_al = []
    tot_gal_al = []
    #tot_noticl_al = []
    tot_unclass_al = []

    #%
    at_test = []
    #%

    xc = xs / 2.
    yc = ys / 2.

    ######################################## MEMORY v
    opath = nfp + '*ol.it*.hdf5'
    itpath = nfp + '*itl.it*.hdf5'
    opathl = glob.glob(opath)
    opathl.sort()

    # Interscale tree lists
    itpathl = glob.glob(itpath)
    itpathl.sort()
    
    for i, ( op, itlp ) in enumerate( zip( opathl, itpathl )):
        print('Iteration %d' %(i), end ='\r')
        
        #ol = d.store_objects.read_ol_from_hdf5(op)
        #itl = d.store_objects.read_itl_from_hdf5(itlp)
        with h5py.File(op, "r") as f1, h5py.File(itlp, "r") as f2:
            gc.collect()
            icl_al = []
            gal_al = []
            noticl_al = []
            unclass_al = []
            for o, it in zip(f1.keys(), f2.keys()):
                ######################################## MEMORY ^

                x_min, y_min, x_max, y_max = f1[o]['bbox'][()]
                image = f1[o]['image'][()]
                det_err_image = f1[o]['det_err_image'][()]
                itm = f2[it]['interscale_maximum']
                xco = itm['x_max'][()]
                yco = itm['y_max'][()]
                lvlo = f1[o]['level'][()]
                sx = x_max - x_min
                sy = y_max - y_min
            
                if kurt_filt == True:
                    k = kurtosis(image.flatten(), fisher=True)
                    if k < 0:
                        im_art[ x_min : x_max, y_min : y_max ] += image
                        continue
        
                # Remove background
                if lvlo >= lvl_sep_max:
                    tot_unclass_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
                    continue
        
                # Only atoms within analysis radius
                dR = np.sqrt( (xc - xco)**2 + (yc - yco)**2 )
                if dR > R:
                    continue
    
                # ICL + BCG
                if (mscstar[xco, yco] != 1) & (mscell[xco, yco] == 1):
    
                    '''# BCG
                    xbcg, ybcg = [ xs, ys ] # pix long, ds9 convention
                    if mscbcg[xco, yco] == 1:
    
                        dr = np.sqrt( (xbcg - xco)**2 + (ybcg - yco)**2 )
                        if (o.level <= 3) & (dr < rc_pix):
    
                            icl[ x_min : x_max, y_min : y_max ] += o.image
                            icl_al.append([o, xco, yco])
    
                        elif o.level >= 4:
                            icl[ x_min : x_max, y_min : y_max ] += o.image
                            icl_al.append([o, xco, yco])'''
                            
                    # BCG
                    if mscbcg[xco, yco] == 1:
                        icl[ x_min : x_max, y_min : y_max ] += image
                        icl_dei[ x_min : x_max, y_min : y_max ] += det_err_image
                        icl_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
                        tot_icl_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
                        continue
    
                    # ICL
                    if (lvlo >= lvl_sep) & (sx >= size_sep_pix) & (sy >= size_sep_pix):

                        #%%%%% Je laisse au cas où %%%%% v
                        coo_spur_halo = []
                        # [ [1615, 1665], [1685, 1480], [530, 260] ] # pix long, ds9 convention
                        flag = False
                        for ygal, xgal in coo_spur_halo:

                            dr = np.sqrt( (xgal - xco)**2 + (ygal - yco)**2 )
                            if (dr <= rc_pix) & (lvlo == 5):
                                flag = True
                        #%%%%%%% ^^^^^^^^^^^^^^^^^^^^^^^^^^
                            
                        if flag == False:
                            icl[ x_min : x_max, y_min : y_max ] += image
                            icl_dei[ x_min : x_max, y_min : y_max ] += det_err_image
                            icl_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
                            tot_icl_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
                            at_test.append([xco, yco])
                            continue
                            
                    else:
                        noticl_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
                else:
                    noticl_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
                    
            # Galaxies
            for j, (image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo) in enumerate(noticl_al):
                
                # Satellites
                if (mscstar[xco, yco] != 1) & (lvlo < lvl_sep) & (msat[xco, yco] == 1):

                    gal[ x_min : x_max, y_min : y_max ] += image
                    gal_dei[ x_min : x_max, y_min : y_max ] += det_err_image
                    gal_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
                    tot_gal_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
                    
                # If not identified as galaxies --> test if BCG again
                else:
                    unclass_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])

            # Test for unclassified atoms --> sometimes extended BCG halo is missed because
            # of the nature of wavsep scheme.
            for j, (image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo) in enumerate(unclass_al):
                
                # Case in which it is possible that it is BCG halo?
                if (lvl_sep > lvl_sep_bcg) & (lvlo >= lvl_sep_bcg) & (mscell[xco, yco] == 1) :
                    icl[ x_min : x_max, y_min : y_max ] += image
                    icl_dei[ x_min : x_max, y_min : y_max ] += det_err_image
                    icl_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
    
                #If not --> unclassified
                else:
                    im_unclass[ x_min : x_max, y_min : y_max ] += image
                    im_unclass_dei[ x_min : x_max, y_min : y_max ] += det_err_image
                    tot_unclass_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])

    # Remove potential foreground star artifacts
    #gal[mscstar == 1.] = 0.
    
    at_test = np.array(at_test)

    if write_fits == True:
        print('\nWS + SF + SS -- ICL+BCG -- LVL_SEP = %d -- write fits as %s*'%(lvl_sep, nfp))

        hdu = fits.PrimaryHDU()
        hdu_icl = fits.ImageHDU(icl, name = 'ICL+BCG')
        hdu_gal = fits.ImageHDU(gal, name = 'SATELLITES')
        tot = gal + icl
        hdu_tot = fits.ImageHDU(tot, name = 'ICL+BCG+SATELLITES')
        hdu_icl_dei = fits.ImageHDU(icl_dei, name = 'ICL+BCG DET. ERR.')
        hdu_gal_dei = fits.ImageHDU(gal_dei, name = 'SAT DET. ERR.')
        tot_dei = icl_dei + gal_dei
        hdu_tot_dei = fits.ImageHDU(tot_dei, name = 'ICL+BCG+SAT DET. ERR.')
        hdu_unclass = fits.ImageHDU(im_unclass, name = 'UNCLASSIFIED')
        hdu_unclass_dei = fits.ImageHDU(im_unclass, name = 'UNCLASSIFIED DET. ERR.')

        
        hdul = fits.HDUList([ hdu, hdu_icl, hdu_gal, hdu_tot, hdu_unclass, hdu_icl_dei, hdu_gal_dei, hdu_tot_dei, hdu_unclass_dei ])
        hdul.writeto( nfp + 'synth.bcgwavsizesepmask_%03d_%03d.fits'%(lvl_sep, size_sep), overwrite = True )

    # Plot vignets
    if plot_vignet == True:

        interval = AsymmetricPercentileInterval(5, 99.5) # meilleur rendu que MinMax or ZScale pour images reconstruites
        fig, ax = plt.subplots(2, 2)
        poim = ax[0][0].imshow(gal, norm = ImageNormalize( gal, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')
        poim = ax[1][0].imshow(icl, norm = ImageNormalize( icl, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')
        #%
        rco = pyr.open(os.path.join(path_data, 'star_flags_polygon_pix_long.reg'))
        patch_list, artist_list = rco.get_mpl_patches_texts()
        for p in patch_list:
            ax[1][0].add_patch(p)
        for a in artist_list:
            ax[1][0].add_artist(a)
        for at in at_test:
            ax[1][0].plot(at[1], at[0], 'b+')
        #%
        poim = ax[0][1].imshow(im_unclass, norm = ImageNormalize( gal, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')
        poim = ax[1][1].imshow(im_art, norm = ImageNormalize( gal, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')

        #plt.show()
        plt.tight_layout()
        plt.savefig( nfp + 'results.bcgwavsizesepmask_%03d_%03d_testspur.png'%(lvl_sep, size_sep), format = 'png' )
        print('Write vignet to' + nfp + 'synth.bcgwavsizesepmask_%03d_%03d_testspur.png'%(lvl_sep, size_sep))
        plt.close('all')

    if measure_PR == True:
        print('start bootstrap')
        start = datetime.now()
        # Measure Fractions and uncertainties
        F_ICL_m, F_ICL_low, F_ICL_up, FICL_det_err, low_FICL_det_err, up_FICL_det_err, out_sed_icl =  selection_error(tot_icl_al, tot_unclass_al+tot_gal_al, M = N_err, percent = per_err, xs = xs, ys = ys, flux_lim = flux_lim, mscsedl = mscsedl)
        F_gal_m, F_gal_low, F_gal_up, Fgal_det_err, low_Fgal_det_err, up_Fgal_det_err, out_sed_gal =  selection_error(tot_gal_al, tot_unclass_al+tot_icl_al, M = N_err, percent = per_err, xs = xs, ys = ys, flux_lim = flux_lim, mscsedl = mscsedl)
        f_ICL_m = F_ICL_m / (F_ICL_m + F_gal_m)
        f_ICL_low = F_ICL_low / (F_ICL_low + F_gal_up)
        f_ICL_up = F_ICL_up / (F_ICL_up + F_gal_low)
        print(datetime.now() - start)
        
        print('\nWS + SF -- ICL+BCG --  z = %d'%lvl_sep)
        print('N = %4d   F_ICL = %f Mjy/sr  err_low = %f Mjy/sr  err_up = %f Mjy/sr'%(len(tot_icl_al), F_ICL_m, F_ICL_low, F_ICL_up))
        print('N = %4d   F_gal = %f Mjy/sr  err_low = %f Mjy/sr  err_up = %f Mjy/sr'%(len(tot_gal_al), F_gal_m, F_gal_low, F_gal_up))
        print('Det. error: deICL = %f ADU  deICL_low = %f  deICLup = %f ADU'%(FICL_det_err, low_FICL_det_err, up_FICL_det_err))
        print('f_ICL = %1.3f    f_ICL_low = %1.3f   f_ICL_up = %1.3f'%(f_ICL_m, f_ICL_low, f_ICL_up))
        
        # Measure Power ratio
        results_PR = PR_with_selection_error(atom_in_list = tot_icl_al, atom_out_list = tot_unclass_al+tot_gal_al, M = N_err, percent = per_err, R = R, xs = xs, ys = ys)
        PR_1_m, PR_1_up, PR_1_low = results_PR[0]
        PR_2_m, PR_2_up, PR_2_low = results_PR[1]
        PR_3_m, PR_3_up, PR_3_low = results_PR[2]
        PR_4_m, PR_4_up, PR_4_low = results_PR[3]
        
        print('PR_1_m = %1.2e    PR_1_low = %1.2e    PR_1_up = %1.2e'%(PR_1_m, PR_1_low, PR_1_up))
        print('PR_2_m = %1.2e    PR_2_low = %1.2e    PR_2_up = %1.2e'%(PR_2_m, PR_2_low, PR_2_up))
        print('PR_3_m = %1.2e    PR_3_low = %1.2e    PR_3_up = %1.2e'%(PR_3_m, PR_3_low, PR_3_up))
        print('PR_4_m = %1.2e    PR_4_low = %1.2e    PR_4_up = %1.2e'%(PR_4_m, PR_4_low, PR_4_up))
    
        out_sed_icl_df = pd.DataFrame( out_sed_icl, columns = [ 'icl_reg_%d_%s'%(i/3, hkw[i%3]) for i in range(len(output[-1]))]) # create df with all SED flux for all regions with correctly numbered column names
        out_sed_gal_df = pd.DataFrame( out_sed_gal, columns = [ 'gal_reg_%d_%s'%(i/3, hkw[i%3]) for i in range(len(output[-1]))]) # create df with all SED flux for all regions with correctly numbered column names

        output_df = pd.DataFrame( [[ nf, filt, sch, R_kpc, R_pix, lvl_sep, size_sep, F_ICL_m, F_ICL_low, F_ICL_up, F_gal_m, F_gal_low, F_gal_up, f_ICL_m, f_ICL_low, f_ICL_up, PR_1_m, PR_1_up, PR_1_low, PR_2_m, PR_2_up, PR_2_low, PR_3_m, PR_3_up, PR_3_low, PR_4_m, PR_4_up, PR_4_low ]], \
                        columns = [ 'nf', 'filter', 'Atom selection scheme', 'R_kpc', 'R_pix', 'lvl_sep', 'size_sep','F_ICL_m', 'F_ICL_low', 'F_ICL_up', 'F_gal_m', 'F_gal_low', 'F_gal_up', 'f_ICL_m', 'f_ICL_low', 'f_ICL_up', 'PR_1_m', 'PR_1_up', 'PR_1_low', 'PR_2_m', 'PR_2_up', 'PR_2_low', 'PR_3_m', 'PR_3_up', 'PR_3_low', 'PR_4_m', 'PR_4_up', 'PR_4_low'  ])
        output_df = pd.concat( [output_df, out_sed_icl_df, out_sed_gal_df], axis = 1)
        
        return output_df

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if __name__ == '__main__':

    # Paths, lists & variables
    
    path_data = '/n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3/'
    path_wavelets = '/n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3/wavelets/out7'
    path_analysis = path_data
    
    '''# Paths, lists & variables
    path_data = '/home/aellien/Euclid_ERO/data/Euclid-NISP-Stack-ERO-Abell2390.DR3'
    path_scripts = '/home/aellien/Euclid_ERO/Euclid_ERO_scripts'
    path_wavelets = '/home/aellien/Euclid_ERO/wavelets/out6/'
    path_plots = '/home/aellien/Euclid_ERO/plots'
    path_analysis = '/home/aellien/Euclid_ERO/analysis/'
    '''
    
        
    nf = sys.argv[1] #[ 'Euclid-NISP-J-ERO-Abell2390-LSB.DR3.crop.fits', 'Euclid-NISP-H-ERO-Abell2390-LSB.DR3.crop.fits', 'Euclid-NISP-Y-ERO-Abell2390-LSB.DR3.crop.fits' ]
    nfp = os.path.join( path_wavelets, nf[:-4] )
    oim_file = os.path.join( path_data, nf )
    hdu = fits.open(oim_file)
    oim = hdu[0].data
    xs, ys = oim.shape
    
    lvl_sep = int(sys.argv[2]) # lvl_sepl = [ 4, 5, 6 ] # wavelet scale separation
    size_sep = int(sys.argv[3]) # size_sepl = [60, 80, 100, 140, 200 ] # size separation [kpc]
    R_kpc = int(sys.argv[4]) # R_kpcl = [ 1000 ] # radius in which quantities are measured [kpc]
    lvl_sep_bcg = int(sys.argv[5]) # lvl_sep_bcg = 6
    n_levels = int(sys.argv[6]) # n_levels = 10
    lvl_sep_max = int(sys.argv[7]) # lvl_sep_max = 9
    rc = int(sys.argv[8]) # rc = 10 # kpc, distance to center to be classified as gal
    
    counter=int(os.environ['SLURM_PROCID'])
    print(f"hello from worker number {counter} - {nf}")
    
    physcale = 3.68 # kpc/"
    pix_scale = 0.3
    mu_lim = 35
    ZP_AB = 30
    flux_lim = 10**( (ZP_AB - mu_lim) / 2.5 )
    rc_pix = rc / physcale / pix_scale # pixels
    size_sep_pix = size_sep / physcale / pix_scale # pixels
    R_pix = R_kpc / physcale / pix_scale # pixels
    
    N_err = 100
    per_err = 0.05

    kurt_filt = True
    plot_vignet = False
    write_fits = True
    measure_PR = True
    write_dataframe = True
    
    # Read galaxy catalog
    # cat_gal = np.loadtxt(os.path.join(path_analysis,'A2390_redMapper_Pmem_gt_0.8.txt'))[:, -2:]

    # Masks
    mscell = fits.getdata(os.path.join(path_analysis,'mscell.fits'))
    mscbcg = fits.getdata(os.path.join(path_analysis,'mscbcg.fits'))
    mscstar = fits.getdata(os.path.join(path_analysis,'mscstar.fits'))
    msat = fits.getdata(os.path.join(path_analysis,'msat.fits'))
    mscsedl = []
    
    # Synthesis
    output_df = synthesis_bcgwavsizesep_with_masks( nfp = nfp, 
                                             lvl_sep = lvl_sep, 
                                             lvl_sep_max = lvl_sep_max, 
                                             lvl_sep_bcg = lvl_sep_bcg, 
                                             size_sep = size_sep, 
                                             size_sep_pix =  size_sep_pix, 
                                             xs = xs, 
                                             ys = ys, 
                                             n_levels = n_levels,
                                             mscstar = mscstar, 
                                             mscell = mscell, 
                                             mscbcg = mscbcg, 
                                             mscsedl = mscsedl,
                                             msat = msat, 
                                             R = R_pix, 
                                             rc_pix = rc_pix,
                                             N_err = N_err, 
                                             per_err = per_err, 
                                             flux_lim = flux_lim, 
                                             kurt_filt = kurt_filt, 
                                             plot_vignet = plot_vignet, 
                                             write_fits = write_fits, 
                                             measure_PR = measure_PR )
    
    ofp = nfp + 'df.csv'
    print('Write results to %s'%ofp)
    results_df.to_csv(ofp)
