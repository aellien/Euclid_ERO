#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 23:44:18 2024

@author: aellien
"""
import os as os
import glob as glob
import dawis as d
import numpy as np
import pandas as pd
import random
import gc
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import *
from scipy.stats import kurtosis
import h5py
from cosmo_calc import cosmo_calc
from photutils.segmentation import SourceCatalog, detect_sources

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
def selection_error(atom_in_list, atom_out_list, M, percent, xs, ys, mscann):
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
    for s, r in zip(size_sample, replace_sample):

        im_s = np.zeros((xs, ys))
        if s < len(atom_in_list):
            flux = 0
            draw = random.sample(atom_in_list, s)

        if s >= len(atom_in_list):
            flux = 0
            draw1 = random.sample(atom_in_list, len(atom_in_list) - r)
            draw2 = random.sample(atom_out_list, s - len(atom_in_list) + r)
            draw = draw1 + draw2

        for (image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo) in draw:            
            im_s[ x_min : x_max, y_min : y_max ] += image
            #flux += np.sum(o.image)

        flux = np.sum(im_s[mscann.astype(bool)])
        flux_sample.append(flux)

    
    flux_sample = np.array(flux_sample)
    mean_flux = np.median(flux_sample)
    up_err = np.percentile(flux_sample, 95)
    low_err = np.percentile(flux_sample, 5)
    

    #plt.figure()
    #plt.hist(flux_sample, bins = 10)
    #plt.show()

    return flux_sample, mean_flux, low_err, up_err


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def synthesis_bcgwavsizesep_with_masks( nfwp, nfap, lvl_sep, lvl_sep_max, lvl_sep_bcg, size_sep, size_sep_pix, xs, ys, n_levels, mscicl, mscbcg, mscann, N_err, per_err, kurt_filt = True, plot_vignet = False, write_fits = True ):
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
    recim = np.zeros( (xs, ys) )

    tot_icl_al = []
    tot_gal_al = []
    #tot_noticl_al = []
    icl_boot_al = []

    mscell = mscicl
    mscstar = np.zeros((xs, ys))
    msat = np.zeros((xs, ys))
    
    #%
    at_test = []
    #%

    xc = xs / 2.
    yc = ys / 2.

    ######################################## MEMORY v
    opath = nfwp + '*ol.it*.hdf5'
    opathl = glob.glob(opath)
    opathl.sort()
    memory = []
    marea = []
    for i, op in enumerate(opathl):
        
        print('Iteration %d' %(i))#, end ='\r')
        #snapshot = tracemalloc.take_snapshot()
        #top_stats = snapshot.statistics('lineno')
        #print(top_stats[0])
        #memory.append(top_stats[0].size/1e9)
        
        
        #ol = d.store_objects.read_ol_from_hdf5(op)
        #itl = d.store_objects.read_itl_from_hdf5(itlp)
        with h5py.File(op, "r") as f1:
            gc.collect()
            icl_al = []
            gal_al = []
            noticl_al = []
            unclass_al = []
            areal = []
            for o in f1.keys():

                x_min, y_min, x_max, y_max = np.copy(f1[o]['bbox'][()])
                image = np.copy(f1[o]['image'][()])
                det_err_image = np.copy(f1[o]['det_err_image'][()])
                lvlo = np.copy(f1[o]['level'][()])
                ######################################## MEMORY ^
            
                sx = x_max - x_min
                sy = y_max - y_min
                areal.append(sx*sy)
                m = detect_sources(image, threshold = 0., npixels=1)
                c = SourceCatalog(image, m)
                xco = int(c.centroid_quad[0][1] + x_min)
                yco = int(c.centroid_quad[0][0] + y_min)
            
                if kurt_filt == True:
                    k = kurtosis(image.flatten(), fisher=True)
                    if k < 0:
                        im_art[ x_min : x_max, y_min : y_max ] += image
                        gc.collect()
                        continue
        
                # Remove background, and add it to bootstrap if center is within ICL ellipse
                if lvlo >= lvl_sep_max:
                    if mscell[xco, yco] == 1:
                        icl_boot_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
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
                        coo_spur_halo =  [] # pix long, ds9 convention
                        flag = False
                        for ygal, xgal in coo_spur_halo:

                            dr = np.sqrt( (xgal - xco)**2 + (ygal - yco)**2 )
                            if (dr <= 100) & (4 < lvlo < 7 ):
                                flag = True
                        #%%%%%%% ^^^^^^^^^^^^^^^^^^^^^^^^^^
                            
                        if flag == False:
                            icl[ x_min : x_max, y_min : y_max ] += image
                            icl_dei[ x_min : x_max, y_min : y_max ] += det_err_image
                            icl_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
                            tot_icl_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
                            at_test.append([xco, yco])
                            continue
                        
                        noticl_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
                            
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

                    # add 'tendencious' galaxy atoms to list for bootstrap
                    if ( lvlo == (lvl_sep - 1)) & (mscell[xco, yco] == 1):
                        icl_boot_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo]) 
                        
                # If not identified as galaxies --> test if BCG again
                else:
                    unclass_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])

            # Test for unclassified atoms --> sometimes extended BCG halo is missed because
            # of the nature of wavsep scheme.
            for j, (image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo) in enumerate(unclass_al):
                
                # Case in which it is possible that it is BCG halo?
                if (lvl_sep > lvl_sep_bcg) & (lvlo >= lvl_sep_bcg) & (np.sqrt( (xc - xco)**2 + (yc - yco)**2 ) < 100) :
                    icl[ x_min : x_max, y_min : y_max ] += image
                    icl_dei[ x_min : x_max, y_min : y_max ] += det_err_image
                    icl_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
    
                #If not --> unclassified 
                else:
                    im_unclass[ x_min : x_max, y_min : y_max ] += image
                    im_unclass_dei[ x_min : x_max, y_min : y_max ] += det_err_image
                    
                    #  add tendencious atoms to list for bootstrap (arbitrary conditions)
                    if ( lvlo == (lvl_sep - 1)) & (mscell[xco, yco] == 1):
                        icl_boot_al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
                        
            marea.append(np.mean(areal))
            print(len(tot_icl_al), len(tot_gal_al), len(icl_boot_al))

            recim[ x_min : x_max, y_min : y_max ] += image
            
    # clear some memory
    gc.collect()
    icl_al.clear()
    gal_al.clear()
    noticl_al.clear()
    unclass_al.clear()

    if write_fits == True:
        print('\nWS + SF + SS -- ICL+BCG -- write fits as %s*'%(nfap))

        # write to fits
        hduo = fits.PrimaryHDU(icl)
        hduo.writeto( nfap + '.synth.icl.bcgwavsizesepmask_%03d_%03d.fits'%(lvl_sep, size_sep), overwrite = True )

    # Measure Fractions and uncertainties
    #wr = np.sum(np.array(wrl)**2)
    det_err = np.sum(icl_dei**2)
    flux_sample, F_ICL_m, F_ICL_low, F_ICL_up = selection_error(tot_icl_al, icl_boot_al, M = N_err, percent = per_err, xs = xs, ys = ys, mscann = mscann)
    
    icl_flux = np.sum(icl[mscann.astype(bool)])
    sel_err_up = F_ICL_up - F_ICL_m
    sel_err_low = F_ICL_m - F_ICL_low
    tot_err_up = np.sqrt( det_err + sel_err_up**2 )
    tot_err_low = np.sqrt( det_err + sel_err_low**2 )

    print('\nWS + SF + SS -- ICL+BCG -- z = %d    sise_sep = %d'%(lvl_sep, size_sep))
    print('N = %4d   F_ICL = %f ADU  err_low = %f ADU  err_up = %f ADU'%(len(icl_al), F_ICL_m, F_ICL_low, F_ICL_up))
    print(tot_err_up, tot_err_low)
    # Plot vignets
    if plot_vignet == True:

        interval = AsymmetricPercentileInterval(5, 99.5) # meilleur rendu que MinMax or ZScale pour images reconstruites
        fig, ax = plt.subplots(2, 3)
        poim = ax[0][0].imshow(im_art, norm = ImageNormalize( im_art, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')
        poim = ax[1][0].imshow(icl, norm = ImageNormalize( icl, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')
        poim = ax[0][1].imshow(im_unclass, norm = ImageNormalize( recim, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')
        poim = ax[1][1].imshow(recim, norm = ImageNormalize( recim, interval = interval, stretch = LogStretch()), cmap = 'binary', origin = 'lower')
        poim = ax[0][2].imshow(icl_err, norm = ImageNormalize( icl_err, interval = MinMaxInterval(), stretch = LinearStretch()), cmap = 'binary', origin = 'lower')
        poim = ax[1][2].hist(flux_sample, bins = 10)
        
        for i in range(0,2):
            for j in range(0,3):
                if (i != 1) & ( j != 2 ):
                    ax[i][j].get_xaxis().set_ticks([])
                    ax[i][j].get_yaxis().set_ticks([])
                
        #plt.show()
        plt.tight_layout()
        
        plt.savefig( nfap + '.results.bcgwavsizesepmask_%03d_%03d.png'%(lvl_sep, size_sep), format = 'png' )
        print('Write vignet to' + nfap + 'synth.bcgwavsizesepmask_%03d_%03d.png'%(lvl_sep, size_sep))
        plt.close('all')
    
    return icl_flux, tot_err_up, tot_err_low

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if __name__ == '__main__':

    # Paths, lists & variables
    path_data = '/n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3/vignets'
    path_wavelets = '/n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3/wavelets/out14/'
    path_analysis = '/n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3/wavelets/out14/'

    # spatial scales related
    pix_scale = 0.3 # "/pix
    physcale = 3.68 # kpc/"
    size_sep = 100 # size separation [kpc]
    size_sep_pix = size_sep / physcale / pix_scale # pixels
    R_kpc = 1000 # radius in which quantities are measured [kpc]
    R_pix = R_kpc / physcale / pix_scale # pixels
    R = R_pix
    rc = 10 # kpc, distance to center to be classified as gal
    rc_pix = rc / physcale / pix_scale # pixels
                
    # wavelet scales related
    lvl_sep = 5 # wavelet scale separation
    lvl_sep_bcg = 4
    lvl_sep_max = 1000
    n_levels = 10
    
    # photometry
    mu_lim = 35
    ZP_AB = 30
    flux_lim = 10**( (ZP_AB - mu_lim) / 2.5 )
    
    # bootstrap
    N_err = 100
    per_err = 0.05
    
    # misc
    kurt_filt = True
    plot_vignet = False
    write_fits = True
    measure_PR = True
    write_dataframe = True
    plot_boot = True
    
    mscicl = fits.getdata('/n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3/icl_mock_circle_2arcmin_vignet.fits')
    
    col_ICL_flux = []
    col_tot_err_up = []
    col_tot_err_low = []
    col_re = []
    col_fICL = []
    col_z = []
    col_cl_name = []
    col_num_vignet = []
    
    for nfp in glob.glob(os.path.join(path_data, '*vignet_?.fits' )):
                
        nf = nfp.split('/')[-1]
        split = nf.split('_')
        num_vignet = split[-1][0]

        print(nfp)
          
        hdu = fits.open(nfp)
        head = hdu[0].header
        oim = hdu[0].data
        
        xs, ys = oim.shape
        xc, yc = xs / 2., ys / 2.

        mscbcg = np.zeros((xs, ys))
        mscann = np.ones((xs, ys))
        
        nf = nf[:-5]
        
        nfwp = os.path.join(path_wavelets, nf)
        
        nfap = path_data

        ficl, tot_err_up, tot_err_low = synthesis_bcgwavsizesep_with_masks( nfwp, nfap, lvl_sep, lvl_sep_max, lvl_sep_bcg,
                                           size_sep, size_sep_pix, xs, ys,
                                           n_levels, mscicl, mscbcg, mscann,
                                           N_err = N_err,
                                           per_err = per_err,
                                           kurt_filt = kurt_filt,
                                           plot_vignet = plot_vignet,
                                           write_fits = write_fits )

        col_ICL_flux.append(ficl)
        col_tot_err_up.append(tot_err_up)
        col_tot_err_low.append(tot_err_low)
        col_num_vignet.append(num_vignet)
        col_cl_name.append(nf)
    
    df = pd.DataFrame(columns = ['cl_name', 'num_vignet', 'ICL_flux', 'ICL_flux_err_hi', 'ICL_flux_err_low'])
    df['cl_name'] = col_cl_name
    df['num_vignet'] = col_num_vignet
    df['ICL_flux'] = col_ICL_flux
    df['ICL_flux_err_hi'] = col_tot_err_up
    df['ICL_flux_err_low'] = col_tot_err_low
 
    print('Write results to %s' %os.path.join(path_data, 'Euclid-NISP-H-ERO-Abell2390-LSB.DR3.mocks.csv'))
    df.to_csv(os.path.join(path_data, 'Euclid-NISP-H-ERO-Abell2390-LSB.DR3.mocks.csv'))
