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
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import *
from scipy.stats import kurtosis
from cosmo_calc import cosmo_calc

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

        for (o, xco, yco) in draw:
            x_min, y_min, x_max, y_max = o.bbox
            
            im_s[ x_min : x_max, y_min : y_max ] += o.image
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
def synthesis_bcgwavsizesep_with_masks( nfwp, nfap, lvl_sep, lvl_sep_max, lvl_sep_bcg, size_sep, size_sep_pix, xs, ys, n_levels, mscicl, mscbcg, mscann, N_err, per_err, kurt_filt = True, plot_vignet = False, write_fits = True, measure_flux = False ):
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
    icl_err = np.zeros( (xs, ys) )
    im_art = np.zeros( (xs, ys) )
    recim = np.zeros( (xs, ys) )
    im_unclass = np.zeros( (xs, ys) )

    icl_al = []
    gal_al = []
    noticl_al = []
    unclass_al = []
    
    wrl = [] # reconstruction error list

    # Read atoms
    ol, itl = d.read_image_atoms( nfwp, file_format = 'hdf5', verbose = True )

    # Kurtosis + ICL+BCG
    for j, o in enumerate(ol):

        x_min, y_min, x_max, y_max = o.bbox
        sx = x_max - x_min
        sy = y_max - y_min
        itm = itl[j].interscale_maximum
        xco = itm.x_max
        yco = itm.y_max
        
        recim[ x_min : x_max, y_min : y_max ] += o.image

        if kurt_filt == True:
            k = kurtosis(o.image.flatten(), fisher=True)
            if k < 0:
                im_art[ x_min : x_max, y_min : y_max ] += o.image
                continue

        # Remove background
        if o.level >= lvl_sep_max:
            continue

        # ICL + BCG
        if mscicl[xco, yco] == 1:

            # BCG
            if mscbcg[xco, yco] == 1:
                icl[ x_min : x_max, y_min : y_max ] += o.image
                icl_err[ x_min : x_max, y_min : y_max ] += o.det_err_image
                icl_al.append([o, xco, yco])
                wrl.append(o.norm_wr)
                
            # ICL
            else:

                if (o.level >= lvl_sep) & (sx >= size_sep_pix) & (sy >= size_sep_pix):

                    icl[ x_min : x_max, y_min : y_max ] += o.image
                    icl_err[ x_min : x_max, y_min : y_max ] += o.det_err_image
                    icl_al.append([o, xco, yco])
                    wrl.append(o.norm_wr)
                    
                else:
                    im_unclass[ x_min : x_max, y_min : y_max ] += o.image
                    if mscann[xco, yco] == 1:
                        noticl_al.append([o, xco, yco])

        else:
            im_unclass[ x_min : x_max, y_min : y_max ] += o.image
            if mscann[xco, yco] == 1:
                noticl_al.append([ o, xco, yco ])

    if write_fits == True:
        print('\nWS + SF + SS -- ICL+BCG -- write fits as %s*'%(nfap))

        # write to fits
        hduo = fits.PrimaryHDU(icl)
        hduo.writeto( nfap + '.synth.icl.bcgwavsizesepmask_%03d_%03d.fits'%(lvl_sep, size_sep), overwrite = True )

    # Measure Fractions and uncertainties
    wr = np.sum(np.array(wrl)**2)
    det_err = np.sum(icl_err**2)
    flux_sample, F_ICL_m, F_ICL_low, F_ICL_up = selection_error(icl_al, noticl_al, M = N_err, percent = per_err, xs = xs, ys = ys, mscann = mscann)
    
    icl_flux = np.sum(icl[mscann.astype(bool)])
    sel_err_up = F_ICL_up - F_ICL_m
    sel_err_low = F_ICL_m - F_ICL_low
    tot_err_up = np.sqrt( wr + det_err + sel_err_up**2 )
    tot_err_low = np.sqrt( wr + det_err + sel_err_low**2 )

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
    path_analysis = '/n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3/'

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
    #mscbcg = fits.getdata('/n03data/ellien/Euclid_ICL/simulations/out2/mscbcg.fits')
    
    col_ICL_flux = []
    col_tot_err_up = []
    col_tot_err_low = []
    col_re = []
    col_fICL = []
    col_z = []
    col_cl_name = []
    col_num_vignet = []
    
    for nfp in glob.glob(os.path.join(path_data, '*_1.fits' )):
                
        nf = nfp.split('/')[-1]
        split = nf.split('_')
        num_vignet = split[-1][0]
          
        hdu = fits.open(nfp)
        head = hdu[0].header
        oim = hdu[0].data
        
        xs, ys = oim.shape
        xc, yc = xs / 2., ys / 2.

        mscbcg = np.zeros((xs, ys))
        
        nf = nfp.split('/')[-1][:-5]
        
        #/n03data/ellien/Euclid_ICL/wavelets/out2/360009933000018/360009933000018_EUC_NIR_W-STK-IMAGE_H_z_0.3_fICL0.15_re_1.0_galsim_swarp_grid_bgsub_vignet_2/
        nfwp = os.path.join(path_wavelets, nf)
        
        # /n03data/ellien/Euclid_ICL/analysis/out2/373000139000019
        nfap = path_data

        ficl, tot_err_up, tot_err_low = synthesis_bcgwavsizesep_with_masks( nfwp, nfap, lvl_sep, lvl_sep_max, lvl_sep_bcg,
                                           size_sep, size_sep_pix, xs, ys,
                                           n_lvl, mscicl, mscbcg, mscann,
                                           N_err = N_err,
                                           per_err = per_err,
                                           kurt_filt = kurt_filt,
                                           plot_vignet = plot_vignet,
                                           write_fits = write_fits,
                                           measure_flux = measure_flux )
                
        col_ICL_flux.append(ficl)
        col_tot_err_up.append(tot_err_up)
        col_tot_err_low.append(tot_err_low)
        col_cl_name.append(cluster_name)
        col_re.append(re)
        col_fICL.append(fICL)
        col_z.append(z)
        col_num_vignet.append(num_vignet)
    
    df = pd.DataFrame(columns = ['cl_name', 'num_vignet', 'ICL_flux', 'ICL_flux_err_hi', 'ICL_flux_err_low'])
    df['cl_name'] = nf
    df['num_vignet'] = col_num_vignet
    df['ICL_flux'] = col_ICL_flux
    df['ICL_flux_err_hi'] = col_tot_err_up
    df['ICL_flux_err_low'] = col_tot_err_low
 
    print('Write results to %s' %os.path.join(path_data, 'Euclid-NISP-H-ERO-Abell2390-LSB.DR3.mocks.csv'))
    df.to_csv(os.path.join(path_data, 'Euclid-NISP-H-ERO-Abell2390-LSB.DR3.mocks.csv'))
