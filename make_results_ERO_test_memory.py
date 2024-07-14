import sys
import dawis as d
import glob as glob
import os
import numpy as np
import pyregion as pyr
import random
import gc
import h5py
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import *
from scipy.stats import kurtosis
from power_ratio import *
from datetime import datetime
from photutils.segmentation import SourceCatalog, detect_sources
import make_results_ERO_noray
import tracemalloc

# Paths, lists & variables
path_data = '/home/aellien/Euclid_ERO/data/Euclid-NISP-Stack-ERO-Abell2390.DR3'
path_scripts = '/home/aellien/Euclid_ERO/Euclid_ERO_scripts'
path_wavelets = '/home/aellien/Euclid_ERO/wavelets/out7/'
path_plots = '/home/aellien/Euclid_ERO/plots'
path_analysis = '/home/aellien/Euclid_ERO/analysis/'

# Input files
nfl = [ 'Euclid-NISP-J-ERO-Abell2390-LSB.DR3.crop.fits' ]#, 
       #'Euclid-NISP-H-ERO-Abell2390-LSB.DR3.crop.fits', 
       #'Euclid-NISP-Y-ERO-Abell2390-LSB.DR3.crop.fits' ]

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
lvl_sep_bcg = 6
lvl_sep_max = 9
n_levels = 10

# photometry
mu_lim = 35
ZP_AB = 30
flux_lim = 10**( (ZP_AB - mu_lim) / 2.5 )

# bootstrap
N_err = 50
per_err = 0.1

# misc
kurt_filt = True
plot_vignet = False
write_fits = True
measure_PR = True
write_dataframe = True
plot_boot = True

# Masks
mscell = fits.getdata(os.path.join(path_analysis,'mscell.fits'))
mscbcg = fits.getdata(os.path.join(path_analysis,'mscbcg.fits'))
mscstar = fits.getdata(os.path.join(path_analysis,'mscstar.fits'))
msat = fits.getdata(os.path.join(path_analysis,'msat.fits'))
mscsedl = []

# Iterate over dictionary list
for nf in nfl:
    
    # Read image data from dic
    filt = nf.split('-')[2]
    print(nf)
    
    #Read image file
    nfp = os.path.join( path_wavelets, nf[:-4] )
    oim_file = os.path.join( path_data, nf )
    hdu = fits.open(oim_file)
    oim = hdu[0].data
    xs, ys = oim.shape

    # Synthesis
    #output_df = make_results_ERO_noray.synthesis_bcgwavsizesep_with_masks( nfp = nfp, 
    #                                         lvl_sep = lvl_sep, 
    #                                         lvl_sep_max = lvl_sep_max, 
    #                                         lvl_sep_bcg = lvl_sep_bcg, 
    #                                         size_sep = size_sep, 
    #                                         size_sep_pix =  size_sep_pix, 
    #                                         xs = xs, 
    #                                         ys = ys, 
    #                                         n_levels = n_levels,
    #                                         mscstar = mscstar, 
    #                                         mscell = mscell, 
    #                                         mscbcg = mscbcg, 
    #                                         mscsedl = mscsedl,
    #                                         msat = msat, 
    #                                         R = R_pix, 
    #                                         rc_pix = rc_pix,
    #                                         N_err = N_err, 
    #                                         per_err = per_err, 
    #                                         flux_lim = flux_lim, 
    #                                         kurt_filt = kurt_filt, 
    #                                         plot_vignet = plot_vignet, 
    #                                         write_fits = write_fits, 
    #                                         measure_PR = measure_PR,
    #                                         plot_boot = plot_boot)
    #if write_dataframe:
    #    ofp = nfp + 'df.csv'
    #    print('Write results to %s'%ofp)
    #    output_df.to_csv(ofp)

    tracemalloc.start()
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

    memory = []
    marea = []
    for i, ( op, itlp ) in enumerate( zip( opathl, itpathl )):
        
        print('Iteration %d' %(i))#, end ='\r')
        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics('lineno')
        print(top_stats[0])
        memory.append(top_stats[0].size/1e9)
        
        
        #ol = d.store_objects.read_ol_from_hdf5(op)
        #itl = d.store_objects.read_itl_from_hdf5(itlp)
        with h5py.File(op, "r") as f1, h5py.File(itlp, "r") as f2:
            gc.collect()
            icl_al = []
            gal_al = []
            noticl_al = []
            unclass_al = []
            areal = []
            for o, it in zip(f1.keys(), f2.keys()):

                x_min, y_min, x_max, y_max = np.copy(f1[o]['bbox'][()])
                image = np.copy(f1[o]['image'][()])
                det_err_image = np.copy(f1[o]['det_err_image'][()])
                itm = np.copy(f2[it]['interscale_maximum'])
                lvlo = np.copy(f1[o]['level'][()])
                #al.append([image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo])
                ######################################## MEMORY ^
            
                #for (image, det_err_image, x_min, y_min, x_max, y_max, xco, yco, lvlo) in al
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
                        coo_spur_halo =  [ [2142, 2216], [1890, 2270] ] # pix long, ds9 convention
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
            marea.append(np.mean(areal))
    # clear some memory
    gc.collect()
    icl_al.clear()
    gal_al.clear()
    noticl_al.clear()
    unclass_al.clear()    
    
    '''
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
        hdul.writeto( nfp + 'synth.bcgwavsizesepmask_%03d_%03d.fits'%(lvl_sep, size_sep), overwrite = True )'''

    '''if measure_PR == True:
        print('start bootstrap')
        start = datetime.now()
        # Measure Fractions and uncertainties
        F_ICL_m, F_ICL_low, F_ICL_up, FICL_det_err, low_FICL_det_err, up_FICL_det_err, out_sed_icl =  make_results_ERO_noray.selection_error(tot_icl_al, tot_unclass_al+tot_gal_al, M = N_err, percent = per_err, xs = xs, ys = ys, flux_lim = flux_lim, mscsedl = mscsedl, write_plots = plot_boot, nfp = nfp)
        F_gal_m, F_gal_low, F_gal_up, Fgal_det_err, low_Fgal_det_err, up_Fgal_det_err, out_sed_gal =  make_results_ERO_noray.selection_error(tot_gal_al, tot_unclass_al+tot_icl_al, M = N_err, percent = per_err, xs = xs, ys = ys, flux_lim = flux_lim, mscsedl = mscsedl)
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
        results_PR = make_results_ERO_noray.PR_with_selection_error(atom_in_list = tot_icl_al, atom_out_list = tot_unclass_al+tot_gal_al, M = N_err, percent = per_err, R = R, xs = xs, ys = ys)
        PR_1_m, PR_1_up, PR_1_low = results_PR[0]
        PR_2_m, PR_2_up, PR_2_low = results_PR[1]
        PR_3_m, PR_3_up, PR_3_low = results_PR[2]
        PR_4_m, PR_4_up, PR_4_low = results_PR[3]
        
        print('PR_1_m = %1.2e    PR_1_low = %1.2e    PR_1_up = %1.2e'%(PR_1_m, PR_1_low, PR_1_up))
        print('PR_2_m = %1.2e    PR_2_low = %1.2e    PR_2_up = %1.2e'%(PR_2_m, PR_2_low, PR_2_up))
        print('PR_3_m = %1.2e    PR_3_low = %1.2e    PR_3_up = %1.2e'%(PR_3_m, PR_3_low, PR_3_up))
        print('PR_4_m = %1.2e    PR_4_low = %1.2e    PR_4_up = %1.2e'%(PR_4_m, PR_4_low, PR_4_up))'''

    snapshot = tracemalloc.take_snapshot()
    top_stats = snapshot.statistics('lineno')
    for stat in top_stats[:20]:
        print(stat)