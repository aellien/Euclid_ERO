import os
import dawis
import sys
import shutil
import glob

indir = '/home/aellien/Euclid_ERO/data/Euclid-VIS-Stack-ERO-Abell2390.DR3/'
infile = 'cirrus_first_pass.fits'
outdir = '/home/aellien/Euclid_ERO/wavelets/out8/'

tau = 0.1   # Relative Threshold
gamma = 1.0   # Attenuation (CLEAN) factor

ceps = 1E-4    # Convergence value for epsilon
scale_lvl_eps = 1 # Scale convergence value with wavelet scale
max_iter = 1500      # Maximum number of iterations

starting_level = 6 # Starting wavelet scale (this is the third scale - Python convention 0 1 2)
n_levels = 11    # Number of wavelet scales
min_span = 1    # Minimum of wavelet scales spanned by an interscale tree (must be >= 1)
max_span = 2    # Maximum number of wavelet scales spanned by an interscale tree
deblend_contrast = 0.01 # contrast parameter for deblending
lvl_deblend = 3 # Scale at which the regions of significant wavelet coefficients are deblended
lvl_sep_big = 5     # Scale at wich mix_span, max_span & gamma are set to 1, and monomodality is enforced
lvl_sep_op = 8  # Scale at which synthesis operator switch from SUM to ADJOINT
rm_gamma_for_big = True # If set to true, the attenuation factor is not applied for scales higher than lvl_sep_big

extent_sep = 0.15   # Ratio n_pix/vignet under which the Haar wavelet is used for restoration
ecc_sep = 0.9       # Eccentricity threshold over which the Haar wavelet is used for restoration
lvl_sep_lin = -1    # Wavelet scale under which the Haar wavelet can be used for restoration

n_sigmas = 3 # Threshold for detection
inpaint_res = True # If set to true high negative value residual will be inpainted by noise
iptd_sigma = 5  # Threshold for noise inpainted values

data_dump = False   # Write data at each iteration /!\ demands lot of space on hardware /!\
gif = True      # Make gifs of the run (need data_dump = True)
conditions = 'prolongation' # Border conditions for wavelet convolution

outfile_format = 'hdf5'
resume = True 
deconv = True

n_cpus = 8

if os.path.isdir( outdir ) == False:
    os.makedirs( outdir, exist_ok = True )

shutil.copyfile( os.path.abspath(__file__), os.path.join( outdir, 'input.dawis.py' ) )

dawis.synthesis_by_analysis( indir = indir, infile = infile, outdir = outdir, n_cpus = n_cpus, starting_level = starting_level, tau = tau, n_levels = n_levels, n_sigmas = n_sigmas, outfile_format = outfile_format,\
                                gamma = gamma, min_span = min_span, max_span = max_span, deblend_contrast = deblend_contrast, lvl_sep_big = lvl_sep_big, rm_gamma_for_big = rm_gamma_for_big, lvl_deblend = lvl_deblend, \
                                extent_sep = extent_sep, ecc_sep = ecc_sep, lvl_sep_op = lvl_sep_op, lvl_sep_lin = lvl_sep_lin, ceps = ceps, scale_lvl_eps = scale_lvl_eps, conditions = conditions, deconv = deconv, \
                                max_iter = max_iter, inpaint_res = inpaint_res, data_dump = data_dump, gif = gif, iptd_sigma = iptd_sigma, resume = resume )
