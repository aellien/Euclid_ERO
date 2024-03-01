import os
import sys
sys.path.append("/home/aellien/dawis")
import dawis
import shutil
import cProfile
from datetime import datetime

indir = '/n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3'
#indir = '/home/aellien/Euclid_ERO/data/Euclid-NISP-Stack-ERO-Abell2390.DR3'
infile = sys.argv[1]
#infile = 'gauss_star_test.fits'#"Euclid-NISP-J-ERO-Abell2390-LSB.DR3.crop.fits"
outdir = '/n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3/wavelets/out14/'
#outdir = '/home/aellien/Euclid_ERO/wavelets/local/run14'
if os.path.isdir( outdir ) == False:
    os.makedirs( outdir, exist_ok = True )
    
tau = 0.8   # Relative Threshold
gamma = 0.2   # Attenuation (CLEAN) factor

ceps = 1E-4    # Convergence value for epsilon
scale_lvl_eps = 1 # Scale convergence value with wavelet scale
max_iter = 1500      # Maximum number of iterations

starting_level = 2 # Starting wavelet scale (this is the third scale - Python convention 0 1 2)
n_levels = 10    # Number of wavelet scales
min_span = 1    # Minimum of wavelet scales spanned by an interscale tree (must be >= 1)
max_span = 2    # Maximum number of wavelet scales spanned by an interscale tree
lvl_deblend = 3 # Scale at which the regions of significant wavelet coefficients are deblended
lvl_sep_big = 5     # Scale at wich mix_span, max_span & gamma are set to 1, and monomodality is enforced
rm_gamma_for_big = False # If set to true, the attenuation factor is not applied for scales higher than lvl_sep_big

extent_sep = 0.15   # Ratio n_pix/vignet under which the Haar wavelet is used for restoration
ecc_sep = 0.9      # Eccentricity threshold over which the Haar wavelet is used for restoration
lvl_sep_lin = -1     # Wavelet scale under which the Haar wavelet can be used for restoration

n_sigmas = 3 # Threshold for detection
inpaint_res = True # If set to true high negative value residual will be inpainted by noise
iptd_sigma = 5  # Threshold for noise inpainted values

data_dump = True    # Write data at each iteration /!\ demands lot of space on hardware /!\
gif = True      # Make gifs of the run (need data_dump = True)
conditions = 'prolongation' # Border conditions for wavelet convolution

n_cpus = 4 # Number of CPUs
size_patch = 100 # Number of objects in parallelized patch

resume = True 

if os.path.isdir( outdir ) == False:
    os.makedirs( outdir, exist_ok = True )

shutil.copyfile( os.path.abspath(__file__), os.path.join( outdir, infile[:-4] + 'input.dawis.py' ) )

dawis.synthesis_by_analysis( indir = indir, infile = infile, outdir = outdir, n_cpus = n_cpus, starting_level = starting_level, tau = tau, n_levels = n_levels, n_sigmas = n_sigmas,\
                                gamma = gamma, min_span = min_span, max_span = max_span, lvl_sep_big = lvl_sep_big, rm_gamma_for_big = rm_gamma_for_big, lvl_deblend = lvl_deblend, \
                                extent_sep = extent_sep, ecc_sep = ecc_sep, lvl_sep_lin = lvl_sep_lin, ceps = ceps, scale_lvl_eps = scale_lvl_eps, conditions = conditions, \
                                max_iter = max_iter, size_patch = size_patch, inpaint_res = inpaint_res, data_dump = data_dump, gif = gif, iptd_sigma = iptd_sigma, resume = resume )


'''start = datetime.now()
cProfile.run('dawis.synthesis_by_analysis( indir = indir, infile = infile, outdir = outdir, n_cpus = n_cpus, n_levels = n_levels, \
                                    tau = tau, gamma = gamma, ceps = ceps, conditions = conditions, min_span = min_span, \
                                    max_span = max_span, lvl_sep_big = lvl_sep_big, monomodality = monomodality, threshold_rel = 0.05, \
                                    max_iter = max_iter, extent_sep = extent_sep, ecc_sep = ecc_sep, size_patch = size_patch, rm_gamma_for_big = rm_gamma_for_big, \
                                    data_dump = data_dump, gif = gif, resume = resume )', '/home/aellien/Euclid_ERO/analysis/profilerstats.txt')
print(datetime.now()-start)
import pstats
from pstats import SortKey
p = pstats.Stats('/home/aellien/Euclid_ERO/analysis/profilerstats.txt')
p.sort_stats(SortKey.CUMULATIVE).print_stats(20)

start = datetime.now()
cProfile.run('dawis.synthesis_by_analysis( indir = indir, infile = infile, outdir = outdir, n_cpus = 1, n_levels = n_levels, \
                                    tau = tau, gamma = gamma, ceps = ceps, conditions = conditions, min_span = min_span, \
                                    max_span = max_span, lvl_sep_big = lvl_sep_big, monomodality = monomodality, threshold_rel = 0.05, \
                                    max_iter = max_iter, extent_sep = extent_sep, ecc_sep = ecc_sep, rm_gamma_for_big = rm_gamma_for_big, \
                                    data_dump = data_dump, gif = gif, resume = resume )', '/home/aellien/Euclid_ERO/analysis/profilerstats2.txt')
print(datetime.now()-start)
p = pstats.Stats('/home/aellien/Euclid_ERO/analysis/profilerstats2.txt')
p.sort_stats(SortKey.CUMULATIVE).print_stats(20)'''