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
#infile = "Euclid-NISP-J-ERO-Abell2390-LSB.DR3.crop.fits"
outdir = '/n03data/ellien/Euclid_ERO/Euclid-NISP-Stack-ERO-Abell2390.DR3/wavelets/out8/'
#outdir = '/home/aellien/Euclid_ERO/wavelets/local/run10'

n_cpus = 8 # Number of CPUs
tau = 0.8   # Relative Threshold
gamma = 0.2   # Attenuation (CLEAN) factor
ceps = 1E-4    # Convergence value for epsilon
n_levels = 9    # Number of wavelet scales
min_span = 1    # Minimum of wavelet scales spanned by an interscale tree (must be >= 1)
max_span = 2    # Maximum number of wavelet scales spanned by an interscale tree
lvl_sep_big = 8     # Scale at wich mix_span, max_span & gamma are set to 1, and monomodality is enforced
extent_sep = 0.1   # Ratio n_pix/vignet under which the Haar wavelet is used for restoration
ecc_sep = 0.05      # Eccentricity threshold over which the Haar wavelet is used for restoration
lvl_sep_lin = -1     # Wavelet scale under which the Haar wavelet can be used for restoration
max_iter = 1500      # Maximum number of iterations
data_dump = True    # Write data at each iteration /!\ demands lot of space on hardware /!\
gif = True      # Make gifs of the run (need data_dump = True)
starting_level = 2 # Starting wavelet scale (this is the third scale - Python convention 0 1 2)
conditions = 'prolongation'
monomodality = True
threshold_rel = 0.15
size_patch = 100
resume = True
rm_gamma_for_big = True

shutil.copyfile( os.path.abspath(__file__), os.path.join( outdir, infile[:-4] + 'input.dawis.py' ) )

dawis.synthesis_by_analysis( indir = indir, infile = infile, outdir = outdir, n_cpus = n_cpus, n_levels = n_levels, \
                                    tau = tau, gamma = gamma, ceps = ceps, conditions = conditions, min_span = min_span, \
                                    max_span = max_span, lvl_sep_big = lvl_sep_big, monomodality = monomodality, threshold_rel = threshold_rel, \
                                    max_iter = max_iter, extent_sep = extent_sep, ecc_sep = ecc_sep, size_patch = size_patch, rm_gamma_for_big = rm_gamma_for_big, \
                                    data_dump = data_dump, gif = gif, resume = resume )

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