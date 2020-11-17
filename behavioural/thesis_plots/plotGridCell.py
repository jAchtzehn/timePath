from os.path import join as opj
from os.path import abspath
import pandas as pd
import numpy as np
import matplotlib
from scipy.io import loadmat
import matplotlib.pyplot as plt
import pylab as pl

# plotting options
matplotlib.rc('text', usetex=True)
pl.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',    # helvetica font
    r'\usepackage{sansmath}',   # math-font matching  helvetica
    r'\sansmath'                # actually tell tex to use it!
    r'\usepackage{siunitx}',    # micro symbols
    r'\sisetup{detect-all}',    # force siunitx to use the fonts
]
matplotlib.rcParams['ytick.labelsize'] = 22
matplotlib.rcParams['xtick.labelsize'] = 22
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['axes.labelsize'] = 28
matplotlib.rcParams['axes.labelpad'] = 12
matplotlib.rcParams['axes.titlesize'] = 32
matplotlib.rcParams['axes.titlepad'] = 20

figsize = (10, 10)
tight_pad = 5
lw = 3
ms = 10

output_dir = abspath('/Users/jachtzehn/Documents/Medizin/thesis/figures/basics/')
source_dir = abspath('/Users/jachtzehn/Documents/Development/privat/Grid_Cell_simulation-master_matlab')

data = loadmat(opj(source_dir, 'simulation.mat'))
autocorr = data['autocorr']
recField = data['rec_field']

#autocorr[autocorr == 0] = np.nan
#recField[recField == 0] = np.nan

fig1, ax1 = plt.subplots(1, 1, figsize=figsize)
im_recField = ax1.imshow(recField, cmap='viridis', interpolation='bicubic', extent=[0, 64, 64, 0])
fig1.colorbar(im_recField, ax=ax1, fraction=0.04675, pad=0.02, orientation='horizontal')

fig2, ax2 = plt.subplots(1, 1, figsize=figsize)
im_autocorr = ax2.imshow(autocorr, cmap='viridis', interpolation='bicubic', extent=[0, 127, 127, 0])
fig2.colorbar(im_autocorr, ax=ax2, fraction=0.04675, pad=0.02, orientation='horizontal')

ax1.set_axis_off()
ax2.set_axis_off()
fig1.tight_layout()

fig2.tight_layout()

fig1.savefig(opj(output_dir, 'recFields.pdf'), format='pdf', dpi=300)
fig2.savefig(opj(output_dir, 'autoCorr.pdf'), format='pdf', dpi=300)

#plt.show()
