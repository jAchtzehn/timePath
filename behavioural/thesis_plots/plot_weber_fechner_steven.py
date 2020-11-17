'''
Plot webers ratio and fechenr law
'''


from os.path import join as opj
from os.path import abspath
import numpy as np
import os
import pandas as pd
import psignifit as ps
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pylab as pl
import math

experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
behavioral_dir = opj(experiment_dir, 'behavioural')
output_dir = abspath('/Users/jachtzehn/Documents/Medizin/thesis/figures/basics/')

# plotting options
matplotlib.rc('text', usetex=True)
pl.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',    # helvetica font
    r'\usepackage{sansmath}',   # math-font matching  helvetica
    r'\sansmath'                # actually tell tex to use it!
    r'\usepackage{siunitx}',    # micro symbols
    r'\sisetup{detect-all}',    # force siunitx to use the fonts
]
matplotlib.rcParams['ytick.labelsize'] = 26
matplotlib.rcParams['xtick.labelsize'] = 26
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['axes.labelsize'] = 32
matplotlib.rcParams['axes.labelpad'] = 12
figsize = (16, 7)
tight_pad = 5
lw = 5

ln_scale = 0.43
pw_scale = 1

def linearFunc(x):
	y = x / 10
	return y

def logFunc(x):
	y = ln_scale * np.log(x)
	return y

def powerFunc(x, b):
	y = pw_scale * x ** b
	return y

fig, ax = plt.subplots(1, 2, figsize=figsize)

x_lin = np.linspace(0, 10.3, num=1000)
x_log = np.linspace(1, 10.3, num=1000)

ax[0].plot(x_lin, linearFunc(x_lin), lw=lw, zorder=3, color='black')
ax[1].plot(x_log, logFunc(x_log), lw=lw, zorder=3, color='black')

ax[0].yaxis.set_major_locator(MaxNLocator(prune='lower'))
ax[1].yaxis.set_major_locator(MaxNLocator(prune='lower'))

ax[0].set_yticks([0.2, 0.4, 0.6, 0.8, 1])
ax[0].set_xticks([2, 4, 6, 8, 10])
ax[1].set_yticks([0.2, 0.4, 0.6, 0.8, 1])
ax[1].set_xticks([2, 4, 6, 8, 10])
ax[0].set_ylim([0, 1.02])
ax[0].set_xlim([0, 10.5])
ax[1].set_ylim([0, 1.02])
ax[1].set_xlim([0, 10.5])

ax[0].set_xlabel('Stimulus intensity I (a.u.)')
ax[1].set_xlabel('Stimulus intensity I (a.u.)')
ax[0].set_ylabel(r'MOD $\Delta$I (a.u.)')
ax[1].set_ylabel('Perceived sensation E (a.u.)')

# draw horizontal lines
for y_loc in np.linspace(0.2, 1, num=5):
	x_len = np.exp(y_loc * 1/ln_scale)
	ax[1].plot([0, x_len], [y_loc, y_loc], '--', zorder=2, lw=lw/3, color='Dimgray')

# draw horizontal lines
x_values = [np.exp(y_loc * 1/ln_scale) for y_loc in np.linspace(0.2, 1., num=5)]
for x_loc in x_values:
	ax[1].plot([x_loc, x_loc], [0, logFunc(x_loc)], '--', zorder=2, lw=lw/3, color='Dimgray')

ax[0].spines['top'].set_visible(False)
ax[1].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)
ax[1].spines['right'].set_visible(False)

plt.tight_layout(w_pad=5)
plt.savefig(opj(output_dir, 'weber-fechner.pdf'), format='pdf', dpi=300, fig=fig)
plt.close()


# stevens
fig, ax = plt.subplots(1, 1, figsize=(13, 8))

x = np.linspace(0, 10, num=1000)
ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))


ax.plot(x, powerFunc(x, 1), lw=lw, label=r'Apparent length ($b \approx 1$)', color='slategray')
ax.plot(x, powerFunc(x, 2.4), lw=lw, label=r'Electric shock ($b > 1$)', color='IndianRed')
ax.plot(x, powerFunc(x, 0.3), lw=lw, label=r'Brightness ($b < 1$)', color='darkseagreen')
ax.plot(x, powerFunc(x, 0.6), lw=lw, label=r'Loudness ($b < 1$)', color='cornflowerblue')

ax.set_ylim([0, 4])
ax.set_xlim([0, 5])

ax.set_xlabel('Stimulus intensity I (a.u.)')
ax.set_ylabel('Perceived sensation E (a.u.)')
plt.legend(loc='lower right', ncol=1, fontsize=24)

plt.tight_layout()
plt.savefig(opj(output_dir, 'steven.pdf'), format='pdf', dpi=300)