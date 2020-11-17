from os.path import join as opj
from os.path import abspath
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pylab as pl


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
figsize = (14, 10)
tight_pad = 6
lw = 5


def likelihood(a):

	N = 10
	m = 6

	y = (np.math.factorial(10) / (np.math.factorial(m) * np.math.factorial(N - m))) * pow(a, m) * pow(1-a, N-m)
	return y


fig, ax = plt.subplots(1, 1, figsize=figsize)
plt.tight_layout(tight_pad)

x = np.linspace(0, 1, 1000)

y = likelihood(x)

ax.plot(x, y, color='IndianRed', linewidth=lw)
ax.set_xlim([0, 1])
ax.set_ylim([0, 0.275])
ax.set_xlabel('Parameter a')
ax.set_ylabel('Likelihood L(a$|$y)')

ax.tick_params(axis="x", direction="out", width=1, length=5)
ax.tick_params(axis="y", direction="out", width=1, length=5)
plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig(opj(output_dir, 'likelihood-param.pdf'), format='pdf', dpi=300)





