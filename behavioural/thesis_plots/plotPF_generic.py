'''
Plot generic pf
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
figsize = (14, 10)
tight_pad = 5
lw = 5


data = np.array([[0.10,   0.0000,   90.0000],
                 [0.15,   10.0000,   90.0000],
                 [0.20,   10.0000,   90.0000],
                 [0.25,   10.0000,   90.0000],
                 [0.30,   52.0000,   90.0000],
                 [0.35,   53.0000,   90.0000],
                 [0.40,   62.0000,   90.0000],
                 [0.45,   64.0000,   90.0000],
                 [0.50,   76.0000,   90.0000],
                 [0.60,   79.0000,   90.0000],
                 [0.70,   88.0000,   90.0000],
                 [0.80,   90.0000,   90.0000],
                 [0.90,   90.0000,   90.0000]])


fig, ax = plt.subplots(1, 1, figsize=figsize)
plt.tight_layout(tight_pad)

plotDiffFunctions = False

options = dict()   # initialize as an empty dict


### norm
options['sigmoidName'] = 'norm'
options['expType'] = 'equalAsymptote'

result = ps.psignifit(data, options)

fit = result['Fit']
data = result['data']
options_res = result['options']

xData = data[:, 0]
yData = data[:, 1] / data[:, 2]
xMin = min(xData)
xMax = max(xData)
xLength = xMax - xMin
x = np.linspace(xMin, xMax, num=1000)
xLow = np.linspace(xMin - .2*xLength, xMin, num=100)
xLow = xLow[xLow >= 0]
xHigh = np.linspace(xMax, xMax + .2*xLength, num=100)
xHigh = xHigh[xHigh <= 1]
fitValues = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](x, fit[0], fit[1]) + fit[3]
fitValuesLow = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](xLow, fit[0], fit[1]) + fit[3]
fitValuesHigh = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](xHigh, fit[0], fit[1]) + fit[3]

if plotDiffFunctions:
    ax.plot(x, fitValues, color='black', linewidth=lw, label='PF: CND')
else:
    ax.plot(x, fitValues, color='black', linewidth=lw, label='Fitted PF')

plt.plot(xLow, fitValuesLow, '--', c='black', linewidth=lw)
plt.plot(xHigh, fitValuesHigh, '--', c='black', linewidth=lw)

# threshhold x
x = [fit[0], fit[0]]
y = [0, fit[3] + (1 - fit[2] - fit[3]) * .5]
plt.plot(x, y, '-', c='DimGray', linestyle='--')

# threshhold y
x = [0, fit[0]]
y = [.5, .5]
plt.plot(x, y, '-', c='DimGray', linestyle='--')

#### logistic
options['sigmoidName'] = 'logistic'
result = ps.psignifit(data, options)

fit = result['Fit']
data = result['data']
options_res = result['options']

xData = data[:, 0]
yData = data[:, 1] / data[:, 2]
xMin = min(xData)
xMax = max(xData)
xLength = xMax - xMin
x = np.linspace(xMin, xMax, num=1000)
xLow = np.linspace(xMin - .2*xLength, xMin, num=100)
xLow = xLow[xLow >= 0]
xHigh = np.linspace(xMax, xMax + .2*xLength, num=100)
xHigh = xHigh[xHigh <= 1]
fitValues = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](x, fit[0], fit[1]) + fit[3]
fitValuesLow = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](xLow, fit[0], fit[1]) + fit[3]
fitValuesHigh = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](xHigh, fit[0], fit[1]) + fit[3]

if plotDiffFunctions:
    ax.plot(x, fitValues, color='SkyBlue', clip_on=False, linewidth=lw, label='PF: Logistic')

    ax.plot(xLow, fitValuesLow, '--', c='SkyBlue', linewidth=lw)
    ax.plot(xHigh, fitValuesHigh, '--', c='SkyBlue', linewidth=lw)

    # threshhold x
    x = [fit[0], fit[0]]
    y = [0, .5]
    ax.plot(x, y, '-', c='SkyBlue', linestyle='--')


#### logistic
options['sigmoidName'] = 'Weibull'
result = ps.psignifit(data, options)

fit = result['Fit']
data = result['data']
options_res = result['options']

xData = data[:, 0]
yData = data[:, 1] / data[:, 2]
xMin = min(xData)
xMax = max(xData)
xLength = xMax - xMin
x = np.linspace(xMin, xMax, num=1000)
xLow = np.linspace(xMin - .2*xLength, xMin, num=100)
xLow = xLow[xLow >= 0]
xHigh = np.linspace(xMax, xMax + .2*xLength, num=100)
xHigh = xHigh[xHigh <= 1]
fitValues = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](x, fit[0], fit[1]) + fit[3]
fitValuesLow = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](xLow, fit[0], fit[1]) + fit[3]
fitValuesHigh = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](xHigh, fit[0], fit[1]) + fit[3]

if plotDiffFunctions:
    ax.plot(x, fitValues, color='Darkseagreen', linewidth=lw, label='PF: Weibull')

    ax.plot(xLow, fitValuesLow, '--', c='Darkseagreen', linewidth=lw)
    ax.plot(xHigh, fitValuesHigh, '--', c='Darkseagreen', linewidth=lw)

    # threshhold x
    x = [np.exp(fit[0]), np.exp(fit[0])]
    y = [0, .5]
    ax.plot(x, y, '-', c='Darkseagreen', linestyle='--')

# data
ax.scatter(data[:, 0], data[:, 1] / data[:, 2], marker='x', s=200, color='IndianRed', linewidth=lw, label='Input data')
# fit


ax.set_ylim([0, 1.05])
ax.set_xlim([0, 1])
ax.set_xlabel('Stimulus intensity difference')
ax.set_ylabel('Proportion of correct responses')

ax.tick_params(axis="x", direction="out", width=1, length=5)
ax.tick_params(axis="y", direction="out", width=1, length=5)
plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.legend(loc='lower right', ncol=1, fontsize=24)


if plotDiffFunctions:
    plt.savefig(opj(output_dir, 'pf-diff-func.pdf'), format='pdf', dpi=300)
else:
    plt.savefig(opj(output_dir, 'pf.pdf'), format='pdf', dpi=300)
