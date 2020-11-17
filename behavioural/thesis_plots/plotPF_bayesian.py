'''
Plot bayesian priors and 2d results for figure in basics
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
from psignifit import utils as _utils
from scipy.signal import convolve as _convn
from psignifit.marginalize import marginalize
from matplotlib.ticker import NullFormatter
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
matplotlib.rcParams['ytick.labelsize'] = 22
matplotlib.rcParams['xtick.labelsize'] = 22
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['axes.labelsize'] = 28
matplotlib.rcParams['axes.labelpad'] = 12
matplotlib.rcParams['axes.titlesize'] = 32
matplotlib.rcParams['axes.titlepad'] = 20
matplotlib.rcParams['axes.titleweight'] = 900

figsize = (14, 10)
tight_pad = 5
lw = 5

ms = 30
plot_prior = False
plot_2d = True
parameters = [0, 1]

colors = ['slategray', 'darkgoldenrod', 'lightcoral', 'cornflowerblue', 'darkseagreen']

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


options = dict()   # initialize as an empty dict


options['sigmoidName'] = 'norm'
options['expType'] = 'equalAsymptote'
options['fixedPars'] = np.array([float('nan'), float('nan'),
                                 float('nan'), float(0),
                                 float('nan')])
options['threshPC'] = .5

result = ps.psignifit(data, options)

fit = result['Fit']
data = result['data']
options_res = result['options']


### prior #####
if plot_prior:

	matplotlib.rcParams['axes.linewidth'] = 1

	fig, ax = plt.subplots(2, 2, figsize=figsize)
	plt.tight_layout(tight_pad)
	fig.subplots_adjust(wspace=0.2)
	fig.subplots_adjust(hspace=0.5)

	if np.size(result['options']['stimulusRange']) <= 1:
		result['options']['stimulusRange'] = np.array([min(data[:, 0]), max(data[:, 0])])
		stimRangeSet = False
	else:
		stimRangeSet = True

	stimRange = result['options']['stimulusRange']
	r = stimRange[1] - stimRange[0]

	if len(np.unique(data[:, 0])) > 1 and not (stimRangeSet):
		widthmin = min(np.diff(np.sort(np.unique(data[:, 0]))))
	else:
		widthmin = 100 * np.spacing(stimRange[1])

	# We use the same prior as we previously used... e.g. we use the factor by
	# which they differ for the cumulative normal function
	Cfactor = (_utils.my_norminv(.95, 0, 1) - _utils.my_norminv(.05, 0, 1)) / \
	          (_utils.my_norminv(1 - result['options']['widthalpha'], 0, 1) - \
	           _utils.my_norminv(result['options']['widthalpha'], 0, 1))
	widthmax = r

	steps = 10000
	theta = np.empty(5)
	for itheta in range(0, 5):
		if itheta == 0:
			x = np.linspace(stimRange[0] - .5 * r, stimRange[1] + .5 * r, steps)
		elif itheta == 1:
			x = np.linspace(min(result['X1D'][itheta]), max(result['X1D'][1], ), steps)
		elif itheta == 2:
			x = np.linspace(0, .5, steps)
		elif itheta == 3:
			x = np.linspace(0, .5, steps)
		elif itheta == 4:
			x = np.linspace(0, 1, steps)

		y = result['options']['priors'][itheta](x)
		theta[itheta] = np.sum(x * y) / np.sum(y)

	if result['options']['expType'] == 'equalAsymptote':
		theta[3] = theta[2]
	if result['options']['expType'] == 'nAFC':
		theta[3] = 1 / result['options']['expN']

	# get limits for the psychometric function plots
	xLimit = [stimRange[0] - .5 * r, stimRange[1] + .5 * r]

	# thresh
	xthresh = np.linspace(xLimit[0], xLimit[1], steps)
	ythresh = result['options']['priors'][0](xthresh)
	wthresh = _convn(np.diff(xthresh), .5 * np.array([1, 1]))
	cthresh = np.cumsum(ythresh * wthresh)

	ax[0, 0].plot(xthresh, ythresh, lw=lw, c='black')
	ax[0, 0].set_xlim(xLimit)
	ax[0, 0].set_title('Threshold')
	ax[0, 0].set_ylabel('Density')

	ax[1, 0].scatter(data[:, 0], np.zeros(data[:, 0].shape), marker='x', s=ms * 2.5, linewidth=lw * .75, color='IndianRed')

	ax[1, 0].set_ylabel('Percent Correct')
	ax[1, 0].set_xlim(xLimit)
	ax[1, 0].set_ylim([-5, 100])
	x = np.linspace(xLimit[0], xLimit[1], steps)
	for idot in range(0, 5):
		if idot == 0:
			xcurrent = theta[0]
			color = colors[idot]
		elif idot == 1:
			xcurrent = min(xthresh)
			color = colors[idot]
		elif idot == 2:
			tix = cthresh[cthresh >= .25].size
			xcurrent = xthresh[-tix]
			color = colors[idot]
		elif idot == 3:
			tix = cthresh[cthresh >= .75].size
			xcurrent = xthresh[-tix]
			color = colors[idot]
		elif idot == 4:
			xcurrent = max(xthresh)
			color = colors[idot]
		y = 100 * (theta[3] + ((1 - theta[2]) - theta[3]) * result['options']['sigmoidHandle'](x, xcurrent, theta[1]))

		ax[1, 0].plot(x, y, '-', lw=lw, c=color)

		ax[0, 0].plot(xcurrent, result['options']['priors'][0](xcurrent), '.', c=color, ms=ms * .75)


	# width
	xwidth = np.linspace(widthmin, 3/Cfactor*widthmax, steps)
	ywidth = result['options']['priors'][1](xwidth)
	wwidth = _convn(np.diff(xwidth), .5*np.array([1,1]))
	cwidth = np.cumsum(ywidth*wwidth)

	ax[0, 1].plot(xwidth, ywidth, lw=lw, c='black')
	# plt.hold(True)
	ax[0, 1].set_xlim([widthmin, 3 / Cfactor * widthmax])
	ax[0, 1].set_title('Width')
	ax[1, 1].scatter(data[:, 0], np.zeros(data[:, 0].size), marker='x', s=ms * 2.5, linewidth=lw * .75, color='IndianRed')
	# plt.hold(True)
	ax[1, 1].set_xlim(xLimit)
	ax[1, 1].set_ylim([-5, 100])

	x = np.linspace(xLimit[0], xLimit[1], steps)

	for idot in range(0, 5):
		if idot == 0:
			xcurrent = theta[1]
			color = colors[idot]
		elif idot == 1:
			xcurrent = min(xwidth)
			color = colors[idot]
		elif idot == 2:
			wix = cwidth[cwidth >= .25].size
			xcurrent = xwidth[-wix]
			color = colors[idot]
		elif idot == 3:
			wix = cwidth[cwidth >= .75].size
			xcurrent = xwidth[-wix]
			color = colors[idot]
		elif idot == 4:
			xcurrent = max(xwidth)
			color = colors[idot]

		y = 100 * (theta[3] + (1 - theta[2] - theta[3]) * result['options']['sigmoidHandle'](x, theta[0], xcurrent))

		ax[1, 1].plot(x, y, '-', lw=lw, c=color)
		ax[0, 1].plot(xcurrent, result['options']['priors'][1](xcurrent), '.', c=color, ms=ms * .75)

	# ax.set_ylim([0, 1.05])
	# ax.set_xlim([0, 1])
	# ax.set_xlabel('Stimulus intensity difference', fontsize=20)
	# ax.set_ylabel('Proportion of correct responses', fontsize=20)
	#
	# ax.tick_params(axis="x", direction="in", width=2, length=5)
	# ax.tick_params(axis="y", direction="in", width=2, length=5)
	# plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))
	#
	# plt.legend(loc='lower right', ncol=1, fontsize=16)

	for curAx in ax:
		curAx[0].tick_params(axis="x", direction="out", width=1, length=5)
		curAx[0].tick_params(axis="y", direction="out", width=1, length=5)
		curAx[1].tick_params(axis="x", direction="out", width=1, length=5)
		curAx[1].tick_params(axis="y", direction="out", width=1, length=5)

		curAx[0].spines['top'].set_visible(False)
		curAx[0].spines['right'].set_visible(False)
		curAx[1].spines['top'].set_visible(False)
		curAx[1].spines['right'].set_visible(False)

	# xlabel
	fig.text(0.5, 0.03, 'Stimulus intensity', ha='center', va='center', fontsize=28)
	fig.text(0.5, 0.52, 'Stimulus intensity', ha='center', va='center', fontsize=28)

	# legend
	# fig.legend(handles=qs, labels=['0%', '25%', '50%', '75%', '100%'], loc="lower center", ncol=4, borderaxespad=0.1, fontsize=24)

	plt.savefig(opj(output_dir, 'pf-bayes-prior.pdf'), format='pdf', dpi=300)
	plt.close()


### 2D ####
if plot_2d:

	matplotlib.rcParams['axes.linewidth'] = 1

	# definitions for the axes
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.65
	bottom_h = left_h = left + width + 0.02
	spacing = 0.007

	rect_img = [left, bottom, width, height]
	rect_margx = [left, bottom + height + spacing, width, 0.2]
	rect_margy = [left + width + spacing, bottom, 0.2, height]

	fig = plt.figure(1, figsize=(13, 13))

	axImg = plt.axes(rect_img)
	axMargx = plt.axes(rect_margx)
	axMargy = plt.axes(rect_margy)

	par1, label1 = _utils.strToDim(str(parameters[0]))
	par2, label2 = _utils.strToDim(str(parameters[1]))
	assert (par1 != par2), 'par1 and par2 must be different numbers to code for the parameters to plot'
	plt.set_cmap('plasma')

	marg, _, _ = marginalize(result, np.array([par1, par2]))
	if par1 > par2:
		marg = marg.T

	e_orig = [result['X1D'][par2][0], result['X1D'][par2][-1],
	          result['X1D'][par1][0], result['X1D'][par1][-1]]
	e = [0.2, 0.8, 0.3, 0.42]

	xbins = (e_orig[1] - e_orig[0]) / len(marg[0])
	ybins = (e_orig[3] - e_orig[2]) / len(marg[0])

	x_low = int((e[0] - e_orig[0])/xbins)
	x_high = len(marg[0]) - int((e_orig[1] - e[1])/xbins)

	y_low = int((e[2] - e_orig[2])/ybins)
	y_high = len(marg[0]) - int((e_orig[3] - e[3])/ybins)

	cfplot = axImg.contour(marg[y_low:y_high, x_low:x_high], extent=e_orig, aspect='auto', interpolation='bicubic', levels=8, cmap='viridis', antialiased=True, linewidths=2.5)
	axImg.clabel(cfplot, inline=1, fontsize=16, fmt='%2.0f')
	#axImg.contourf(marg[y_low:y_high, x_low:x_high], extent=e_orig, aspect='auto', interpolation='bicubic', antialiased=True, cmap='viridis', levels=8)
	#axImg.imshow(marg[y_low:y_high, x_low:x_high], extent=e_orig, aspect='auto', interpolation='bicubic', cmap='viridis')
	axImg.set_ylabel(label1)
	axImg.set_xlabel(label2)

	axImg.tick_params(direction='out', right=False, top=False)
	#for side in ['top', 'right']:
	#	axImg.spines[side].set_visible(False)

	axMargx.tick_params(direction='out')
	axMargy.tick_params(direction='out')
	# plot marginals
	for par in parameters:
		if par == 1:
			ax = axMargx
			orientation = 'vertical'
		else:
			ax = axMargy
			orientation = 'horizontal'

		x = result['marginalsX'][par]
		marginal = result['marginals'][par]
		CI = np.hstack(result['conf_Intervals'][par].T)
		Fit = result['Fit'][par]

		# CI
		xCI = np.array([CI[0], CI[1], CI[1], CI[0]])
		xCI = np.insert(xCI, 1, x[np.logical_and(x >= CI[0], x <= CI[1])])
		yCI = np.array([np.interp(CI[0], x, marginal), np.interp(CI[1], x, marginal), 0, 0])
		yCI = np.insert(yCI, 1, marginal[np.logical_and(x >= CI[0], x <= CI[1])])

		if orientation == 'vertical':
			pat = plt.Polygon(np.array([xCI, yCI]).T, facecolor='slategray', alpha=0.3)
		else:
			pat = plt.Polygon(np.array([yCI, xCI]).T, facecolor='slategray', alpha=0.3)

		ax.add_patch(pat)

		# prior
		xprior = np.linspace(min(x), max(x), 1000)
		if orientation == 'vertical':
			ax.plot(xprior, result['options']['priors'][par](xprior), '--', c='DimGray', linewidth=1)
		else:
			ax.plot(result['options']['priors'][par](xprior), xprior,  '--', c='DimGray', linewidth=1)

		# posterior
		if orientation == 'vertical':
			ax.plot(x, marginal, lw=2, color='black')
		else:
			ax.plot(marginal, x, lw=2, color='black')

		# point estimate
		if orientation == 'vertical':
			ax.plot([Fit, Fit], [0, np.interp(Fit, x, marginal)], 'slategray')
		else:
			ax.plot([0, np.interp(Fit, x, marginal)], [Fit, Fit], 'slategray')

	#axMargx.set_xlim(axImg.get_xlim())
	axMargx.set_xlim([e[0], e[1]])
	axMargx.set_ylim([0, 5.75])
	axMargx.set_yticks([0, 1, 2, 3, 4, 5])
	axMargx.set_ylabel('Density')

	#axMargy.set_ylim(axImg.get_ylim())
	axMargy.set_ylim([e[2], e[3]])
	axMargy.set_xlim([0, 23])
	axMargy.set_xticks([0, 5, 10, 15, 20])
	axMargy.set_xlabel('Density')

	axMargx.xaxis.set_major_formatter(NullFormatter())
	axMargy.yaxis.set_major_formatter(NullFormatter())
	plt.tight_layout()
	plt.savefig(opj(output_dir, 'pf-bayes-marg.pdf'), format='pdf', dpi=300)
