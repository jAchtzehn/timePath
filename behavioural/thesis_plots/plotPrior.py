'''
Plot standard priors and adjusted priors for methods section
'''

from os.path import join as opj
from os.path import abspath
import numpy as np
import psignifit as ps
import matplotlib
import matplotlib.pyplot as plt
from psignifit import utils as _utils
from scipy.signal import convolve as _convn
from utilFunctions import widthPrior
import pylab as pl

experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
behavioral_dir = opj(experiment_dir, 'behavioural')
output_dir = abspath('/Users/jachtzehn/Documents/Medizin/thesis/figures/methods/')

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

figsize = (16, 10)
tight_pad = 5
lw = 5
ms = 30

colors = ['slategray', 'darkgoldenrod', 'lightcoral', 'cornflowerblue', 'darkseagreen']

data = np.array([[0.1,   0.0000,   90.0000],
                 [0.2,   10.0000,   90.0000],
                 [0.3,   10.0000,   90.0000],
                 [0.4,   10.0000,   90.0000],
                 [0.5,   10.0000,   90.0000],
                 [0.6,   52.0000,   90.0000],
                 [0.7,   53.0000,   90.0000],
                 [0.8,   62.0000,   90.0000],
                 [0.9,   64.0000,   90.0000],
                 [1.0,   76.0000,   90.0000],
                 [1.1,   76.0000,   90.0000],
                 [1.2,   79.0000,   90.0000],
                 [1.3,   79.0000,   90.0000],
                 [1.4,   88.0000,   90.0000],
                 [1.5,   88.0000,   90.0000],
                 [1.6,   90.0000,   90.0000],
                 [1.7,   90.0000,   90.0000],
                 [1.8,   90.0000,   90.0000],
                 [1.9,   90.0000,   90.0000],
                 [2.0,   90.0000,   90.0000]
                 ])


def plotPriors(result, output_name, mode):
	fig, ax = plt.subplots(2, 3, figsize=figsize)
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
			#x = np.linspace(min(result['X1D'][itheta]), max(result['X1D'][1]), steps)
			x = np.linspace(widthmin, 3 / Cfactor * widthmax, steps)

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

	ax[0, 0].plot(xthresh, ythresh, lw=lw, c='Black')
	ax[0, 0].set_xlim(xLimit)
	ax[0, 0].set_title('Threshold')
	ax[0, 0].set_ylabel('Density')

	ax[1, 0].scatter(data[:, 0], np.zeros(data[:, 0].shape), marker='x', s=ms, linewidth=lw * .5, color='IndianRed')

	ax[1, 0].set_ylabel('p(M)')
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
	if mode == 'std':
		xwidth = np.linspace(widthmin, 3 / Cfactor * widthmax, steps)
	else:
		xwidth = np.linspace(widthmin, 3 / Cfactor * widthmax, steps)
	ywidth = result['options']['priors'][1](xwidth)
	wwidth = _convn(np.diff(xwidth), .5 * np.array([1, 1]))
	cwidth = np.cumsum(ywidth * wwidth)

	ax[0, 1].plot(xwidth, ywidth, lw=lw, c='Black')
	# plt.hold(True)
	if mode == 'std':
		ax[0, 1].set_xlim([widthmin, 3 / Cfactor * widthmax])
	else:
		ax[0, 1].set_xlim([widthmin, 3 / Cfactor * widthmax])
	ax[0, 1].set_title('Width')
	ax[1, 1].scatter(data[:, 0], np.zeros(data[:, 0].size), marker='x', s=ms, linewidth=lw * .5, color='IndianRed')
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

	# lapse
	xlapse = np.linspace(0, .5, steps)
	ylapse = result['options']['priors'][2](xlapse)
	wlapse = _convn(np.diff(xlapse), .5 * np.array([1, 1]))
	clapse = np.cumsum(ylapse * wlapse)

	ax[0, 2].plot(xlapse, ylapse, lw=lw, c='Black')
	# plt.hold(True)
	ax[0, 2].set_xlim([0, .5])
	if mode == 'std':
		ax[0, 2].set_title('Lapse, Guess and Eta')
	else:
		ax[0, 2].set_title('Eta')

	ax[1, 2].scatter(data[:, 0], np.zeros(data[:, 0].size), marker='x', s=ms, linewidth=lw * .5, color='IndianRed')
	ax[1, 2].set_xlim(xLimit)
	ax[1, 2].set_ylim([-5, 100])

	x = np.linspace(xLimit[0], xLimit[1], steps)
	for idot in range(0, 5):
		if idot == 0:
			xcurrent = theta[2]
			color = colors[idot]
		elif idot == 1:
			xcurrent = 0
			color = colors[idot]
		elif idot == 2:
			lix = clapse[clapse >= .25].size
			xcurrent = xlapse[-lix]
			color = colors[idot]
		elif idot == 3:
			lix = clapse[clapse >= .75].size
			xcurrent = xlapse[-lix]
			color = colors[idot]
		elif idot == 4:
			xcurrent = .5
			color = colors[idot]
		y = 100 * (theta[3] + (1 - xcurrent - theta[3]) * result['options']['sigmoidHandle'](x, theta[0], theta[1]))
		ax[1, 2].plot(x, y, '-', lw=lw, c=color)
		ax[0, 2].plot(np.array(xcurrent), result['options']['priors'][2](np.array(xcurrent)), '.', c=color, ms=ms*.75)

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
		curAx[2].tick_params(axis="x", direction="out", width=1, length=5)
		curAx[2].tick_params(axis="y", direction="out", width=1, length=5)
		curAx[0].spines['top'].set_visible(False)
		curAx[0].spines['right'].set_visible(False)
		curAx[1].spines['top'].set_visible(False)
		curAx[1].spines['right'].set_visible(False)
		curAx[2].spines['top'].set_visible(False)
		curAx[2].spines['right'].set_visible(False)

	# xlabel
	fig.text(0.5, 0.03, 'Normalised comparison value', ha='center', va='center', fontsize=28)
	fig.text(0.5, 0.52, 'Parameter value', ha='center', va='center', fontsize=28)

	# legend
	# fig.legend(handles=qs, labels=['0%', '25%', '50%', '75%', '100%'], loc="lower center", ncol=4, borderaxespad=0.1, fontsize=24)
	plt.tight_layout(tight_pad)

	plt.savefig(opj(output_dir, output_name + '.pdf'), format='pdf', dpi=300)
	plt.close()

options = dict()   # initialize as an empty dict

# standard priors
options['sigmoidName'] = 'norm'
options['expType'] = 'equalAsymptote'
options['widthalpha'] = 0.2
result_std = ps.psignifit(data, options)

plotPriors(result_std, 'std-prior', 'std')

# custom priors
options_custom = result_std['options']
# get priors from that fit
options['priors'] = ps.priors.getStandardPriors(result_std['data'], result_std['options'])
# adjust priors
widthmax = (options_custom['stimulusRange'][1] - options_custom['stimulusRange'][0])
widthmin = options_custom['widthmin']
xspread = options_custom['stimulusRange'][1]-options_custom['stimulusRange'][0]

#  width
Cfactor = (ps.utils.my_norminv(.95, 0, 1) - ps.utils.my_norminv(.05, 0, 1)) / \
          (ps.utils.my_norminv(1 - options['widthalpha'], 0, 1) - ps.utils.my_norminv(
	          options['widthalpha'], 0, 1))

options['priors'][1] = lambda x: widthPrior(x, Cfactor, widthmin, widthmax)

result_custom = ps.psignifit(data, options)
plotPriors(result_custom, 'custom-prior', 'custom')


