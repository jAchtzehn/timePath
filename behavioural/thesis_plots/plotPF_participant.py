'''
Plot PF of a specific participant, modality and magnitude
Plot width of the PF in a shaded area beneath the curve
For Precision judgements in results folder
'''


from os.path import join as opj
from os.path import abspath
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import psignifit as ps
import pylab as pl

from utilFunctions import getAvgRespValues, getPsiSignInputData_p, widthPrior


# ------------ options ------------
subject = 3
condition = 'dots'
magnitude = 'large'
annotate_pos = (1.08, 0.4)
clean_up_data = True                # remove RTs lower than 0.3s and higher than 2s and remove for participants 4, 5, 6 lumin in the first session
boundaries_RT = [0.3, 2.1]
priors = 'custom'                   # std or custom possible

# fit options
options = dict()  # initialize as an empty dictionary
options['sigmoidName'] = 'norm'  # choose a cumulative Gauss as the sigmoid
options['expType'] = 'equalAsymptote'
options['threshPC'] = .5
options['useGPU'] = 1
options['confP'] = [.95]
options['stepN'] = [40, 40, 50, 50, 20]
options['widthalpha'] = 0.2

# threshold, width, lapse, guess, eta
options['fixedPars'] = np.array([float('nan'), float('nan'),
                                 float(.01), float(0),
                                 float('nan')])

std_values = {'time': [2.8, 4.8],
              'dist': [11.5, 19.7],
              'dots': [45.0, 77.0],
              'lumin': [0.16, 0.28]}


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


styles = {'time': ['v-', 'SkyBlue', 'full'],
          'dist': ['o:', 'IndianRed', 'full'],
          'dots': ['x--', 'GoldenRod', 'full'],
          'lumin': ['D-.', 'DimGray', 'full']}

# ------------ File I/O ------------
experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
behavioral_dir = opj(experiment_dir, 'behavioural')
output_dir = abspath('/Users/jachtzehn/Documents/Charite/talks/lab_meeting_08102020/figures')

data_all = pd.read_csv(opj(behavioral_dir, 'RT_data_all_subj.tsv'), delimiter='\t')

# clean up data
if clean_up_data:
	# delete RTs below 300 ms
	data_all = data_all[(data_all.RT >= boundaries_RT[0]) & (data_all.RT <= boundaries_RT[1])]

	# delete session 1 lumin trials for subjects 4, 5 and 6
	data_all = data_all[np.invert((data_all.subject == 4) & (data_all.session == 1) & (data_all.trial_type == 'lumin'))]
	data_all = data_all[np.invert((data_all.subject == 5) & (data_all.session == 1) & (data_all.trial_type == 'lumin'))]
	data_all = data_all[np.invert((data_all.subject == 6) & (data_all.session == 1) & (data_all.trial_type == 'lumin'))]

	# delete participant 15
	data_all = data_all[data_all.subject != 15]


def plotPF(subject, condition, magnitude, options, annotate_pos):

	if magnitude == 'large':
		idx_compVal = 1
	else:
		idx_compVal = 0

	# resp values
	data_compVal = data_all[(data_all.trial_type == condition) & (data_all['std_' + condition] == std_values[condition][idx_compVal]) & (
					data_all.subject == subject)]

	unique_comp_values = np.sort(data_compVal.compVal.unique())

	# now calculate % of 1's for each comparison value
	avg_resp_values = getAvgRespValues(unique_comp_values, data_compVal)
	psign_input = getPsiSignInputData_p(unique_comp_values, avg_resp_values)

	# calculate the marker size and generate marker size scale array
	max_freq = max(avg_resp_values[:, 1])

	if condition == 'lumin':
		label_cond = 'control'
	elif condition == 'dist':
		label_cond = 'space'
	elif condition == 'dots':
		label_cond = 'numerosity'
	else:
		label_cond = condition

	if condition == 'dist' or condition == 'dots':
		x_ref_val = std_values[condition][idx_compVal] / 2
	else:
		x_ref_val = std_values[condition][idx_compVal]

	# before fitting, reference on the mean of the standard values
	unique_comp_values = unique_comp_values / x_ref_val
	psign_input[:, 0] = psign_input[:, 0] / x_ref_val

	# fit first
	result_std_priors = ps.psignifit(psign_input, options)
	options_temp = result_std_priors['options']

	# get priors from that fit
	options['priors'] = ps.priors.getStandardPriors(result_std_priors['data'], result_std_priors['options'])
	# adjust priors
	widthmax = (options_temp['stimulusRange'][1] - options_temp['stimulusRange'][0])
	widthmin = options_temp['widthmin']
	Cfactor = (ps.utils.my_norminv(.95, 0, 1) - ps.utils.my_norminv(.05, 0, 1)) / (
			ps.utils.my_norminv(1 - options['widthalpha'], 0, 1) - ps.utils.my_norminv(options['widthalpha'], 0, 1))

	options['priors'][1] = lambda x: widthPrior(x, Cfactor, widthmin, widthmax)

	# fit again with tweaked priors
	if priors == 'custom':
		result = ps.psignifit(psign_input, options)
	else:
		result = result_std_priors

	# plot fit
	fit = result['Fit']
	data = result['data']
	options_res = result['options']

	xData = data[:, 0]
	yData = data[:, 1] / data[:, 2]
	xMin = min(xData)
	xMax = max(xData)
	xLength = xMax - xMin
	x = np.linspace(xMin, xMax, num=1000)

	# add higher and lower borders
	xLow = np.linspace(xMin - .25 * xLength, xMin, num=100)
	xLow = xLow[xLow >= 0]
	xHigh = np.linspace(xMax, xMax + .25 * xLength, num=100)
	xHigh = xHigh[xHigh <= 3]

	fitValues = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](x, fit[0], fit[1]) + fit[3]
	fitValuesLow = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](xLow, fit[0], fit[1]) + fit[3]
	fitValuesHigh = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](xHigh, fit[0], fit[1]) + fit[3]

	[width, wi_l, wi_u] = [result['Fit'][1], result['conf_Intervals'][1][0], result['conf_Intervals'][1][1]]

	fig, ax = plt.subplots(1, 1, figsize=figsize)
	ax.plot(x, fitValues, color='black', linewidth=lw, label='Fitted PF')
	ax.plot(xLow, fitValuesLow, '--', linewidth=lw, color='black')
	ax.plot(xHigh, fitValuesHigh, '--', linewidth=lw, color='black')

	# plot shaded area
	#ax.fill_between(x, 0, fitValues, where=(fitValues <= 0.8), alpha=1, facecolor='slategray')
	#ax.fill_between(x, 0, fitValues, where=(fitValues <= 0.2), alpha=1, facecolor='white', edgecolor='white', linewidth=2)
	#ax.fill_between(x, 0, fitValues, where=(fitValues <= 0.12), alpha=1, facecolor='white', edgecolor='White', linewidth=4)

	ax.scatter(unique_comp_values, avg_resp_values[:, 0], marker='x', s=200, linewidth=lw, color='IndianRed', label='Input data')
	# labels
	ax.set_xlabel('Normalised comparison value')
	ax.set_ylabel('p(M)')
	ax.set_ylim([0, 1.05])

	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)

	# annotate
	ax.annotate('Width={:.2f}'.format(width), annotate_pos, fontsize=26, color='slategray')

	fig.savefig(opj(output_dir, str(subject) + '-' + condition + '-' + magnitude + '.pdf'), dpi=300)


plotPF(subject, condition, magnitude, options, annotate_pos)