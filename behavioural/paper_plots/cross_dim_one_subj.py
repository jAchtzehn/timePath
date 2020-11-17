'''
Cross-dimensional influence
'''


from os.path import join as opj
from os.path import abspath
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import psignifit as ps
from utilFunctions import getAvgRespValues, getPsiSignInputData_p, fitPF
import pylab as pl

# ------------ options ------------
subjects = [2]

conditions = ['time', 'dist', 'dots']
clean_up_data = True            # remove RTs lower than 0.3s and higher than 2s and remove for participants 4, 5, 6 lumin in the first session
boundaries_RT = [0.3, 2.1]
priors = 'custom'                  # std or custom possible
output_filename = 'pse_data_cross_dim_individual_norm'

calc_crossDim_pse = True
plot_crossDim_pse = False

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

create_html = True

# ------------ File I/O ------------
experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
behavioral_dir = opj(experiment_dir, 'behavioural')
output_dir = abspath('/Users/jachtzehn/Documents/Medizin/thesis/figures/results/behavioural')


if calc_crossDim_pse:
	# read in RT data
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


	# init necessary list
	out_data = {'subject': [],
	            'condition_rel': [],
	            'condition_irrel': [],
	            'pse_irrel_small': [],
	            'pse_irrel_small_ci_l': [],
	            'pse_irrel_small_ci_u': [],
	            'pse_irrel_large': [],
	            'pse_irrel_large_ci_l': [],
	            'pse_irrel_large_ci_u': [],
	            'pse_diff': [],
	            'slope_irell_small': [],
	            'slope_irell_large': []
	            }

	img_list = []
	for subject in subjects:
		subj_str = 'sub-' + str(subject).zfill(2)
		print('Subject: %s' % subj_str)

		os.system('mkdir -p %s' % opj(behavioral_dir, subj_str, 'results'))
		for condition_rel in conditions:
			for condition_irrel in [x for x in conditions if x != condition_rel]:
				pse_data_irrel = []

				fig, ax = plt.subplots(1, 1, figsize=figsize)

				for idx_compVal_irrel in range(2):

					psign_input_rel = []
					for idx_compVal_rel in range(2):
						data_compVal = data_all[(data_all.trial_type == condition_rel) &
						                        (data_all['std_' + condition_rel] == std_values[condition_rel][idx_compVal_rel]) &
						                        (data_all['std_' + condition_irrel] == std_values[condition_irrel][idx_compVal_irrel]) &
						                        (data_all.subject == subject)]

						unique_comp_values = np.sort(data_compVal.compVal.unique())

						# now calculate % of 1's for each comparison value
						avg_resp_values = getAvgRespValues(unique_comp_values, data_compVal)
						psign_input_tmp = getPsiSignInputData_p(unique_comp_values, avg_resp_values)

						if condition_rel == 'dist' or condition_rel == 'dots':
							x_ref_val = std_values[condition_rel][idx_compVal_rel] / 2
						else:
							x_ref_val = std_values[condition_rel][idx_compVal_rel]

						# before fitting, reference on the mean of the standard values
						unique_comp_values = unique_comp_values / x_ref_val
						psign_input_tmp[:, 0] = psign_input_tmp[:, 0] / x_ref_val

						psign_input_rel.append(psign_input_tmp)

					psign_input_unsort = np.concatenate((psign_input_rel[0], psign_input_rel[1]), axis=0)
					psign_input = psign_input_unsort[np.argsort(psign_input_unsort[:, 0])]
					result, output_vars = fitPF(psign_input, options, priors='custom')
					[thresh, th_l, th_u] = [result['Fit'][0], result['conf_Intervals'][0][0],
					                        result['conf_Intervals'][0][1]]
					slope = ps.getSlopePC(result, 0.5)

					#ps.psigniplot.plotPsych(result, title='Rel cond: {}, irrel cond: {}, index irrel {}, '.format(condition_rel, condition_irrel, idx_compVal_irrel))
					pse_data_irrel.append([thresh, th_l, th_u, slope])

					# plot
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
					xLow = np.linspace(0.25, xMin, num=100)
					xLow = xLow[xLow >= 0]
					xHigh = np.linspace(xMax, 1.75, num=100)
					xHigh = xHigh[xHigh <= 3]

					fitValues = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](x, fit[0], fit[1]) + fit[3]
					fitValuesLow = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](xLow, fit[0], fit[1]) + fit[3]
					fitValuesHigh = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](xHigh, fit[0], fit[1]) + fit[3]

					[width, wi_l, wi_u] = [result['Fit'][1], result['conf_Intervals'][1][0],
					                       result['conf_Intervals'][1][1]]

					if idx_compVal_irrel == 0:
						ax.plot(x, fitValues, color='GoldenRod', linewidth=lw, label='after short distances (low)')
						ax.plot(xLow, fitValuesLow, '--', linewidth=lw, color='GoldenRod')
						ax.plot(xHigh, fitValuesHigh, '--', linewidth=lw, color='GoldenRod')

					else:
						ax.plot(x, fitValues, color='indianred', linewidth=lw, label='after long distances (high)')
						ax.plot(xLow, fitValuesLow, '--', linewidth=lw, color='indianred')
						ax.plot(xHigh, fitValuesHigh, '--', linewidth=lw, color='indianred')

					# threshhold x
					x = [fit[0], fit[0]]
					y = [0, .5]
					ax.plot(x, y, '-', c='slategrey', linestyle='--')

				ax.set_ylim([0, 1.05])
				ax.set_xlim([0.25, 1.75])
				ax.spines['top'].set_visible(False)
				ax.spines['right'].set_visible(False)
				ax.set_xlabel('Normalized comparison value')
				ax.set_ylabel('p(M)')
				plt.legend(loc='top left', ncol=1, fontsize=24)

				ax.annotate('PSE difference={:.3f}'.format(pse_data_irrel[1][0] - pse_data_irrel[0][0]), (1.18, 0.25), fontsize=26, color='slategray')
				fig.savefig(opj('/Users/jachtzehn/Documents/Medizin/thesis/figures/results/behavioural', 'corrDim-rel-{}-irell-{}.pdf'.format(condition_rel, condition_irrel)), dpi=300)
				plt.close()

				out_data['subject'].append(subject)
				out_data['condition_rel'].append(condition_rel)
				out_data['condition_irrel'].append(condition_irrel)
				out_data['pse_irrel_small'].append(pse_data_irrel[0][0])
				out_data['pse_irrel_small_ci_l'].append(pse_data_irrel[0][1])
				out_data['pse_irrel_small_ci_u'].append(pse_data_irrel[0][2])
				out_data['slope_irell_small'].append(pse_data_irrel[0][3])
				out_data['pse_irrel_large'].append(pse_data_irrel[1][0])
				out_data['pse_irrel_large_ci_l'].append(pse_data_irrel[1][1])
				out_data['pse_irrel_large_ci_u'].append(pse_data_irrel[1][2])
				out_data['slope_irell_large'].append(pse_data_irrel[1][3])

				out_data['pse_diff'].append(pse_data_irrel[1][0] - pse_data_irrel[0][0])

				print('Relevant condition: {}, irrelevant condition: {}, pse_diff = {:.3f}'.format(
					condition_rel, condition_irrel, pse_data_irrel[1][0] - pse_data_irrel[0][0]))

	pse_data = pd.DataFrame.from_dict(out_data)
	pse_data.to_csv(opj(behavioral_dir, output_filename + '.tsv'), sep='\t', columns=pse_data.columns, index=False)

if plot_crossDim_pse:
	data = pd.read_csv(opj(behavioral_dir, output_filename + '.tsv'), delimiter='\t')

	p_values = []
	for condition_rel in conditions:
		for condition_irrel in [x for x in conditions if x != condition_rel]:

			data_condition = data[(data.condition_rel == condition_rel) & (data.condition_irrel == condition_irrel)]

			[res_ks, p_ks] = stats.normaltest(data_condition.pse_diff)
			mean = np.mean(data_condition.pse_diff)
			[t, p] = stats.ttest_1samp(data_condition.pse_diff, 0)

			[w, p_w] = stats.wilcoxon(data_condition.pse_diff)

			#print('Norm-Test: Condition_rel: {}, Condition_irrel: {}, k-s={:.3f}, p={:.4f}'.format(
			# 	condition_rel, condition_irrel, res_ks, p_ks))

			p_values.append(p_w)

			# print('T-Test: Condition_rel: {}, Condition_irrel: {}, mean={:.3f}, t={:.4f}, p={:.4f}'.format(
			# 	condition_rel, condition_irrel, mean, t, p))

			print('W-Test: Condition_rel: {}, Condition_irrel: {}, mean={:.3f}, t={:.4f}, p={:.4f}'.format(
				condition_rel, condition_irrel, mean, w, p_w))

	result_multitest = multipletests(p_values, alpha=0.05, method='hs')
	print(result_multitest[1].T)
