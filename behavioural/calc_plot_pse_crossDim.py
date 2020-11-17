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
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import psignifit as ps
from utilFunctions import getAvgRespValues, getPsiSignInputData_p, fitPF

# ------------ options ------------
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]

conditions = ['time', 'dist', 'dots', 'lumin']
clean_up_data = True            # remove RTs lower than 0.3s and higher than 2s and remove for participants 4, 5, 6 lumin in the first session
boundaries_RT = [0.3, 2.1]
priors = 'custom'                  # std or custom possible
output_filename = 'pse_data_cross_dim'

calc_crossDim_pse = True
plot_crossDim_pse = True

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

# update font sizes
matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams['ytick.labelsize'] = 18
matplotlib.rcParams['xtick.labelsize'] = 18

styles = {'time': ['v-', 'SkyBlue', 'full'],
          'dist': ['o:', 'IndianRed', 'full'],
          'dots': ['x--', 'GoldenRod', 'full'],
          'lumin': ['D-.', 'DimGray', 'full']}

create_html = False

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
	            'pse_diff': []
	            }

	img_list = []
	for subject in subjects:
		subj_str = 'sub-' + str(subject).zfill(2)
		print('Subject: %s' % subj_str)

		os.system('mkdir -p %s' % opj(behavioral_dir, subj_str, 'results'))
		for condition_rel in conditions:
			for condition_irrel in [x for x in conditions if x != condition_rel]:

				pse_data_rel = []
				for idx_compVal_rel in range(2):
					pse_data_irrel = []
					for idx_compVal_irrel in range(2):
						# resp values
						data_compVal = data_all[(data_all.trial_type == condition_rel) &
						                        (data_all['std_' + condition_rel] == std_values[condition_rel][idx_compVal_rel]) &
						                        (data_all['std_' + condition_irrel] == std_values[condition_irrel][idx_compVal_irrel]) &
						                        (data_all.subject == subject)]

						unique_comp_values = np.sort(data_compVal.compVal.unique())

						# now calculate % of 1's for each comparison value
						avg_resp_values = getAvgRespValues(unique_comp_values, data_compVal)
						psign_input = getPsiSignInputData_p(unique_comp_values, avg_resp_values)

						if condition_rel == 'dist' or condition_rel == 'dots':
							x_ref_val = np.mean([std_values[condition_rel][0] / 2, std_values[condition_rel][1] / 2])
						else:
							x_ref_val = np.mean(std_values[condition_rel])

						# before fitting, reference on the mean of the standard values
						unique_comp_values = unique_comp_values / x_ref_val
						psign_input[:, 0] = psign_input[:, 0] / x_ref_val

						result, output_vars = fitPF(psign_input, options, priors='custom')
						ps.psigniplot.plotPsych(result, title='Relevant: {}, irrelevant: {}, num: {}'.format(condition_rel, condition_irrel, np.sum(psign_input[:, 2])))
						plt.show()
						pse_data_irrel.append(output_vars['thresh'][0])

					pse_data_rel.append(pse_data_irrel[1] - pse_data_irrel[0])

				out_data['subject'].append(subject)
				out_data['condition_rel'].append(condition_rel)
				out_data['condition_irrel'].append(condition_irrel)
				out_data['pse_diff'].append(np.mean(pse_data_rel))

				print('Relevant condition: {}, irrelevant condition: {}, pse_diff = {:.3f}'.format(
					condition_rel, condition_irrel, np.mean(pse_data_rel)))

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

			# print('Condition_rel: {}, Condition_irrel: {}, k-s={:.3f}, p={:.4f}'.format(
			# 	condition_rel, condition_irrel, res_ks, p_ks))

			p_values.append(p_w)

			# print('T-Test: Condition_rel: {}, Condition_irrel: {}, mean={:.3f}, t={:.4f}, p={:.4f}'.format(
			# 	condition_rel, condition_irrel, mean, t, p))

			print('W-Test: Condition_rel: {}, Condition_irrel: {}, mean={:.3f}, t={:.4f}, p={:.4f}'.format(
				condition_rel, condition_irrel, mean, w, p_w))

	result_multitest = multipletests(p_values, alpha=0.05, method='hs')
	print(result_multitest[1].T)
