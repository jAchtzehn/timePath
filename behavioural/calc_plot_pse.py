'''
Plotting of all PFs for each participants for visual inspection
'''


from os.path import join as opj
from os.path import abspath
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import psignifit as ps
from matplotlib.colors import to_rgb

from utilFunctions import getAvgRespValues, getPsiSignInputData_p, create_html_report, fitPF

# ------------ options ------------
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]

conditions = ['time', 'dist', 'dots', 'lumin']
clean_up_data = True            # remove RTs lower than 0.3s and higher than 2s and remove for participants 4, 5, 6 lumin in the first session
boundaries_RT = [0.3, 2.1]
priors = 'custom'                  # std or custom possible
output_filename = 'pse_data_all_common_martin'

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

create_html = True

# ------------ File I/O ------------
experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
behavioral_dir = opj(experiment_dir, 'behavioural')
output_dir = abspath('/Users/jachtzehn/Documents/Medizin/thesis/figures/results/behavioural')

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


# originally in my thesis
def calc_plot_diff_magn():
	# init necessary list
	out_data = {'subject': [],
	            'condition': [],
	            'magnitude': [],
	            'thresh': [],
	            'th_lower': [],
	            'th_upper': [],
	            'width': [],
	            'wi_lower': [],
	            'wi_upper': [],
	            'slope': [],
	            'eta': [],
	            'et_lower': [],
	            'et_upper': []
	            }

	img_list = []
	for subject in subjects:
		subj_str = 'sub-' + str(subject).zfill(2)
		print('Subject: %s' % subj_str)

		os.system('mkdir -p %s' % opj(behavioral_dir, subj_str, 'results'))
		for condition in conditions:
			for idx_compVal in range(2):

				print('Condition: %s, modality: %s' % (condition, str(idx_compVal)))

				# resp values
				data_compVal = data_all[(data_all.trial_type == condition) &
				                        (data_all['std_' + condition] == std_values[condition][idx_compVal]) &
				                        (data_all.subject == subject)]

				unique_comp_values = np.sort(data_compVal.compVal.unique())

				# now calculate % of 1's for each comparison value
				avg_resp_values = getAvgRespValues(unique_comp_values, data_compVal)
				psign_input = getPsiSignInputData_p(unique_comp_values, avg_resp_values)

				# calculate the marker size and generate marker size scale array
				max_freq = max(avg_resp_values[:, 1])
				marker_size = np.array((avg_resp_values[:, 1] / max_freq) * 200)

				color = styles[condition][1]

				if condition == 'lumin':
					label_cond = 'control'
				elif condition == 'dist':
					label_cond = 'space'
				elif condition == 'dots':
					label_cond = 'numerosity'
				else:
					label_cond = condition

				# if condition == 'dist' or condition == 'dots':
				# 	x_ref_val = std_values[condition][idx_compVal] / 2
				# else:
				x_ref_val = std_values[condition][idx_compVal]

				# before fitting, reference on the mean of the standard values
				unique_comp_values = unique_comp_values / x_ref_val
				psign_input[:, 0] = psign_input[:, 0] / x_ref_val

				# fit first
				result, output_vars = fitPF(psign_input, options, priors='custom')

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
				xLow = np.linspace(xMin - .4 * xLength, xMin, num=100)
				xLow = xLow[xLow >= 0]
				xHigh = np.linspace(xMax, xMax + .4 * xLength, num=100)
				xHigh = xHigh[xHigh <= 3]

				fitValues = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](x, fit[0], fit[1]) + fit[3]
				fitValuesLow = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](xLow, fit[0], fit[1]) + fit[3]
				fitValuesHigh = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](xHigh, fit[0], fit[1]) + fit[3]

				# plot prior
				plt.figure(figsize=(12, 8))
				ps.psigniplot.plotPrior(result)
				plt.savefig(opj(behavioral_dir, subj_str, 'results', 'sub-' + subj_str + '_condition-' +
				                condition + '_modality-' + str(idx_compVal) + '_prior.png'), dpi=300)
				plt.close()
				# plot prior

				fig, ax = plt.subplots(1, 2, figsize=(15, 8))
				fig.subplots_adjust(hspace=0.5)

				ax[0].scatter(unique_comp_values, avg_resp_values[:, 0], marker='x', color='black', label=label_cond)
				ax[0].plot(x, fitValues, color=color, clip_on=False)
				ax[0].plot(xLow, fitValuesLow, '--', color=color)
				ax[0].plot(xHigh, fitValuesHigh, '--', color=color)

				[thresh, th_l, th_u] = [result['Fit'][0], result['conf_Intervals'][0][0], result['conf_Intervals'][0][1]]
				[width, wi_l, wi_u] = [result['Fit'][1], result['conf_Intervals'][1][0], result['conf_Intervals'][1][1]]
				[eta, et_l, et_u] = [result['Fit'][4], result['conf_Intervals'][4][0], result['conf_Intervals'][4][1]]
				slope = ps.getSlopePC(result, 0.5)
				ps.plotMarginal(result, dim=0, axisHandle=ax[1], showImediate=False, lineColor=to_rgb('IndianRed'), priorColor=to_rgb('red'))
				ps.plotMarginal(result, dim=1, axisHandle=ax[1], showImediate=False, lineColor=to_rgb('SkyBlue'), priorColor=to_rgb('blue'))

				ax[0].set_xlabel('Normalised stimulus intensity')
				ax[0].set_ylabel('p(Comparison higher than perceived value)')
				ax[1].set_xlabel('Normalised stimulus intensity')

				if idx_compVal == 0:
					fig.suptitle('{}, small {}, thresh: {:.2f} ({:.2f}, {:.2f}), width: {:.2f} ({:.2f}, {:.2f}), slope: {:.2f})'.
					             format(subj_str, condition, thresh, th_l[0], th_u[0], width, wi_l[0], wi_u[0], slope))
				else:
					fig.suptitle('{}, large {}, thresh: {:.2f} ({:.2f}, {:.2f}), width: {:.2f} ({:.2f}, {:.2f}), slope: {:.2f})'.
					             format(subj_str, condition, thresh, th_l[0], th_u[0], width, wi_l[0], wi_u[0], slope))

				fig.savefig(opj(behavioral_dir, subj_str, 'results', 'sub-' + subj_str + '_condition-' +
				                condition + '_modality-' + str(idx_compVal) + '.png'), dpi=300)
				plt.close()

				img_list.append(opj(behavioral_dir, subj_str, 'results', 'sub-' + subj_str + '_condition-' +
				                condition + '_modality-' + str(idx_compVal) + '.png'))

				# append data
				out_data['subject'].append(subject)
				out_data['condition'].append(condition)
				if idx_compVal == 0:
					out_data['magnitude'].append('small')
				else:
					out_data['magnitude'].append('large')
				out_data['thresh'].append(thresh)
				out_data['th_lower'].append(float(th_l))
				out_data['th_upper'].append(float(th_u))
				out_data['width'].append(width)
				out_data['wi_lower'].append(float(wi_l))
				out_data['wi_upper'].append(float(wi_u))
				out_data['eta'].append(eta)
				out_data['et_lower'].append(float(et_l))
				out_data['et_upper'].append(float(et_u))
				out_data['slope'].append(float(slope))

	pse_data = pd.DataFrame.from_dict(out_data)
	pse_data.to_csv(opj(behavioral_dir, output_filename + '.tsv'), sep='\t', columns=pse_data.columns, index=False)

	if create_html:
		create_html_report(img_list, opj(behavioral_dir, output_filename + '.html'), '')


# for martin
def calc_plot_common_magn():
	# init necessary list
	out_data = {'subject': [],
	            'condition': [],
	            'thresh': [],
	            'th_lower': [],
	            'th_upper': [],
	            'width': [],
	            'wi_lower': [],
	            'wi_upper': [],
	            'slope': [],
	            'eta': [],
	            'et_lower': [],
	            'et_upper': []
	            }

	img_list = []
	for subject in subjects:
		subj_str = 'sub-' + str(subject).zfill(2)
		print('Subject: %s' % subj_str)

		os.system('mkdir -p %s' % opj(behavioral_dir, subj_str, 'results'))
		for condition in conditions:

			print('Condition: %s' % condition)

			# get data from both magnitudes

			psign_input_both_magn = []
			for idx_compVal in range(2):
				# resp values
				data_compVal = data_all[(data_all.trial_type == condition) &
				                        (data_all['std_' + condition] == std_values[condition][idx_compVal]) &
				                        (data_all.subject == subject)]

				unique_comp_values = np.sort(data_compVal.compVal.unique())

				# now calculate % of 1's for each comparison value
				avg_resp_values = getAvgRespValues(unique_comp_values, data_compVal)
				psign_input = getPsiSignInputData_p(unique_comp_values, avg_resp_values)

				# calculate the marker size and generate marker size scale array
				max_freq = max(avg_resp_values[:, 1])
				marker_size = np.array((avg_resp_values[:, 1] / max_freq) * 200)

				color = styles[condition][1]

				if condition == 'lumin':
					label_cond = 'control'
				elif condition == 'dist':
					label_cond = 'space'
				elif condition == 'dots':
					label_cond = 'numerosity'
				else:
					label_cond = condition

				# if condition == 'dist' or condition == 'dots':
				# 	x_ref_val = std_values[condition][idx_compVal] / 2
				# else:
				x_ref_val = std_values[condition][idx_compVal]

				# before fitting, reference on the mean of the standard values
				unique_comp_values = unique_comp_values / x_ref_val
				psign_input[:, 0] = psign_input[:, 0] / x_ref_val
				psign_input_both_magn.append(psign_input)

			# concatenate the two arrays
			psign_input = np.concatenate((psign_input_both_magn[0], psign_input_both_magn[1]))
			psign_input_sorted = psign_input[np.argsort(psign_input[:,0])]
			# fit first
			result, output_vars = fitPF(psign_input_sorted, options, priors='custom')

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
			xLow = np.linspace(xMin - .4 * xLength, xMin, num=100)
			xLow = xLow[xLow >= 0]
			xHigh = np.linspace(xMax, xMax + .4 * xLength, num=100)
			xHigh = xHigh[xHigh <= 3]

			fitValues = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](x, fit[0], fit[1]) + fit[3]
			fitValuesLow = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](xLow, fit[0], fit[1]) + fit[3]
			fitValuesHigh = (1 - fit[2] - fit[3]) * options_res['sigmoidHandle'](xHigh, fit[0], fit[1]) + fit[3]

			# plot prior
			plt.figure(figsize=(12, 8))
			ps.psigniplot.plotPrior(result)
			plt.savefig(opj(behavioral_dir, subj_str, 'results', 'sub-' + subj_str + '_condition-' +
			                condition + '_prior.png'), dpi=300)
			plt.close()
			# plot prior

			fig, ax = plt.subplots(1, 2, figsize=(15, 8))
			fig.subplots_adjust(hspace=0.5)

			ax[0].scatter(psign_input_sorted[:, 0], psign_input_sorted[:, 1] / psign_input_sorted[:, 2], marker='x', color='black', label=label_cond)
			ax[0].plot(x, fitValues, color=color, clip_on=False)
			ax[0].plot(xLow, fitValuesLow, '--', color=color)
			ax[0].plot(xHigh, fitValuesHigh, '--', color=color)

			[thresh, th_l, th_u] = [result['Fit'][0], result['conf_Intervals'][0][0],
			                        result['conf_Intervals'][0][1]]
			[width, wi_l, wi_u] = [result['Fit'][1], result['conf_Intervals'][1][0], result['conf_Intervals'][1][1]]
			[eta, et_l, et_u] = [result['Fit'][4], result['conf_Intervals'][4][0], result['conf_Intervals'][4][1]]
			slope = ps.getSlopePC(result, 0.5)
			ps.plotMarginal(result, dim=0, axisHandle=ax[1], showImediate=False, lineColor=to_rgb('IndianRed'),
			                priorColor=to_rgb('red'))
			ps.plotMarginal(result, dim=1, axisHandle=ax[1], showImediate=False, lineColor=to_rgb('SkyBlue'),
			                priorColor=to_rgb('blue'))

			ax[0].set_xlabel('Normalised stimulus intensity')
			ax[0].set_ylabel('p(Comparison higher than perceived value)')
			ax[1].set_xlabel('Normalised stimulus intensity')

			fig.suptitle(
				'{}, {}, thresh: {:.2f} ({:.2f}, {:.2f}), width: {:.2f} ({:.2f}, {:.2f}), slope: {:.2f})'.
				format(subj_str, condition, thresh, th_l[0], th_u[0], width, wi_l[0], wi_u[0], slope))

			fig.savefig(opj(behavioral_dir, subj_str, 'results', 'sub-' + subj_str + '_condition-' +
			                condition + '.png'), dpi=300)
			plt.close()

			img_list.append(opj(behavioral_dir, subj_str, 'results', 'sub-' + subj_str + '_condition-' +
			                    condition + '.png'))

			# append data
			out_data['subject'].append(subject)
			out_data['condition'].append(condition)
			out_data['thresh'].append(thresh)
			out_data['th_lower'].append(float(th_l))
			out_data['th_upper'].append(float(th_u))
			out_data['width'].append(width)
			out_data['wi_lower'].append(float(wi_l))
			out_data['wi_upper'].append(float(wi_u))
			out_data['eta'].append(eta)
			out_data['et_lower'].append(float(et_l))
			out_data['et_upper'].append(float(et_u))
			out_data['slope'].append(float(slope))

	pse_data = pd.DataFrame.from_dict(out_data)
	pse_data.to_csv(opj(behavioral_dir, output_filename + '.tsv'), sep='\t', columns=pse_data.columns, index=False)

	if create_html:
		create_html_report(img_list, opj(behavioral_dir, output_filename + '.html'), '')


calc_plot_common_magn()