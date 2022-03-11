def readSubjData(subjects, dataType, conditions):

	from os.path import join as opj
	from os.path import abspath
	import pandas as pd
	import numpy as np

	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
	behavioral_dir = opj(experiment_dir, 'behavioural')

	output_data_subj_wise = []
	subj_list = []
	for subjNr in subjects:
		sub_str = 'sub-' + str(subjNr).zfill(2)  # for readability

		subj_logfile = pd.read_csv(opj(behavioral_dir, sub_str, sub_str + '_merged_events.tsv'), delimiter='\t')

		# init dictionary
		subj_data = {key: [] for key in conditions}

		for condition in conditions:
			for runNr in range(1, 9):
				run_mask = (subj_logfile['run'] == runNr) & (subj_logfile['trial'] == condition) & (subj_logfile['RT'] >= 0.2)
				run_data = np.array(subj_logfile[dataType][run_mask])
				subj_data[condition].append(run_data)
		subj_list.append([subjNr] * len(subj_data))

		output_data_subj_wise.append(subj_data)

	return output_data_subj_wise


def readVASdata(subjects, conditions):
	from os.path import join as opj
	from os.path import abspath
	import pandas as pd
	import numpy as np

	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
	behavioral_dir = opj(experiment_dir, 'behavioural')

	output_data = []

	for subjNr in subjects:
		sub_str = 'sub-' + str(subjNr).zfill(2)  # for readability
		subj_logfile = pd.read_csv(opj(behavioral_dir, sub_str, sub_str + '_vas_responses.tsv'), delimiter='\t')

		subj_data = {}
		for condition in conditions:
			subj_data[condition] = np.array(subj_logfile['VAS_' + condition])

		output_data.append(subj_data)

	return output_data

def getAvgRespValues(unique_comp_values, data_compVal):
	import numpy as np

	avg_resp_values = np.empty([len(unique_comp_values), 2])
	for i, comp_val in enumerate(unique_comp_values):
		resp_data_compVal = data_compVal[data_compVal.compVal == comp_val]
		avg_resp_values[i] = np.array([1 - np.mean(resp_data_compVal.resp), len(resp_data_compVal)])

	return avg_resp_values


def getPsiSignInputData(unique_comp_values, data_compVal, condition):
	import numpy as np

	output_list = np.empty([len(unique_comp_values), 3])
	for i, comp_val in enumerate(unique_comp_values):

		n_correct_resp_values = 0
		resp_data_compVal = data_compVal[data_compVal.compVal == comp_val]
		for idx, row in resp_data_compVal.iterrows():
			if row['std_' + condition] < row.compVal and row.resp == 0:
				n_correct_resp_values += 1
			elif row['std_' + condition] > row.compVal and row.resp == 1:
				n_correct_resp_values += 1

		comp_list = np.array([comp_val, n_correct_resp_values, len(resp_data_compVal)])
		output_list[i] = comp_list

	return output_list


def getPsiSignInputData_p(unique_comp_values, avgValues):
	import numpy as np

	output_list = np.empty([len(unique_comp_values), 3])
	for i, comp_val in enumerate(unique_comp_values):

		comp_list = np.array([comp_val, avgValues[i][0] * avgValues[i][1], avgValues[i][1]])
		output_list[i] = comp_list

	return output_list


def create_html_report(image_list, output_filename, title=''):
	import base64
	"""
	Currently does not work in python3
	"""

	f = open(output_filename, 'w')
	f.write('<html>\n')
	f.write('<h1 style="font-family:helvetica;">')
	if title != '':
		f.write('<head><title>' + title + '</title></head>\n')
	f.write('<body><p><font size="14">' + title + '</font></p></body>\n')
	if title != '':
		f.write('<hr>')

	for img in image_list:
		image = base64.b64encode(open(img, 'rb').read()).decode().replace('\n', '')
		f.write('<img src="data:image/png;base64,{0}" width="1650">'.format(image))
		f.write('<hr>')
		f.write('\n')

	f.write('</h1>')
	f.write('</html>\n')
	f.close()


def widthPrior(x, Cfactor, wmin, wmax):

	import numpy as np

	# r = ((x * Cfactor) >= wmin) * ((x * Cfactor) <= 2 * wmin) * (
	# 			.5 - .5 * np.cos(np.pi * ((x * Cfactor) - wmin) / wmin)) \
	#     + ((x * Cfactor) > 2 * wmin) * ((x * Cfactor) < wmax) \
	#     + ((x * Cfactor) >= wmax) * ((x * Cfactor) <= 3 * wmax) * (
	# 			    .5 + .5 * np.cos(np.pi / 2 * (((x * Cfactor) - wmax) / wmax)))

	r = ((x * Cfactor) >= wmin * 5) * ((x * Cfactor) <= 5 * wmin) * (
			.5 - .5 * np.cos(np.pi / 2 * ((x * Cfactor) - wmin * 5) / wmin)) \
	    + ((x * Cfactor) > 4 * wmin) * ((x * Cfactor) < wmax * 5) \
	    + ((x * Cfactor) >= wmax * 5) * ((x * Cfactor) <= 6 * wmax) * (
			    .5 + .5 * np.cos(np.pi / 0.9 * (((x * Cfactor) - wmax * 5) / wmax)))

	return r


def fitPF(input_data, options, priors='std'):

	import psignifit as ps

	result_std_priors = ps.psignifit(input_data, options)
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
		result = ps.psignifit(input_data, options)
	else:
		result = result_std_priors

	[thresh, th_l, th_u] = [result['Fit'][0], result['conf_Intervals'][0][0], result['conf_Intervals'][0][1]]
	[width, wi_l, wi_u] = [result['Fit'][1], result['conf_Intervals'][1][0], result['conf_Intervals'][1][1]]
	[eta, et_l, et_u] = [result['Fit'][4], result['conf_Intervals'][4][0], result['conf_Intervals'][4][1]]
	slope = ps.getSlopePC(result, 0.5)

	output_var = {'thresh': [thresh, th_l, th_u],
	              'width': [width, wi_l, wi_u],
	              'eta': [eta, et_l, et_u],
	              'slope': slope}

	return result, output_var
