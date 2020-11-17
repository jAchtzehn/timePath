from os.path import join as opj
from os.path import abspath
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
from seaborn import regplot


def plt_regression(x, y, label_1, label_2, title):

	reg_plot = regplot(x=x, y=y, fit_reg=True)
	plt.xlabel(label_1)
	plt.ylabel(label_2)
	plt.suptitle(title)
	plt.show()

# ------------ File I/O ------------
experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
nilearn_dir = opj(experiment_dir, 'nilearn')
rsa_dir = opj(experiment_dir, 'rsa')
behavioral_dir = opj(experiment_dir, 'behavioural')

# ------------ options ------------
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
conditions = ['time', 'dist', 'dots']
roi_names = ['ips_lh', 'ips_rh', 'ifg_lh', 'ifg_rh', 'hc_lh', 'hc_rh', 'V5_rh', 'V5_lh', 'insula_rh', 'insula_lh']

# read in behavioural data
behav_data = pd.read_csv(opj(behavioral_dir, 'pse_data_cross_dim_individual_norm' + '.tsv'), delimiter='\t')

# 1. correlation between RSA/MVPA and crossDim influence
rsa_data = pd.read_csv(opj(rsa_dir, 'group_results_space-MNI152NLin2009cAsym', 'RDM_space-MNI152NLin2009cAsym_searchlight.csv'), delimiter='\t')
mvpa_data = pd.read_csv(opj(nilearn_dir, 'group_results_space-MNI152NLin2009cAsym', 'group_level_solver-linRegr_cv_scores.csv'), delimiter='\t')

for condition_pair in itertools.permutations(conditions, 2):

	cD_dm_data = {}
	cD_mvpa_data = {}
	for roi_init in roi_names:
		cD_dm_data[roi_init] = np.zeros([len(subjects), 2])
		cD_mvpa_data[roi_init] = np.zeros([len(subjects), 2])

	for s, subjectId in enumerate(subjects):

		# crossDim = np.abs(behav_data['pse_diff'][(behav_data.subject == subjectId) &
		#                                          (((behav_data.condition_rel == condition_pair[0]) & (behav_data.condition_irrel == condition_pair[1])) |
		#                                           ((behav_data.condition_rel == condition_pair[1]) & (behav_data.condition_irrel == condition_pair[0])))].values)
		# crossDim = np.mean(crossDim)

		crossDim = np.abs(behav_data['pse_diff'][(behav_data.subject == subjectId) &
		                                         (behav_data.condition_rel == condition_pair[0]) & (
					                                         behav_data.condition_irrel == condition_pair[1])].values[0])

		for roi in roi_names:
			dm = rsa_data['dissimilarity'][((rsa_data.roi == roi) & (rsa_data.subject == subjectId)) &
										   (((rsa_data.condition_1 == condition_pair[0]) & (rsa_data.condition_2 == condition_pair[1])) |
											((rsa_data.condition_1 == condition_pair[1]) & (rsa_data.condition_2 == condition_pair[0])))].values[0]

			cD_dm_data[roi][s, 0] = crossDim
			cD_dm_data[roi][s, 1] = dm

			cv_score = mvpa_data['accuracy'][((mvpa_data.roi == roi) & (mvpa_data.subject == subjectId)) &
											 (((mvpa_data.condition_1 == condition_pair[0]) & (mvpa_data.condition_2 == condition_pair[1])) |
											  ((mvpa_data.condition_1 == condition_pair[1]) & (mvpa_data.condition_2 == condition_pair[0])))].values[0]

			cD_mvpa_data[roi][s, 0] = crossDim
			cD_mvpa_data[roi][s, 1] = cv_score

	# go through rois and test

	p_vals = []
	for roi_corr in roi_names:

		# RSA
		data_rsa = cD_dm_data[roi_corr]
		[r_rsa, p_rsa] = stats.spearmanr(data_rsa[:, 0], data_rsa[:, 1])
		p_vals.append(p_rsa)

		if p_rsa <= .05:
			print('RSA: ROI: {}, condition pair: {} on {}, r = {:.3f}, p = {:.4f} ***'.format(roi_corr, condition_pair[1], condition_pair[0], r_rsa, p_rsa))
		else:
			print(
				'RSA: ROI: {}, condition pair: {} on {}, r = {:.3f}, p = {:.4f}'.format(roi_corr, condition_pair[1],
				                                                                            condition_pair[0], r_rsa,
				                                                                            p_rsa))
			#plt_regression(data_rsa[:, 0], data_rsa[:, 1], 'crossDim', 'DM', 'ROI: {}, {} on {}'.format(roi_corr, condition_pair[1], condition_pair[0]))

		# MVPA
		# data_mvpa = cD_mvpa_data[roi_corr]
		# [r_mvpa, p_mvpa] = stats.spearmanr(data_mvpa[:, 0], data_mvpa[:, 1])
		# if p_mvpa <= 0.05:
		# 	print(
		# 		'MVPA: ROI: {}, condition pair: {} on {}, r = {:.3f}, p = {:.4f}'.format(roi_corr, condition_pair[1],
		# 		                                                                         condition_pair[0], r_mvpa,
		# 		                                                                         p_mvpa))

			#plt_regression(data_mvpa[:, 0], data_mvpa[:, 1], 'crossDim', 'Accuracy', 'ROI {}, {} on {}'.format(roi_corr, condition_pair[1], condition_pair[0]))

	#p_corr = list(multipletests(p_vals, alpha=0.05, method='hs')[1])
	#print(p_corr)