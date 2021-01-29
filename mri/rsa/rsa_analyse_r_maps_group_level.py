import os
from os.path import join as opj
from os.path import abspath
from sys import platform
import pandas as pd
import itertools
from nilearn.image import load_img, index_img, new_img_like
import numpy as np
from nilearn.input_data import NiftiMasker
from nilearn.masking import _extrapolate_out_mask
import scipy.stats as stats
from tqdm import tqdm

# ------------ File I/O ------------
# experiment_dir = abspath('/home/achtzehnj/data/timePath/')
experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath')

rsa_dir = opj(experiment_dir, 'derivatives', 'rsa')
behav_file = opj(experiment_dir, 'derivatives', 'behavioural', 'pse_data_cross_dim.tsv')
behav_data = pd.read_csv(behav_file, delimiter='\t')
mask_img = load_img(opj(experiment_dir, 'derivatives', 'rsa', 'group_mask_wb_binarized.nii.gz'))

# ------------ options ------------
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
conditions = ['time', 'dist', 'dots']

cbar = tqdm(itertools.permutations(conditions, 2), leave=False)

masker = NiftiMasker(mask_img=mask_img, standardize=False, detrend=False,
                     memory=abspath('/Users/jachtzehn/data/fMRI/nilearn_cache'), memory_level=1)

for condition_pair in cbar:
	cbar.set_description('Processing condition pair {} on {}'.format(condition_pair[1], condition_pair[0]))

	# load crossDim values for all 24 subjects into ndarray
	crossDim = np.abs(behav_data['pse_diff'][(behav_data.condition_rel == condition_pair[0]) &
	                                  (behav_data.condition_irrel == condition_pair[1])].values)

	if (condition_pair[0] == 'time' and condition_pair[1] == 'dist') or (condition_pair[1] == 'time' and condition_pair[0] == 'dist'):
		dimId = 0
	elif (condition_pair[0] == 'time' and condition_pair[1] == 'dots') or (condition_pair[1] == 'time' and condition_pair[0] == 'dots'):
		dimId = 1
	else:
		dimId = 2

	# now load up each participant's RDM map and correlate
	rdm_dim_data = np.zeros((len(subjects), 64291))

	for s, subjectID in enumerate(subjects):
		rdm_map = load_img(opj(rsa_dir, 'sub-' + str(subjectID).zfill(2), 'space-MNI152NLin2009cAsym', 'results',
		                       'sub-' + str(subjectID).zfill(2) + '_space-MNI152NLin2009cAsym_wb_rdm_values.nii.gz'), )

		# extract right 4D volume from nifti file (0 = time vs. space, 1 = time vs. num, 2 = space vs. num

		rdm_map_dim = index_img(rdm_map, dimId)

		rdm_dim_data[s, :] = masker.fit_transform(rdm_map_dim)

	# average over all subjects
	rdm_dim_data_mean = np.mean(rdm_dim_data, axis=0)

	corr_data = np.zeros((1, 64291))

	a = rdm_dim_data.shape[1]

	vbar = tqdm(range(rdm_dim_data.shape[1]), leave=False, mininterval=1)
	for vx in vbar:

		[vx_corr_r, vx_corr_p] = stats.spearmanr(rdm_dim_data[:, vx], crossDim.T)

		if vx_corr_p <= 0.05:
			corr_data[0, vx] = vx_corr_r

	corr_img = masker.inverse_transform(corr_data)
	corr_img.to_filename(opj(rsa_dir, 'corrImg_rel-{}_irrel-{}.nii.gz'.format(condition_pair[0], condition_pair[1])))

