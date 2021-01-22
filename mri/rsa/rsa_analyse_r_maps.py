import os
from os.path import join as opj
from os.path import abspath
from sys import platform
import pandas as pd
import itertools
from nilearn.image import load_img, index_img
import numpy as np
from nilearn.input_data import NiftiMasker
from nilearn.masking import _extrapolate_out_mask

# ------------ File I/O ------------
# experiment_dir = abspath('/home/achtzehnj/data/timePath/')
experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath')

rsa_dir = opj(experiment_dir, 'derivatives', 'rsa')
behav_file = opj(experiment_dir, 'derivatives', 'behavioural', 'pse_data_cross_dim.tsv')
behav_data = pd.read_csv(behav_file, delimiter='\t')

# ------------ options ------------
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
conditions = ['time', 'dist', 'dots']

for condition_pair in itertools.permutations(conditions, 2):

	# load crossDim values for all 24 subjects into ndarray
	crossDim = behav_data['pse_diff'][(behav_data.condition_rel == condition_pair[0]) &
	                                  (behav_data.condition_irrel == condition_pair[1])].values

	if (condition_pair[0] == 'time' and condition_pair[1] == 'dist') or (condition_pair[1] == 'time' and condition_pair[0] == 'dist'):
		dimId = 0
	elif (condition_pair[0] == 'time' and condition_pair[1] == 'dots') or (condition_pair[1] == 'time' and condition_pair[0] == 'dots'):
		dimId = 1
	else:
		dimId = 2

	# now load up each participant's RDM map and correlate
	rdm_dim_data = np.zeros((24, 64291))

	for s, subjectID in enumerate(subjects):
		rdm_map = load_img(opj(rsa_dir, 'sub-' + str(subjectID).zfill(2), 'space-MNI152NLin2009cAsym', 'results',
		                       'sub-' + str(subjectID).zfill(2) + '_space-MNI152NLin2009cAsym_wb_rdm_values.nii.gz'), )

		# extract right 4D volume from nifti file (0 = time vs. space, 1 = time vs. num, 2 = space vs. num

		rdm_map_dim = index_img(rdm_map, dimId)


		masker = NiftiMasker(standardize=False, detrend=False, memory=abspath('/Users/jachtzehn/data/fMRI/nilearn_cache'))
		rdm_dim_data = masker.fit_transform(rdm_map_dim)
		print('Subject: {}, shape {}'.format(str(s), rdm_dim_data.shape))

	print('done')