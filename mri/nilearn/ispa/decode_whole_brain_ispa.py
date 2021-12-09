"""
Decodes various labels (defined in options) of the timePath 2 data
"""

from sklearn import svm
from sklearn.model_selection import LeaveOneGroupOut
from nilearn.image import index_img, mean_img, new_img_like, concat_imgs, load_img, clean_img, smooth_img
from nilearn._utils.niimg_conversions import _safe_get_data
import numpy as np
from sklearn.linear_model import LogisticRegression
from nilearn.decoding import SearchLight
from os.path import join as opj
from os.path import abspath
from sys import platform
import pandas as pd
import os
from tqdm import tqdm
import datetime

# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
	nilearn_cache_dir = abspath('/home/achtzehnj/code/nilearn_utilities/cache')
else:
	experiment_dir = abspath('/home/achtzehnj/data/timePath/derivatives')
	nilearn_cache_dir = abspath('/home/achtzehnj/code/nilearn_utilities/cache')

nilearn_dir = opj(experiment_dir, 'nilearn')
results_dir = opj(nilearn_dir, 'group_results_ispa_std')

# ------------ options ------------
settings = {'subjects': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
            'conditions_to_decode': [['speed']],
            'speed_conditions_to_decode': ['high', 'low'],
            'n_proc': -1,
            'roi': 'wb',
            'standardize': True,
            'solver': 'lbfgs',
            'verbose_decode': 0,
            'smooth_img': False,
            'fwhm': 6}

# ------------ code ------------
# os.system('clear')
print('\nWhole brain ISPA decoding \n')

# print settings in console
for key in settings.keys():
	print('Option: {}, value: {}'.format(key, settings[key]))
print('\n')
print('Output dir: {}\n'.format(results_dir))

start_time = datetime.datetime.now().replace(microsecond=0)

cbar = tqdm(settings['conditions_to_decode'], leave=False)

# mask image
wb_mask_img = load_img(opj(nilearn_dir, 'group_masks', 'space-MNI152NLin2009cAsym', 'group_mask_' + settings['roi'] + '_binarized.nii.gz'))

for conditions_to_decode in cbar:

	if len(conditions_to_decode) > 1:
		cbar.set_description('Processing conditions {} vs. {}'.format(conditions_to_decode[0], conditions_to_decode[1]))
		# 0. setup results folder
		condition_results_dir = opj(results_dir, '{}_vs_{}'.format(conditions_to_decode[0], conditions_to_decode[1]))
	else:
		cbar.set_description('Processing conditions {}'.format(conditions_to_decode))
		# 0. setup results folder
		condition_results_dir = opj(results_dir, '{}'.format(conditions_to_decode[0]))
	os.system('mkdir -p {}'.format(condition_results_dir))

	# 1. load behavioral data into one df, create list of functional files and create functional data
	behav_data = pd.DataFrame()
	func_data_list = []
	for subj in settings['subjects']:
		subj_str = 'sub-' + str(subj).zfill(2)  # create the subject string (just for easier handling)
		subj_behav_data = pd.read_csv(opj(nilearn_dir, subj_str, 'space-MNI152NLin2009cAsym', 'betas',
		                                  subj_str + '_merged_events.tsv'), delimiter='\t')
		behav_data = behav_data.append(subj_behav_data, ignore_index=True)

		subj_func_data = clean_img(opj(nilearn_dir, subj_str, 'space-MNI152NLin2009cAsym', 'betas', subj_str + '_betas_merged.nii.gz'),
		                           standardize=settings['standardize'], ensure_finite=True, detrend=False, mask_img=None)
		func_data_list.append(subj_func_data)

	fmri_data = concat_imgs(func_data_list)

	# 2. mask functional and behavioral data according to conditions of interest
	if len(conditions_to_decode) > 1:
		conditions_mask = behav_data.trial.isin(conditions_to_decode)
	else:
		if conditions_to_decode[0] is not 'speed':
			conditions_mask = behav_data.trial.isin(['time', 'dist', 'dots', 'lumin'])
		else:
			conditions_mask = np.logical_and(behav_data.trial.isin(['time', 'dist', 'dots', 'lumin']), behav_data.speed.isin(settings['speed_conditions_to_decode']))

	fmri_data = index_img(fmri_data, conditions_mask)
	behav_data = behav_data[conditions_mask]
	mean_fmri = mean_img(fmri_data)
	behav_data = behav_data.reset_index(drop=True)

	# smooth
	if settings['smooth_img']:
		fmri_data = smooth_img(fmri_data, settings['fwhm'])
	# create decoding labels
	if len(conditions_to_decode) > 1:
		decode_labels = behav_data.trial
	else:
		decode_labels = behav_data[conditions_to_decode]

	# 3. setup leave-one-subject-out CV
	loso = LeaveOneGroupOut()
	n_splits = loso.get_n_splits(groups=behav_data.subject)

	# 4. setup classification
	svc = svm.NuSVC(kernel='linear')
	clf = LogisticRegression(C=0.1, penalty='l2', solver=settings['solver'])
	process_mask = _safe_get_data(wb_mask_img).astype(np.int)

	cv_bar = tqdm(loso.split(behav_data.subject, behav_data.subject, behav_data.subject), leave=False)
	for split_idx, (train_idx, test_idx) in enumerate(cv_bar):
		cv_bar.set_description('Processing split {:02d}/{:02d}, nr. of voxels {}/{}'.format(
			split_idx + 1, n_splits, np.count_nonzero(process_mask), np.count_nonzero(_safe_get_data(wb_mask_img).astype(np.int))))

		single_split = [(train_idx, test_idx)]
		searchlight = SearchLight(
			mask_img=wb_mask_img,
			radius=6, n_jobs=settings['n_proc'],
			verbose=settings['verbose_decode'], cv=single_split,
			estimator=clf, scoring='balanced_accuracy')
		searchlight.fit(fmri_data, decode_labels)

		if len(conditions_to_decode) > 1:
			img_name = 'ispa_sl_wb_decoding-{}-vs-{}_split-{:02d}.nii'.format(
				conditions_to_decode[0], conditions_to_decode[1], split_idx)
		else:
			if conditions_to_decode[0] is not 'speed':
				img_name = 'ispa_sl_wb_decoding-{}_split-{:02d}.nii'.format(
					conditions_to_decode[0], split_idx)
			else:
				img_name = 'ispa_sl_wb_decoding-speed-{}-vs-{}_split-{:02d}.nii'.format(
					settings['speed_conditions_to_decode'][0], settings['speed_conditions_to_decode'][1], split_idx)
		searchlight_img = new_img_like(mean_fmri, searchlight.scores_ - 0.5)
		searchlight_img.to_filename(opj(condition_results_dir, img_name))

end_time = datetime.datetime.now().replace(microsecond=0)
print('Finished in {}'.format(end_time - start_time))
