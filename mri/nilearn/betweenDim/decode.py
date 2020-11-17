"""
Decodes various labels (defined in options) of the timePath 2 data
Original script that decodes one trial type vs. the other (e.g. time trials vs. space trials) across different
anatomical ROIs.
"""

import _pickle as cPickle
from sklearn import svm
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import permutation_test_score
from sklearn.linear_model import LogisticRegression
import numpy as np
from os.path import join as opj
from os.path import abspath
from sys import platform
import pandas as pd
from sklearn.model_selection import LeaveOneGroupOut
import itertools
import os
from tqdm import tqdm
from nilearn.image import math_img
from nilearn.input_data import NiftiMasker
import datetime


# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/mnt/work/achtzehnj/data/')

nilearn_dir = opj(experiment_dir, 'nilearn')

# ------------ options ------------
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
# [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]       # all except 15
conditions_to_decode = ['time', 'dist', 'dots', 'lumin']
masks = ['ips_lh', 'ips_rh', 'ifg_lh', 'ifg_rh', 'hc_lh', 'hc_rh', 'V5_rh', 'V5_lh', 'insula_rh', 'insula_lh']
scoring = 'balanced_accuracy'                                               # see https://scikit-learn.org/stable/modules/model_evaluation.html#scoring-parameter
write_data = True                                                           # should the data be written into .pkl files
n_proc = 40                                                                 # how many processors should be used? n_proc = -1 will use all available
n_perm = 1000                                                               # number of permutations

file_suffix = 'solver-' + 'linRegr'                                         # without "_"!
space = 'MNI152NLin2009cAsym'                                               # MNI152NLin2009cAsym

# ------------ code ------------
os.system('clear')
print('\nDecoding for subjects: %s, file_suffix: %s\n' % (subjects, file_suffix))

start_time = datetime.datetime.now().replace(microsecond=0)

score_data = {}
sbar = tqdm(subjects, leave=True)
for subj in sbar:
	
	sbar.set_description('Processing subject %s' % str(subj))
	subj_data = {}     # dict in which data is stored for subj
	subj_str = 'sub-' + str(subj).zfill(2)      # create the subject string (just for easier handling)
	os.system('mkdir -p %s' % opj(nilearn_dir, subj_str, 'space-' + space, 'results'))        # make sure the results dir exists

	mbar = tqdm(masks, leave=False)
	for mask in mbar:

		mbar.set_description('Processing mask %s' % mask)
		# make results dir
		if not os.path.exists(opj(nilearn_dir, subj_str, 'space-' + space, 'results')):
			os.mkdir(opj(nilearn_dir, subj_str, 'space-' + space, 'results'))

		anat_data_name = opj(experiment_dir, 'fmriprep', subj_str, 'anat', subj_str + '_T1w_preproc.nii.gz')    # T1 image
		func_data_name = opj(nilearn_dir, subj_str, 'space-' + space, 'betas', subj_str + '_betas_merged.nii.gz')                 # betas
		if mask == 'V5_rh' or mask == 'V5_lh':
			mask_filename = opj(nilearn_dir, 'group_masks', 'space-' + space, 'group_mask_' + mask + '_binarized.nii.gz')   # mask filename
		else:
			mask_filename = opj(nilearn_dir, subj_str, 'space-' + space, 'masks',
			                    subj_str + '_' + mask + '_mask_binarized.nii.gz')  # mask filename
		behav_filename = opj(nilearn_dir, subj_str, 'space-' + space, 'betas', subj_str + '_merged_events.tsv')                   # events file

		masker = NiftiMasker(mask_img=mask_filename, standardize=True, detrend=False)       # create nifti masker (2D array ready) from mask file
		func_data = math_img("np.nan_to_num(img)", img=func_data_name)                      # before continuing, NaNs in the image must be replaced with 0

		# load behavioral information
		behavioral = pd.read_csv(behav_filename, delimiter='\t')        # read in the data of the file
		conditions_all = behavioral['trial']                            # these are the conditions, i.e. trial types

		conditions_data = {}                    # dict in which data is stored for conditions of interest
		cbar = tqdm(list(itertools.combinations(conditions_to_decode, 2)), leave=False)
		for conditions_of_interest in cbar:

			cbar.set_description('Processing condition %s vs. %s' % (conditions_of_interest[0], conditions_of_interest[1]))

			conditions_of_interest = list(conditions_of_interest)               # itertools.combinations give a tuple, we need a list here
			conditions_mask = conditions_all.isin(conditions_of_interest)       # mask to restrict the conditions to the ones specified by options
			loa_label = behavioral['class_chunk'][conditions_mask]              # mask that specifies the runs masked by the conditions

			# apply mask to fMRI data and conditions
			fmri_masked = masker.fit_transform(func_data)                       # mask the betas according to anatomical rois
			fmri_masked = fmri_masked[conditions_mask]                          # mask according to conditions of interest
			conditions = conditions_all[conditions_mask]                        # mask behavioral data

			# get the number of samples
			n_samples = [conditions.tolist().count(conditions_of_interest[0]), conditions.tolist().count(conditions_of_interest[1])]

			# decode
			svc = svm.NuSVC(kernel='linear')
			clf = LogisticRegression(C=0.1, penalty='l2', solver='lbfgs')

			clf.fit(fmri_masked, conditions)
			
			cv_score = cross_val_score(clf, fmri_masked, conditions, cv=LeaveOneGroupOut(), groups=loa_label, n_jobs=n_proc, scoring=scoring)       # using leaveonegroupout
			score, permutation_scores, pvalue = permutation_test_score(clf, fmri_masked, conditions, cv=LeaveOneGroupOut(), groups=loa_label,
			                                                           n_jobs=n_proc, n_permutations=n_perm, scoring=scoring)

			conditions_data[conditions_of_interest[0] + '-' + conditions_of_interest[1]] = [np.mean(cv_score), float(np.std(cv_score)),
			                                                                                float(np.mean(permutation_scores)), pvalue, [clf.coef_, masker],
			                                                                                [n_samples[0], n_samples[1]]]      # save data for conditions
			
			cbar.set_postfix({'Condition': conditions_of_interest, 'CV score': score, 'p': pvalue})

		subj_data[mask] = conditions_data       # save data for subject

	score_data[subj_str] = subj_data            # save data for all subjects

	if write_data:
		sbar.set_postfix_str(s='Saving individual data...')
		# individual level
		with open(opj(nilearn_dir, subj_str, 'space-' + space, 'results', subj_str + '_' + file_suffix + '_classData.pkl'), 'wb+') as f:
			cPickle.dump(subj_data, f)
			f.close()

end_time = datetime.datetime.now().replace(microsecond=0)
print('Finished in {}'.format(end_time - start_time))