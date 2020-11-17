"""
Decoding different modalities with separate training and testing data sets
"""

import _pickle as cPickle
from sklearn import svm
from os.path import join as opj
from os.path import abspath
from sys import platform
import pandas as pd
import itertools
import os
from tqdm import tqdm
from nilearn.image import math_img
from nilearn.input_data import NiftiMasker

# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/media/sf_data/fMRI/timePath/')

nilearn_dir = opj(experiment_dir, 'nilearn')

# ------------ options ------------
subjects = range(1, 26)
# [4, 6, 7, 9, 12, 14, 17, 19, 23, 24]          # 10 best subjects (TSNR >= 50)
conditions_to_decode = ['high', 'low']          # conditions names?
train_trials_all = [['time', 'duration']]
predict_trials_all = [['dots', 'nrdots']]

masks = ['ips_lh', 'ips_rh', 'hc_lh', 'hc_rh', 'ifg_lh', 'ifg_rh', 'pcg_lh', 'pcg_rh']
scoring = 'balanced_accuracy'                                               # see https://scikit-learn.org/stable/modules/model_evaluation.html#scoring-parameter
write_data = True                                                          # should the data be written into .pkl files
n_proc = 4                                                                 # how many processors should be used? n_proc = -1 will use all available
n_perm = 1000                                                               # number of permutations

space = 'T1w'                                                   # MNI152NLin2009cAsym

# ------------ code ------------
os.system('clear')
print('\nDecoding for subjects: %s\n' % subjects)

comb_bar = tqdm(range(len(train_trials_all)), leave=True)

for comb_idx in comb_bar:
	
	train_trials = train_trials_all[comb_idx]
	predict_trials = predict_trials_all[comb_idx]
	
	file_suffix = 'score-' + 'balanced-accuracy' + '_train-' + train_trials[0] + '-' + train_trials[1] \
	              + '_predict-' + predict_trials[0] + '-' + predict_trials[1]
	
	comb_bar.set_description('Train: %s; Predict: %s' % (train_trials, predict_trials))
	
	score_data = {}
	sbar = tqdm(subjects, leave=False)
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
			mask_filename = opj(nilearn_dir, subj_str, 'space-' + space, 'masks', subj_str + '_' + mask + '_mask_binarized.nii.gz')   # mask filename
			behav_filename = opj(nilearn_dir, subj_str, 'space-' + space, 'betas', subj_str + '_merged_events.tsv')                   # events file
	
			masker = NiftiMasker(mask_img=mask_filename, standardize=True, detrend=False)       # create nifti masker (2D array ready) from mask file
			func_data = math_img("np.nan_to_num(img)", img=func_data_name)                      # before continuing, NaNs in the image must be replaced with 0
	
			# load behavioral information
			behavioral = pd.read_csv(behav_filename, delimiter='\t')        # read in the data of the file
			conditions_all = behavioral['trial']                            # these are the conditions, i.e. trial types
			train_all = behavioral[train_trials[1]]                         # trial duration
			predict_all = behavioral[predict_trials[1]]                     # trial distance, either 'short' or 'long'
		
			conditions_data = {}                    # dict in which data is stored for conditions of interest
			cbar = tqdm(list(itertools.combinations(conditions_to_decode, 2)), leave=False)
			for conditions_of_interest in cbar:
	
				cbar.set_description('Processing condition %s vs. %s' % (conditions_of_interest[0], conditions_of_interest[1]))
	
				conditions_mask_train = conditions_all.isin([train_trials][0])           # mask all the train trial types (e.g. all time trials)
				conditions_mask_predict = conditions_all.isin([predict_trials][0])       # mask all the train trial types (e.g. all dist trials)
				
				# apply mask to fMRI data and conditions
				fmri_masked = masker.fit_transform(func_data)                           # mask the betas according to anatomical rois
				fmri_masked_train = fmri_masked[conditions_mask_train]                  # mask according to conditions of interest
				fmri_masked_predict = fmri_masked[conditions_mask_predict]              # mask according to conditions of interest
				conditions_train = train_all[conditions_mask_train]                     # mask behavioral data
				conditions_predict = predict_all[conditions_mask_predict]               # mask behavioral data
				
				n_samples = [len(conditions_train), len(conditions_predict)]
	
				# decode
				svc = svm.NuSVC(kernel='linear')
	
				svc.fit(fmri_masked_train, conditions_train)
				
				prediction = svc.predict(fmri_masked_predict)
				score = (prediction == conditions_predict).sum() / float(len(conditions_predict))
	
				conditions_data[conditions_of_interest[0] + '-' + conditions_of_interest[1]] = [score, [n_samples]]      # save data for conditions
	
				mbar.set_postfix({'CV score': score})
	
			subj_data[mask] = conditions_data       # save data for subject
	
		score_data[subj_str] = subj_data            # save data for all subjects
	
		if write_data:
			sbar.set_postfix_str(s='Saving individual data...')
			# individual level
			with open(opj(nilearn_dir, subj_str, 'space-' + space, 'results', subj_str + '_' + file_suffix + '_classData.pkl'), 'wb+') as f:
				cPickle.dump(subj_data, f)
				f.close()
			
			sbar.set_postfix_str(s=' ')
		