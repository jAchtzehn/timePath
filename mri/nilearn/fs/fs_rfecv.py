import cPickle
import os
import numpy as np
from os.path import join as opj
from os.path import abspath
from sys import platform
import pandas as pd
from sklearn import svm
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.feature_selection import RFECV, RFE
import itertools
from nilearn.image import math_img
from nilearn.input_data import NiftiMasker
import matplotlib.pyplot as plt
import matplotlib
from util_functions import create_html_report
from tqdm import tqdm

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams['ytick.labelsize'] = 14
matplotlib.rcParams['xtick.labelsize'] = 14

# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/media/sf_data/fMRI/timePath/')
nilearn_dir = opj(experiment_dir, 'nilearn')

# ------------ options ------------
subjects = [4]
conditions_to_decode = ['time', 'dist', 'lumin', 'dots']
masks = ['ips_lh', 'ips_rh', 'hc_lh', 'hc_rh', 'insula_lh', 'insula_rh', 'ba6_lh', 'ba6_rh', 'ifg_lh', 'ifg_rh']

conditions_combinations = list(itertools.combinations(conditions_to_decode, 2))         # create all possible unique combinations

image_list = []                 # this stores the filenames of the images so a html report can be created
os.system('clear')

sbar = tqdm(subjects)
for subj in sbar:
	
	sbar.set_description('Processing subject: %s' % subj)
	
	subj_str = 'sub-' + str(subj).zfill(2)      # create the subject string (just for easier handling)
	os.system('mkdir -p %s' % opj(nilearn_dir, subj_str, 'results'))        # make sure the results dir exists

	fig, ax = plt.subplots(len(masks)/2, 2, figsize=(35, 20))
	plt.tight_layout(5)
	fig.subplots_adjust(hspace=0.5)
	plt_enum = 0
	row_enum = 0
	
	mbar = tqdm(masks)
	for mask in mbar:
		mbar.set_description('Processing mask: %s' % mask)
		
		anat_data_name = opj(experiment_dir, 'fmriprep', subj_str, 'anat', subj_str + '_T1w_preproc.nii.gz')    # T1 image
		func_data_name = opj(nilearn_dir, subj_str, 'betas', subj_str + '_betas_merged.nii.gz')                 # betas
		mask_filename = opj(nilearn_dir, subj_str, 'masks', subj_str + '_' + mask + '_mask_binarized.nii.gz')   # mask filename
		behav_filename = opj(nilearn_dir, subj_str, 'betas', subj_str + '_merged_events.tsv')                   # events file
		
		# behav data
		behavioral = pd.read_csv(behav_filename, delimiter='\t')        # read in the data of the file
		conditions_all = behavioral['trial']                            # these are the conditions, i.e. trial types
		
		conditions_data = {}
		
		cbar = tqdm(conditions_combinations)
		for conditions_of_interest in cbar:
			cbar.set_description('Processing conditions: %s' % conditions_of_interest[0] + '-' + conditions_of_interest[1])
			
			conditions_of_interest = list(conditions_of_interest)               # itertools.combinations give a tuple, we need a list here
			conditions_mask = conditions_all.isin(conditions_of_interest)       # mask to restrict the conditions to the ones specified by options
			loa_label = behavioral['class_chunk'][conditions_mask]              # mask that specifies the runs masked by the conditions
			conditions = conditions_all[conditions_mask]                        # mask behavioral data
			
			masker = NiftiMasker(mask_img=mask_filename, standardize=True, detrend=False)  # create nifti masker (2D array ready) from mask file
			func_data = math_img("np.nan_to_num(img)", img=func_data_name)                  # before continuing, NaNs in the image must be replaced with 0
			
			fmri_masked = masker.fit_transform(func_data)                                   # mask the betas according to anatomical rois
			fmri_masked = fmri_masked[conditions_mask]                          # mask according to conditions of interest
			
			# --- classification
			svc = svm.NuSVC(kernel='linear')
			cv = LeaveOneGroupOut()
			
			# rfe + cv
			rfecv = RFECV(estimator=svc, step=1, cv=cv.split(fmri_masked, conditions, groups=loa_label), n_jobs=-1,
			              scoring='accuracy')
			rfecv.fit(fmri_masked, conditions)
			
			#print("Subject: %s, mask: %s, condition: %s, optimal number of features: %d" % (subj, mask, conditions_of_interest, rfecv.n_features_))
			
			conditions_data[conditions_of_interest[0] + '-' + conditions_of_interest[1]] = rfecv.grid_scores_
			
		if plt_enum % 2 == 0:       # even
			subplot_row = plt_enum - row_enum
			subplot_col = 0
		else:                       # odd
			row_enum += 1
			subplot_row = plt_enum - row_enum
			subplot_col = 1
		
		for cond_name, cond_data in zip(conditions_data.keys(), conditions_data.values()):
			msk_plt = ax[subplot_row, subplot_col].plot(range(1, len(cond_data) + 1), cond_data, linewidth=4, label=cond_name)
			
		ax[subplot_row, subplot_col].set_ylim([0.3, 1.0])
		ax[subplot_row, subplot_col].set_title(mask, fontsize=18)
		ax[subplot_row, subplot_col].set_ylabel('CV accuracy')
		ax[subplot_row, subplot_col].set_xlabel('Features')
		
		plt_enum += 1
	
	fig.legend(loc="lower center", borderaxespad=1, fontsize=20,
	           ncol=len([x[0] + '-' + x[1] for x in list(itertools.combinations(conditions_to_decode, 2))]))
	fig.suptitle('RFE with CV for subject %s' % subj, fontsize=24)
	
	plt.savefig(opj(nilearn_dir, subj_str, 'results', subj_str + '_rfe_cv.png'), dpi=300)
	plt.close()
	image_list.append(opj(nilearn_dir, subj_str, 'results', subj_str + '_rfe_cv.png'))
	sbar.set_postfix_str(s='File: %s' % opj('results', subj_str + '_rfe_cv.png'))
	
# create html report
create_html_report(image_list, opj(nilearn_dir, 'group_rfe_cv.html'), 'RFECV results')
