import os
from os.path import join as opj
from os.path import abspath
from sys import platform
from util_functions import createSubjectlist
import pandas as pd

# ------------ options ------------
subjects = range(1, 2)
output_space = 'MNI152NLin2009cAsym'           # T1w or MNI152NLin2009cAsym are possible
ses_list = range(1, 3)                  # sessions that should be included
run_list = range(1, 5)                  # runs that should be included
trial_list = range(3, 66)               # trials that should be included (to get equal label numbers, remove the first two trials from each run)

# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/mnt/work/achtzehnj/data/')

nipype_dir = opj('output_nipype', 'tw_beta_estimation')  # root folder of beta images
results_dir = 'nilearn'  # name of the results dir

subj_list = createSubjectlist(subjects)  # which subjects should be analysed

# create nilearn data folder
if results_dir not in os.listdir(opj(experiment_dir)):
	os.system('mkdir %s' % opj(experiment_dir, results_dir))
	print('Created folder: %s' % opj(experiment_dir, results_dir))

# ------------ 1. get the first beta image of every trial and create a single .nii image ------------

for subj in subj_list:
	# create subj folder
	os.system('mkdir -p %s' % opj(experiment_dir, results_dir, subj, 'space-' + output_space, 'betas'))
	
	# create base command to merge betas into on .nii.gz file
	filename_merged_betas = opj(experiment_dir, results_dir, subj, 'space-' + output_space, 'betas', subj + '_betas_merged_warning.nii')
	merge_cmd = 'fslmerge -t ' + filename_merged_betas
	
	print('Copying for subject %s...' % subj)
	for ses in ses_list:
		for run in run_list:
			for trial in trial_list:
				
				# get trial info (to exclude thumb responses, i.e. the participant has forgotten the trial)
				trial_info = pd.read_csv(opj(experiment_dir, nipype_dir, subj, 'space-' + output_space, 'ses-' + str(ses), ('run-' + str(run)), ('trial-' + str(trial)), 'trial_info.tsv'), delimiter='\t')
				
				# remove for every odd run the last trial and for every even run the last two trials
				if run % 2 == 0 and trial == 65:
					pass
				else:
					if not trial_info['trial_type'][0] == 'thumb':
						os.system('cp %s %s' % (opj(experiment_dir, nipype_dir, subj, 'space-' + output_space, 'ses-' + str(ses), ('run-' + str(run)), ('trial-' + str(trial)), 'beta_0001.nii'),
						opj(experiment_dir, results_dir, subj, 'space-' + output_space, 'betas', (subj + '_ses-' + str(ses).zfill(2) + '_run-' + str(run).zfill(2) + '_trial-' + str(trial).zfill(2) + '_beta.nii'))))
						
						# expand merge cmd
						merge_cmd = merge_cmd + ' ' + opj(experiment_dir, results_dir, subj, 'space-' + output_space, 'betas', (subj + '_ses-' + str(ses).zfill(2) + '_run-' + str(run).zfill(2) + '_trial-' + str(trial).zfill(2) + '_beta.nii'))
	
	# merge files
	print('Merging for subject %s...' % subj)
	os.system(merge_cmd)
	
	# delete individual beta files
	for niiFile in os.listdir(opj(experiment_dir, results_dir, subj, 'space-' + output_space, 'betas')):
		if '_beta.nii' in niiFile:
			os.system('rm -rf %s'%opj(experiment_dir, results_dir, subj, 'space-' + output_space, 'betas', niiFile))