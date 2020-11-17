"""
Decodes various labels (defined in options) of the timePath 2 data
"""

from sklearn import svm
from sklearn.model_selection import LeaveOneGroupOut
from nilearn.image import index_img, mean_img, new_img_like, clean_img
from sklearn.linear_model import LogisticRegression
from nilearn.decoding import SearchLight
from os.path import join as opj
from os.path import abspath
from sys import platform
import pandas as pd
import os
from tqdm import tqdm

# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/mnt/work/achtzehnj/data/')

nilearn_dir = opj(experiment_dir, 'nilearn')

# ------------ options ------------
settings = {'subjects': range(1, 26),
            'conditions_to_decode': [['time', 'dist'], ['time', 'dots'], ['dist', 'dots'], ['time', 'lumin'], ['dist', 'lumin'], ['dots', 'lumin']],
            'n_proc': -1,
            'roi': 'wb',
            'standardize': True,
            'solver': 'lbfgs',
            'verbose_decode': 0,
            'file_suffix': 'score-balanced-accuracy_' + 'solver-lbfgs'}

# ------------ code ------------
# print settings in console
for key in settings.keys():
	print('Option: {}, value: {}'.format(key, settings[key]))
print('\n')

score_data = {}

sbar = tqdm(settings['subjects'], leave=True)
for subj in sbar:

	sbar.set_description('Processing subject %s' % str(subj))

	subj_data = {}     # dict in which data is stored for subj
	subj_str = 'sub-' + str(subj).zfill(2)      # create the subject string (just for easier handling)
	result_dir = opj(nilearn_dir, subj_str, 'space-' + 'MNI152NLin2009cAsym', 'results')
	os.system('mkdir -p %s' % result_dir)        # make sure the results dir exists

	# make results dir
	if not os.path.exists(opj(nilearn_dir, subj_str, 'space-' + 'MNI152NLin2009cAsym', 'results')):
		os.mkdir(opj(nilearn_dir, subj_str, 'space-' + 'MNI152NLin2009cAsym', 'results'))

	anat_data_name = opj(experiment_dir, 'fmriprep', subj_str, 'anat', subj_str + '_T1w_preproc.nii.gz')    # T1 image
	func_data_name = opj(nilearn_dir, subj_str, 'space-' + 'MNI152NLin2009cAsym', 'betas', subj_str + '_betas_merged.nii.gz')                 # betas
	mask_filename = opj(nilearn_dir, subj_str, 'space-' + 'MNI152NLin2009cAsym', 'masks', subj_str + '_' + 'wb' + '_mask_binarized.nii.gz')
	behav_filename = opj(nilearn_dir, subj_str, 'space-' + 'MNI152NLin2009cAsym', 'betas', subj_str + '_merged_events.tsv')                   # events file

	func_data = clean_img(func_data_name, standardize=settings['standardize'], ensure_finite=True, detrend=False, mask_img=None)

	# load behavioral information
	behavioral = pd.read_csv(behav_filename, delimiter='\t')        # read in the data of the file
	conditions_all = behavioral['trial']                            # these are the conditions, i.e. trial types

	conditions_data = {}                    # dict in which data is stored for conditions of interest
	cbar = tqdm(settings['conditions_to_decode'], leave=False)
	for conditions_of_interest in cbar:

		cbar.set_description('Processing conditions %s vs. %s' %(conditions_of_interest[0], conditions_of_interest[1]))

		conditions_mask = conditions_all.isin(conditions_of_interest)       # mask to restrict the conditions to the ones specified by options
		loa_label = behavioral['class_chunk'][conditions_mask]              # mask that specifies the runs masked by the conditions

		# apply mask to fMRI data and conditions
		fmri_masked = index_img(func_data, conditions_mask)            # mask the betas according to labels
		conditions = conditions_all[conditions_mask]                        # mask behavioral data
		mean_fmri = mean_img(fmri_masked)

		# get the number of samples
		n_samples = [conditions.tolist().count(conditions_of_interest[0]), conditions.tolist().count(conditions_of_interest[1])]

		# setup searchlight
		svc = svm.NuSVC(kernel='linear')
		clf = LogisticRegression(C=0.1, penalty='l2', solver=settings['solver'])
		searchlight = SearchLight(
			mask_img=mask_filename,
			radius=6, n_jobs=settings['n_proc'],
			verbose=0, cv=LeaveOneGroupOut(),
			estimator=clf, scoring='balanced_accuracy')

		# fit data of non-randomised data
		searchlight.fit(fmri_masked, conditions, groups=loa_label)

		# save .nii.gz for non-randomized data
		img_name = subj_str + '_sl-whole-brain_' + conditions_of_interest[0] + '_vs_' + conditions_of_interest[1] + '.nii.gz'
		searchlight_img = new_img_like(mean_fmri, searchlight.scores_)
		searchlight_img.to_filename(opj(result_dir, img_name))