"""
This file includes various utility functions that are needed during
the trialwise first level analysis nipype workflow

(c) Johannes Achtzehn, 06.05.2019
"""

def createSubjectlist(subjects):

	subjList = []
	
	for subj in subjects:
		subjList.append('sub-' + str(subj).zfill(2))
	
	return subjList
				
def cleanWorkingOutputDir(output_space, subj_id, ses_id, run_id, trial_id, working_dir, datasink_out):
	
	import os
	from os.path import join as opj
	from time import sleep
	import glob
	
	# working dir
	# construct correct working dir folder and delete subdirs (nodedirs)
	subdir = '_output_space_%s_run_id_%s_ses_id_%s_subj_id_%s_trial_id_%s' % (str(output_space), str(run_id), str(ses_id), str(subj_id), str(trial_id))
	
	node_dirs = os.listdir(opj(working_dir, subdir))
	sleep(5)
	for node_dir in node_dirs:
		if 'n_clean' not in node_dir:
			os.system('rm -rf %s' % opj(working_dir, subdir, node_dir))
			
	# output dir, delete additional beta files that are not needed to save space
	datasink_out = datasink_out[0]
	output_dir = datasink_out[:-datasink_out[::-1].find('/')]   # remove the last file name from the given path (points to a specific file rather than the folder)
	
	betafiles = glob.glob(output_dir + 'beta*')
	
	for i in range(2, len(betafiles) + 1):
		os.system('rm -rf %s' % opj(output_dir, 'beta_' + str(i).zfill(4) + '.nii'))
	
	return []


def plotDesignMatrix(matFile):
	import numpy as np
	from matplotlib import pyplot as plt
	from scipy.io import loadmat
	from os.path import join as opj
	import os
	
	fig_filename_full = opj(os.getcwd(), 'spm_mat.png')      # create filename
	
	spmmat = loadmat(matFile,                           # Using scipy's loadmat function we can access SPM.mat
                 struct_as_record=False)
	
	# normalize data
	designMatrix = spmmat['SPM'][0][0].xX[0][0].X
	names = [i[0] for i in spmmat['SPM'][0][0].xX[0][0].name[0]]
	
	normed_design = designMatrix / np.abs(designMatrix).max(axis=0)
	
	# save figure
	fig_f, ax_f = plt.subplots(figsize=(20, 20))
	plt.imshow(normed_design, aspect='auto', cmap='gray', interpolation='none')
	ax_f.set_ylabel('Volume id')
	ax_f.set_xticks(np.arange(len(names)))
	ax_f.set_xticklabels(names, rotation=90)
	fig_f.savefig(opj(fig_filename_full))
	plt.close(fig_f)
	
	return fig_filename_full

def plotTrialInfo(trial_info, regr_info):
	import pandas
	from os.path import join as opj
	import os
	
	fname_trial = opj(os.getcwd(), 'trial_info.tsv')
	trial_info.to_csv(fname_trial, sep='\t', index=False)
	
	fname_regr = opj(os.getcwd(), 'regr_info.tsv')
	regr_info.to_csv(fname_regr, sep='\t', index=False)
	
	return fname_trial, fname_regr

def getVolNr(behav_dir, subj_id, ses_id, run_id):
	import pandas as pd
	from os.path import join as opj
	
	if ses_id == 2:
		run_id = run_id + 4
	else:
		pass
	
	logfile = pd.read_csv(opj(behav_dir, 'tp2_stana_Nscans.csv'), delimiter=';')    # read in file
	idx_runs = logfile[logfile['run'] == run_id].index.tolist()                     # find all indices of the current run
	idx_subj = idx_runs[int(subj_id[4:]) - 1]                                       # find the idx for the current subject
	nrVol = logfile['volumes'].iloc[idx_subj]                                       # 0 based number of volumes
	
	return nrVol