from os.path import join as opj
from os.path import abspath
from sys import platform
import os
import pandas as pd
import nilearn as nl

# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/mnt/work/achtzehnj/data/')

nilearn_dir = opj(experiment_dir, 'nilearn')
rsa_dir = opj(experiment_dir, 'rsa')

os.system('mkdir -p %s' % rsa_dir)  # create the mvpa2 folder for data storage

subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]

combinations = ['dots-dist', 'dist-time', 'dots-time']


for combination in combinations:

	file_list = []
	for subj in subjects:
		file_list.append(opj(rsa_dir, 'sub-' + str(subj).zfill(2), 'space-MNI152NLin2009cAsym',
		                     'results',
		                     'sub-' + str(subj).zfill(2) + '_space-MNI152NLin2009cAsym_wb_similarity_' + combination + '_searchlight.nii.gz'))

	mean_img = nl.image.mean_img(file_list)
	mean_img.to_filename(opj(rsa_dir, 'group_results_space-MNI152NLin2009cAsym', 'rsa_' + combination + '.nii.gz'))