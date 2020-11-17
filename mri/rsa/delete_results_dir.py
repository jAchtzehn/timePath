from os.path import join as opj
from os.path import abspath
from sys import platform
import os


# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/media/sf_data/fMRI/timePath/')

rsa_dir = opj(experiment_dir, 'rsa')

subjects = range(1, 26)
space = 'T1w'

for subj in subjects:
	subj_str = 'sub-' + str(subj).zfill(2)
	
	os.system('rm -rf %s' % opj(rsa_dir, subj_str, 'space-' + space, 'results'))