from os.path import join as opj
from os.path import abspath
from sys import platform
import os
from scipy import ndimage
from nilearn import plotting, image
from nilearn.image import load_img, index_img, iter_img
from nilearn._utils.niimg_conversions import _safe_get_data as get_data
import matplotlib.pyplot as plt


# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/media/sf_data/fMRI/timePath/')

nilearn_dir = opj(experiment_dir, 'nilearn')


subjects = [4]


for subj in subjects:
	
	subj_str = 'sub-' + str(subj).zfill(2)      # create the subject string (just for easier handling)
	subj_folder = opj(nilearn_dir, subj_str, 'betas')
	func_fname = opj(subj_folder, subj_str + '_betas_merged.nii.gz')
	
	# load images
	func_img = load_img(func_fname)
	
	mean_img_vals = []
	for img in iter_img(func_img):
		mean_img_val = ndimage.mean(get_data(img, ensure_finite=True))
		mean_img_vals.append(mean_img_val)

	fig, ax = plt.subplots()
	ax.plot(mean_img_vals)
	plt.show()
