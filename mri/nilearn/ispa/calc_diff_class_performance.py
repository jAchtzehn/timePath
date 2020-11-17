from nilearn.image import index_img, mean_img, new_img_like, concat_imgs, load_img, clean_img, smooth_img, math_img
from nilearn.input_data import NiftiMasker
from os.path import join as opj
from os.path import abspath
from sys import platform
import pandas as pd
import os
import itertools
from scipy.stats import ttest_ind, wilcoxon
from statsmodels.stats.multitest import multipletests
import numpy as np

# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/mnt/work/achtzehnj/data/')

nilearn_dir = opj(experiment_dir, 'nilearn')
results_dir = opj(nilearn_dir, 'group_results_ispa_std')

# ------------ options ------------
settings = {'splits': 24,
            'conditions_to_decode': ['time-vs-dist', 'time-vs-dots', 'dist-vs-dots']}


for combination in list(itertools.combinations(settings['conditions_to_decode'], 2)):

	imglist_1 = {'conditions': combination[0],
	             'imgs': load_img(opj(results_dir, combination[0], 'ispa_sl_*.nii'))}
	imglist_2 = {'conditions': combination[1],
	             'imgs': load_img(opj(results_dir, combination[1], 'ispa_sl_*.nii'))}

	diff_results_dir = opj(results_dir, combination[0] + '_' + combination[1])
	os.system('mkdir -p %s' % diff_results_dir)

	# create images of differences
	for imgNr in range(0, 24):
		imgDiff = math_img('np.abs((img1 + 0.5) - (img2 + 0.5))', img1=index_img(imglist_1['imgs'], imgNr), img2=index_img(imglist_2['imgs'], imgNr))

		imgDiff.to_filename(opj(diff_results_dir, 'diff_image_' + combination[0] + '_' + combination[1] + '_split-' + str(imgNr).zfill(2) + '.nii'))


	# masker = NiftiMasker(standardize=False, detrend=False)
	# imglist_1['masked_imgs'] = masker.fit_transform(imglist_1['imgs'])
	# imglist_2['masked_imgs'] = masker.fit_transform(imglist_2['imgs'])
	#
	# p_result = np.empty([1, imglist_1['masked_imgs'].shape[1]])
	# t_result = np.empty([1, imglist_1['masked_imgs'].shape[1]])
	# mean_result = np.empty([1, imglist_1['masked_imgs'].shape[1]])
	# for voxel_idx in range(imglist_1['masked_imgs'].shape[1]):
	# 	stats, p = wilcoxon(imglist_1['masked_imgs'][:, voxel_idx], imglist_2['masked_imgs'][:, voxel_idx])
	# 	diff_mean = np.abs(np.mean(imglist_1['masked_imgs'][:, voxel_idx]) - np.mean(imglist_2['masked_imgs'][:, voxel_idx]))
	#
	# 	if p <= 0.05:
	# 		p_result[0, voxel_idx] = p
	# 	else:
	# 		p_result[0, voxel_idx] = 0
	#
	# 	t_result[0, voxel_idx] = stats
	# 	mean_result[0, voxel_idx] = diff_mean
	#
	# # create images
	#
	# p_result_img = masker.inverse_transform(p_result)
	# t_result_img = masker.inverse_transform(t_result)
	# mean_result_img = masker.inverse_transform(mean_result)
	#
	# p_result_img.to_filename(opj(results_dir, combination[0] + '_' + combination[1] + '_p.nii'))
	# t_result_img.to_filename(opj(results_dir, combination[0] + '_' + combination[1] + '_t.nii'))
	# mean_result_img.to_filename(opj(results_dir, combination[0] + '_' + combination[1] + '_mean.nii'))

