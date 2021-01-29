import os
from os.path import join as opj
from os.path import abspath
from sys import platform
import pandas as pd
import itertools
from nilearn.image import load_img, index_img, mean_img
from nilearn.glm import threshold_stats_img
import numpy as np
from nilearn.input_data import NiftiMasker
from nilearn.masking import _extrapolate_out_mask
import scipy.stats as stats
from tqdm import tqdm
import matplotlib.pyplot as plt
from nilearn.plotting import plot_glass_brain
import seaborn as sns
sns.set_theme(style="darkgrid")
from multiprocessing import Pool


# ------------ File I/O ------------
experiment_dir = abspath('/home/achtzehnj/data/timePath/')
# experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath')

rsa_dir = opj(experiment_dir, 'derivatives', 'rsa')
behav_file = opj(experiment_dir, 'derivatives', 'behavioural', 'pse_data_cross_dim_individual_norm.tsv')
behav_data = pd.read_csv(behav_file, delimiter='\t')
mask_img = load_img(opj(experiment_dir, 'derivatives', 'rsa', 'group_mask_wb_binarized.nii.gz'))

# ------------ options ------------
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
conditions = ['time', 'dist', 'dots']

masker = NiftiMasker(mask_img=mask_img, standardize=False, detrend=False,
                     memory=abspath('/home/achtzehnj/code/nilearn_utilities/cache'), memory_level=1)

# compute mean consistency map
cons_imgs = []
for subjectID in subjects:
	cons_imgs.append(opj(rsa_dir, 'sub-' + str(subjectID).zfill(2),
	                     'space-MNI152NLin2009cAsym', 'results',
	                     'sub-' + str(subjectID).zfill(
		                     2) + '_space-MNI152NLin2009cAsym_wb_cons_values.nii.gz'))

cons_mean_img = mean_img(cons_imgs)
cons_mean_img.to_filename(opj(rsa_dir, 'mean_consistency.nii.gz'))

#f, ax = plt.subplots(6, 1, figsize=(6, 20))
#axNr = 0


def calc_crossDim(condition_pair):
	# load crossDim values for all 24 subjects into ndarray
	crossDim = np.abs(behav_data['pse_diff'][(behav_data.condition_rel == condition_pair[0]) &
	                                         (behav_data.condition_irrel == condition_pair[1])].values)

	if (condition_pair[0] == 'time' and condition_pair[1] == 'dist') or (
			condition_pair[1] == 'time' and condition_pair[0] == 'dist'):
		dimId = 0
	elif (condition_pair[0] == 'time' and condition_pair[1] == 'dots') or (
			condition_pair[1] == 'time' and condition_pair[0] == 'dots'):
		dimId = 1
	else:
		dimId = 2

	# now load up each participant's RDM map and correlate
	rdm_dim_data = np.zeros((len(subjects), 64291))


	for s, subjectID in enumerate(subjects):
		rdm_map = load_img(opj(rsa_dir, 'sub-' + str(subjectID).zfill(2), 'space-MNI152NLin2009cAsym', 'results',
		                       'sub-' + str(subjectID).zfill(2) + '_space-MNI152NLin2009cAsym_wb_rdm_values.nii.gz'), )

		# extract right 4D volume from nifti file (0 = time vs. space, 1 = time vs. num, 2 = space vs. num
		rdm_map_dim = index_img(rdm_map, dimId)

		rdm_dim_data[s, :] = masker.fit_transform(rdm_map_dim)

	corr_data = np.zeros((1, 64291))
	corr_data_p = np.zeros((1, 64291))

	# load mean consistency map
	cons_data = masker.fit_transform(cons_mean_img)
	cons_data_norm = cons_data / np.max(cons_data)

	vbar = tqdm(range(rdm_dim_data.shape[1]), leave=False, mininterval=1)
	for vx in vbar:
		[vx_corr_r, vx_corr_p] = stats.spearmanr(rdm_dim_data[:, vx], crossDim.T)
		corr_data[0, vx] = vx_corr_r
		if vx_corr_p <= 0.05 and cons_data[0, vx] > 0.3 and vx_corr_r < 0:
			corr_data_p[0, vx] = vx_corr_r
	corr_data = corr_data * np.abs(cons_data_norm)

	# sns.regplot(x=rdm_dim_data[:, corr_data_p.argmin()], y=crossDim.T, ax=ax[axNr])
	# ax[axNr].set_xlabel('RDM')
	# ax[axNr].set_ylabel('CrossDim')
	# ax[axNr].set_title('rel-{}_irrel-{}'.format(condition_pair[0], condition_pair[1]))
	# axNr = axNr + 1

	# zscore data
	corr_data_z = (corr_data - np.mean(corr_data)) / np.std(corr_data)

	corr_img = masker.inverse_transform(corr_data)
	corr_img_p = masker.inverse_transform(corr_data_p)
	corr_img_z = masker.inverse_transform(corr_data_z)

	[corr_img_thresh, th] = threshold_stats_img(corr_img_z, cluster_threshold=5, alpha=0.05, height_control='fdr')

	plot_glass_brain(corr_img_p, colorbar=True, cmap='viridis', plot_abs=False, threshold=0, vmax=1, symmetric_cbar=False)
	plt.savefig(opj(rsa_dir, 'corrImg_rel-{}_irrel-{}_p.pdf'.format(condition_pair[0], condition_pair[1])), format='pdf')
	plt.close()

	# corr_img_thresh.to_filename(opj(rsa_dir, 'corrImg_rel-{}_irrel-{}_thresh.nii.gz'.format(condition_pair[0], condition_pair[1])))
	corr_img.to_filename(opj(rsa_dir, 'corrImg_rel-{}_irrel-{}.nii.gz'.format(condition_pair[0], condition_pair[1])))
	corr_img_p.to_filename(
		opj(rsa_dir, 'corrImg_rel-{}_irrel-{}_p.nii.gz'.format(condition_pair[0], condition_pair[1])))

pool = Pool(6)
pool.map(calc_crossDim, itertools.permutations(conditions, 2))