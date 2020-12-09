"""
Calculates TSNR of different ROIs (defined by <masks>) and compares them

"""


from os.path import join as opj
from os.path import abspath
from sys import platform
from tqdm import tqdm
import matplotlib as mpl
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import numpy as np
from nilearn.masking import apply_mask
import _pickle as cPickle

# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/media/sf_data/fMRI/timePath/')

nilearn_dir = opj(experiment_dir, 'nilearn')

# ------------ options ------------
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
# [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
masks = ['ips_lh', 'ips_rh', 'hc_lh', 'hc_rh', 'ifg_lh', 'ifg_rh', 'pcg_lh', 'pcg_rh', 'V5_rh', 'V5_lh', 'insula_rh', 'insula_lh']
# ['ips_lh', 'ips_rh', 'hc_lh', 'hc_rh', 'ifg_lh', 'ifg_rh', 'pcg_lh', 'pcg_rh', 'V5_rh', 'V5_lh', 'insula_rh', 'insula_lh']
space = 'MNI152NLin2009cAsym'
input_data = 'epi'            # betas or epi
calc_tsnr = True
plot_results = True


if calc_tsnr:
	# ------------ run ----------------
	tsnr_data = {}      # dict in which to store the data of each ROI (key) as lists
	for mask in masks:
		tsnr_data[mask] = []

	print('Input data: %s' % input_data)

	sbar = tqdm(subjects, leave=True)
	for subj in sbar:
		sbar.set_description('Processing subject %s' % str(subj))
		subj_str = 'sub-' + str(subj).zfill(2)      # create the subject string (just for easier handling)

		mbar = tqdm(masks, leave=False)
		for mask in mbar:
			mbar.set_description('Processing mask %s' % str(mask))

			if input_data == 'betas':
				func_data_name = opj(nilearn_dir, subj_str, 'space-' + space, 'betas', subj_str + '_betas_merged.nii.gz')                   # betas
				if mask == 'V5_rh' or mask == 'V5_lh':
					mask_filename = opj(nilearn_dir, 'group_masks', 'space-' + space,
					                    'group_mask_' + mask + '_binarized.nii.gz')  # mask filename
				else:
					mask_filename = opj(nilearn_dir, subj_str, 'space-' + space, 'masks',
					                    subj_str + '_' + mask + '_mask_binarized.nii.gz')  # mask filename

				masked_data = apply_mask(func_data_name, mask_filename)     # this will give us an array of the size [timepoints, voxels]
				masked_data = masked_data.std(axis=0).nonzero()[0]     # this will give us an array of the size [timepoints, voxels]
				masked_data = np.all(np.isfinite(masked_data), axis=0)     # this will give us an array of the size [timepoints, voxels]

				mean_masked_data = np.mean(masked_data, axis=0)             # mean along time axis
				std_masked_data = np.std(masked_data, axis=0)               # std along time axis

				# sanity check to see if any voxel is 0 for all volumes
				idx_zero = np.where(std_masked_data == 0)[0]
				# if any zeros where found, take the mean of the adjacent voxels as values for this voxel
				if idx_zero.size != 0:
					print('Found zeros!')
					# for idx in idx_zero:
					# 	masked_data[:, idx] = (masked_data[:, idx - 1] + masked_data[:, idx + 1]) / 2
					#
					# # recompute mean and std
					# mean_masked_data = np.mean(masked_data, axis=0)             # mean along time axis
					# std_masked_data = np.std(masked_data, axis=0)               # std along time axis

				tsnr_masked_data = mean_masked_data / std_masked_data       # calculate tSNR by: mean/std of each voxel
				tsnr_masked_data_mean = np.mean(tsnr_masked_data)           # calculate the mean over all voxels to get tSNR of ROI

				tsnr_data[mask].append(tsnr_masked_data_mean)

			# if epis are to be used, the runs have to be calculated individually
			elif input_data == 'epi':
				# create list with functionl and mask files
				tsnr_masked_data_mean_runs = []

				for ses in range(1, 3):
					for run in range(1, 5):
						func_data_name = opj(experiment_dir, 'fmriprep', subj_str, 'ses-' + str(ses).zfill(2), 'func',
						                     subj_str + '_ses-' + str(ses).zfill(2) + '_task-class_run-' +
						                     str(run).zfill(2) + '_bold_space-' + space + '_preproc.nii.gz')

						if mask == 'V5_rh' or mask == 'V5_lh':
							mask_filename = opj(nilearn_dir, 'group_masks', 'space-' + space,
							                    'group_mask_' + mask + '_binarized.nii.gz')  # mask filename
						else:
							mask_filename = opj(nilearn_dir, subj_str, 'space-' + space, 'masks',
							                    subj_str + '_' + mask + '_mask_binarized.nii.gz')  # mask filename

						masked_data = apply_mask(func_data_name, mask_filename, ensure_finite=True)     # this will give us an array of the size [timepoints, voxels]
						masked_data = masked_data[:, masked_data.std(axis=0).nonzero()[0]]  # remove invariant features (voxels)

						mean_masked_data = np.mean(masked_data, axis=0)             # mean along time axis
						std_masked_data = np.std(masked_data, axis=0)               # std along time axis

						# sanity check to see if any voxel is 0 for all volumes
						idx_zero = np.where(std_masked_data == 0)[0]
						# if any zeros where found, take the mean of the adjacent voxels as values for this voxel
						# if idx_zero.size != 0:
						# 	for idx in idx_zero:
						# 		masked_data[:, idx] = (masked_data[:, idx - 1] + masked_data[:, idx + 1]) / 2
						#
						# 	# recompute mean and std
						# 	mean_masked_data = np.mean(masked_data, axis=0)             # mean along time axis
						# 	std_masked_data = np.std(masked_data, axis=0)               # std along time axis

						tsnr_masked_data = mean_masked_data / std_masked_data       # calculate tSNR by: mean/std of each voxel
						tsnr_masked_data_mean_runs.append(np.mean(tsnr_masked_data))           # calculate the mean over all voxels to get tSNR of ROI

				tsnr_data[mask].append(np.mean(tsnr_masked_data_mean_runs))

	with open(opj(experiment_dir, 'data_quality_plots', 'tSNR_data.pkl'), 'wb+') as f:
		cPickle.dump(tsnr_data, f)
		f.close()

if plot_results:

	with open(opj(experiment_dir, 'data_quality_plots', 'tSNR_data.pkl'), 'rb') as f:
		tsnr_data = cPickle.load(f)
		f.close()

	# convert dict to nested array
	boxData = []
	for roi in tsnr_data.keys():
		boxData.append(tsnr_data[roi])

	fig, ax1 = plt.subplots(figsize=(15, 8))
	bp = ax1.boxplot(boxData, notch=False, vert=True, whis=1.5)

	plt.tight_layout(5)
	numBoxes = len(boxData)
	medians = list(range(numBoxes))

	# lim
	maxVal = np.max(boxData)
	minVal = np.min(boxData)
	top = maxVal + 0.1 * abs(maxVal)
	bottom = minVal - 0.2 * abs(minVal)
	ax1.set_ylim(bottom, top)

	plt.setp(bp['boxes'], color='black')
	plt.setp(bp['medians'], color='#2ecc71')
	plt.setp(bp['whiskers'], color='black')
	plt.setp(bp['fliers'], color='red', marker='o')\

	# grid
	ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
				   alpha=0.5)
	ax1.set_axisbelow(True)

	# box color and median values
	cmap = mpl.cm.get_cmap('plasma')
	normalize = mpl.colors.Normalize(vmin=minVal, vmax=maxVal)
	colors = [cmap(normalize(value)) for value in [minVal, maxVal]]

	for i in range(numBoxes):
		box = bp['boxes'][i]
		boxX = []
		boxY = []
		for j in range(5):
			boxX.append(box.get_xdata()[j])
			boxY.append(box.get_ydata()[j])
		boxCoords = list(zip(boxX, boxY))

		# get median
		med = bp['medians'][i]
		medians[i] = med.get_ydata()[0]

		boxPolygon = Polygon(boxCoords, facecolor=cmap(normalize(med.get_ydata()[0])))
		ax1.add_patch(boxPolygon)

	# label
	ax1.set_xlabel('ROIs')
	ax1.set_ylabel('TSNR')
	ax1.set_xticklabels(tsnr_data.keys())

	# xtick above
	pos = np.arange(numBoxes) + 1
	upperLabels = [str(np.round(s, 2)) for s in medians]
	weights = ['bold', 'semibold']
	for tick, label in zip(range(numBoxes), ax1.get_xticklabels()):
		k = tick % 2
		ax1.text(pos[tick], top, upperLabels[tick],
			 horizontalalignment='center', fontsize=12, weight=weights[k],
		         color=cmap(normalize(medians[tick])))

	plt.savefig(opj(experiment_dir, 'data_quality_plots', 'tsnr_both_hemispheres_roi_' + input_data + '.png'))
