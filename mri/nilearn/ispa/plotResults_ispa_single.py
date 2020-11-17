from os.path import join as opj
from os.path import abspath
from nilearn.plotting import plot_glass_brain, plot_surf_stat_map
from nilearn.image import mean_img
from nilearn import datasets
from nilearn import surface
from scipy import ndimage
from sys import platform
import matplotlib.pyplot as plt
import matplotlib as mpl
import pylab as pl

from nilearn.image import load_img, new_img_like, threshold_img
from nilearn._utils.niimg_conversions import _safe_get_data
import numpy as np

if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/mnt/work/achtzehnj/data')

nilearn_dir = opj(experiment_dir, 'nilearn')
results_dir = opj(nilearn_dir, 'group_results_ispa_std')
thesis_dir = abspath('/Users/jachtzehn/Documents/Medizin/thesis/figures/results/fmri')

# ---- options ----
colormap = 'viridis'
output_name = '_group_level_mvpa_masked'
output_format = 'png'
dpi = 300
conditions_to_decode = [['time', 'dist'], ['time', 'dots'], ['dist', 'dots'], ['time', 'lumin'], ['dist', 'lumin'], ['dots', 'lumin']]
mask_img = True     # mask image with cluster mask?
plot = False
nrSplits = 24

for condition_pair in conditions_to_decode:

	if len(condition_pair) > 1:
		condition_dir = opj(results_dir, condition_pair[0] + '-vs-' + condition_pair[1], 'snpm_batch')
		output_fname = condition_pair[0] + '-vs-' + condition_pair[1] + output_name
		split_dir = opj(results_dir, condition_pair[0] + '-vs-' + condition_pair[1])

	else:
		condition_dir = opj(results_dir, condition_pair[0], 'snpm_batch')
		split_dir = opj(results_dir, condition_pair[0])
		output_fname = condition_pair[0] + output_name

	img_fname = opj(condition_dir, 'SnPM_filtered.nii')
	mask_fname = opj(condition_dir, 'SnPM_filtered_binary_mask_cluster.nii')

	img = load_img(img_fname)
	img_acc = mean_img(opj(split_dir, 'ispa_sl_wb_decoding-' + condition_pair[0] + '-vs-' + condition_pair[1] + '_split-*.nii'))

	mask = load_img(mask_fname)

	# 1. plot significant cluster
	fig = plt.figure(figsize=(15, 15), constrained_layout=False)
	gs = fig.add_gridspec(3, 2)
	ax_glass = fig.add_subplot(gs[0, :])  # top row
	ax_surf_top_l = fig.add_subplot(gs[1, 0], projection='3d')
	ax_surf_top_r = fig.add_subplot(gs[1, 1], projection='3d')
	ax_surf_bottom_l = fig.add_subplot(gs[2, 0], projection='3d')
	ax_surf_bottom_r = fig.add_subplot(gs[2, 1], projection='3d')

	if mask_img:
		try:
			img = threshold_img(img, 0, mask_img=mask)
			img_acc = threshold_img(img_acc, 0, mask_img=mask)
			# write out masked image
			img.to_filename(opj(condition_dir, 'SnPM_filtered_masked.nii'))
			img_acc.to_filename(opj(condition_dir, 'SnPM_filtered_masked_accuracy.nii'))
		except ValueError:
			img = new_img_like(img, np.zeros(shape=img.shape))
			print('Mask was empty, not masking')
	else:
		img = threshold_img(img, 0)

	if plot:
		# find maximum for plot thresholding
		img_data = _safe_get_data(img, ensure_finite=True)
		max_val = ndimage.maximum(
			img_data)  # multiply this value with the percent threshold to get absolute value of threshold
		if max_val > 0:
			min_val = ndimage.minimum(img_data[img_data > 0])
		else:
			min_val = max_val

		plot_glass_brain(img, display_mode='lyrz', colorbar=False, plot_abs=False,
		                 vmin=min_val, vmax=max_val, cmap=colormap, figure=fig, axes=ax_glass)

		# 3D plot
		fsaverage = datasets.fetch_surf_fsaverage('fsaverage')
		# plot in one fig

		# textures and bg
		bg_map_l = fsaverage.sulc_left
		mesh_l = fsaverage.infl_left
		texture_l = surface.vol_to_surf(img, fsaverage.pial_left)
		bg_map_r = fsaverage.sulc_right
		mesh_r = fsaverage.infl_right
		texture_r = surface.vol_to_surf(img, fsaverage.pial_right)

		# top
		plot_surf_stat_map(mesh_l, texture_l, hemi='left', colorbar=False, view='lateral', cmap=colormap,
		                   bg_map=bg_map_l, threshold=0.001, figure=fig, axes=ax_surf_top_l)
		plot_surf_stat_map(mesh_r, texture_r, hemi='right', colorbar=False, view='lateral', cmap=colormap,
		                   bg_map=bg_map_r, threshold=0.001, figure=fig, axes=ax_surf_top_r)
		# bottom
		plot_surf_stat_map(mesh_l, texture_l, hemi='left', colorbar=False, view='medial', cmap=colormap,
		                   bg_map=bg_map_l, threshold=0.001, figure=fig, axes=ax_surf_bottom_l)
		plot_surf_stat_map(mesh_r, texture_r, hemi='right', colorbar=False, view='medial', cmap=colormap,
		                   bg_map=bg_map_r, threshold=0.001, figure=fig, axes=ax_surf_bottom_r)

		# cbar_ax2 = fig_sp.add_axes([0.25, 0.5, 0.5, 0.05])
		# cb2 = mpl.colorbar.ColorbarBase(cbar_ax2, cmap=cmap_orig,
		#                                 norm=norm,
		#                                 orientation='horizontal')

		# construct cmap
		norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
		cmap_orig = mpl.cm.get_cmap(colormap)

		# fig.subplots_adjust(top=0.9)
		cbar_ax = fig.add_axes([0.15, 0.925, 0.7, 0.05])
		cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap_orig,
		                                norm=norm,
		                                orientation='horizontal')

		plt.tight_layout()
		plt.savefig(opj(results_dir, output_fname + '.' + output_format), fig=fig, dpi=dpi, format=output_format)
		print('Written image: ' + str(opj(results_dir, output_fname + '.' + output_format)))
		plt.close()

		# 2. plot average classification accuracy
		if len(condition_pair) > 1:
			output_fname = condition_pair[0] + '-vs-' + condition_pair[1] + output_name + '_accuracies'
		else:
			output_fname = condition_pair[0] + output_name + '_accuracies'

		img = mean_img(opj(split_dir, 'ispa_sl_wb_decoding-' + condition_pair[0] + '-vs-' + condition_pair[1] + '_split-*.nii'))

		# 1. plot significant cluster
		fig = plt.figure(figsize=(15, 15), constrained_layout=False)
		gs = fig.add_gridspec(3, 2)
		ax_glass = fig.add_subplot(gs[0, :])  # top row
		ax_surf_top_l = fig.add_subplot(gs[1, 0], projection='3d')
		ax_surf_top_r = fig.add_subplot(gs[1, 1], projection='3d')
		ax_surf_bottom_l = fig.add_subplot(gs[2, 0], projection='3d')
		ax_surf_bottom_r = fig.add_subplot(gs[2, 1], projection='3d')

		if mask_img:
			try:
				img = threshold_img(img, 0, mask_img=mask)
			except ValueError:
				img = new_img_like(img, np.zeros(shape=img.shape))
				print('Mask was empty, not masking')
		else:
			img = threshold_img(img, 0)

		# find maximum for plot thresholding
		img_data = _safe_get_data(img, ensure_finite=True)

		max_val = ndimage.maximum(
			img_data)  # multiply this value with the percent threshold to get absolute value of threshold
		if max_val > 0:
			min_val = ndimage.minimum(img_data[img_data > 0])
		else:
			min_val = max_val

		plot_glass_brain(img, display_mode='lyrz', colorbar=False, plot_abs=False,
		                 vmin=min_val, vmax=max_val, cmap=colormap, figure=fig, axes=ax_glass)

		# 3D plot
		fsaverage = datasets.fetch_surf_fsaverage('fsaverage')

		# textures and bg
		bg_map_l = fsaverage.sulc_left
		mesh_l = fsaverage.infl_left
		texture_l = surface.vol_to_surf(img, fsaverage.pial_left)
		bg_map_r = fsaverage.sulc_right
		mesh_r = fsaverage.infl_right
		texture_r = surface.vol_to_surf(img, fsaverage.pial_right)

		# top
		plot_surf_stat_map(mesh_l, texture_l, hemi='left', colorbar=False, view='lateral', cmap=colormap,
		                   bg_map=bg_map_l, threshold=0.001, figure=fig, axes=ax_surf_top_l)
		plot_surf_stat_map(mesh_r, texture_r, hemi='right', colorbar=False, view='lateral', cmap=colormap,
		                   bg_map=bg_map_r, threshold=0.001, figure=fig, axes=ax_surf_top_r)
		# bottom
		plot_surf_stat_map(mesh_l, texture_l, hemi='left', colorbar=False, view='medial', cmap=colormap,
		                   bg_map=bg_map_l, threshold=0.001, figure=fig, axes=ax_surf_bottom_l)
		plot_surf_stat_map(mesh_r, texture_r, hemi='right', colorbar=False, view='medial', cmap=colormap,
		                   bg_map=bg_map_r, threshold=0.001, figure=fig, axes=ax_surf_bottom_r)

		# construct cmap
		norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
		cmap_orig = mpl.cm.get_cmap(colormap)

		# fig.subplots_adjust(top=0.9)
		cbar_ax = fig.add_axes([0.15, 0.925, 0.7, 0.05])
		cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap_orig,
		                                norm=norm,
		                                orientation='horizontal')

		plt.tight_layout()
		plt.savefig(opj(results_dir, output_fname + '.' + output_format), fig=fig, dpi=dpi, format=output_format)
		print('Written image: ' + str(opj(results_dir, output_fname + '.' + output_format)))
		plt.close()