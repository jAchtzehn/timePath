from os.path import join as opj
from os.path import abspath
from nilearn.plotting import plot_glass_brain, plot_surf_stat_map, plot_stat_map
from nilearn import datasets
from nilearn import surface
from scipy import ndimage
from sys import platform
import matplotlib.pyplot as plt
import matplotlib as mpl
import pylab as pl
from mpl_toolkits.mplot3d import axes3d
import argparse
from nilearn.glm import threshold_stats_img
from multiprocessing import Pool

from nilearn.image import load_img, threshold_img
from nilearn._utils.niimg_conversions import _safe_get_data
import numpy as np

if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/home/achtzehnj/data/timePath/')

output_dir = opj(experiment_dir, 'derivatives', 'output_nipype', '2ndlevel')
thesis_dir = abspath('/home/Documents/thesis/figures/results/fmri')
thesis_dir_appendix = abspath('/home/Documents/thesis/figures/appendix')
decode_dir = opj(experiment_dir, 'derivatives', 'nilearn', 'group_results_ispa_std')

colormap = 'viridis'
output_format = 'pdf'
dpi = 400
threshold = 6.39
fs_average = 'fsaverage'

def get_args():
	parser = argparse.ArgumentParser(description='Plot fmri stat maps', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('--type', '-t', metavar='TYPE', help="Which type of plot is to be plotted", required=True, nargs=1)

	parser.add_argument('--contrasts', '-c', metavar='CONTRASTS', nargs='+', help='Contrasts to be plotted', required=False)

	return parser.parse_args()

def plotContrastSurf(contrastIds):

	for x, con in enumerate(contrastIds):

		conDir = opj(output_dir, 'space-' + 'MNI152NLin2009cAsym', 'contrast-' + con)

		# init figure with gridspec
		fig = plt.figure(figsize=(14, 8))
		gs = fig.add_gridspec(2, 4)
		gs.update(wspace=-0.5, hspace=-0.2)  # set the spacing between axes.

		ax_surf_top_l = fig.add_subplot(gs[0, 0], projection='3d')
		ax_surf_top_ml = fig.add_subplot(gs[0, 1], projection='3d')
		ax_surf_top_mr = fig.add_subplot(gs[0, 2], projection='3d')
		ax_surf_top_r = fig.add_subplot(gs[0, 3], projection='3d')
		ax_surf_bottom_l = fig.add_subplot(gs[1, 0], projection='3d')
		ax_surf_bottom_ml = fig.add_subplot(gs[1, 1], projection='3d')
		ax_surf_bottom_mr = fig.add_subplot(gs[1, 2], projection='3d')
		ax_surf_bottom_r = fig.add_subplot(gs[1, 3], projection='3d')

		for ax in fig.get_axes():
			ax.patch.set_alpha(0)

		img = load_img(opj(conDir, 'spmT_0001_thr.nii'))
		img_data = _safe_get_data(img, ensure_finite=True)

		max_val = ndimage.maximum(
			img_data)  # multiply this value with the percent threshold to get absolute value of threshold
		if max_val > 0:
			min_val = 0
		else:
			min_val = max_val

		# 3D plots
		fsaverage = datasets.fetch_surf_fsaverage(fs_average, data_dir='/home/achtzehnj/code/nilearn_utilities/atlases/')

		# textures and bg
		bg_map_l = fsaverage.sulc_left
		mesh_l = fsaverage.infl_left
		texture_l = surface.vol_to_surf(img, fsaverage.pial_left)
		bg_map_r = fsaverage.sulc_right
		mesh_r = fsaverage.infl_right
		texture_r = surface.vol_to_surf(img, fsaverage.pial_right)

		# top
		plot_surf_stat_map(mesh_l, texture_l, hemi='left', colorbar=False, view='lateral', cmap=colormap,
		                   bg_map=bg_map_l, threshold=threshold, vmin=threshold, figure=fig, axes=ax_surf_top_l)
		plot_surf_stat_map(mesh_l, texture_l, hemi='left', colorbar=False, view='lateral', cmap=colormap,
		                   bg_map=bg_map_l, threshold=threshold, vmin=threshold, figure=fig, axes=ax_surf_top_ml)
		plot_surf_stat_map(mesh_r, texture_r, hemi='right', colorbar=False, view='lateral', cmap=colormap,
		                   bg_map=bg_map_r, threshold=threshold, vmin=threshold, figure=fig, axes=ax_surf_top_mr)
		plot_surf_stat_map(mesh_r, texture_r, hemi='right', colorbar=False, view='lateral', cmap=colormap,
		                   bg_map=bg_map_r, threshold=threshold, vmin=threshold, figure=fig, axes=ax_surf_top_r)
		# bottom
		plot_surf_stat_map(mesh_l, texture_l, hemi='left', colorbar=False, view='medial', cmap=colormap,
		                   bg_map=bg_map_l, threshold=threshold, vmin=threshold, figure=fig, axes=ax_surf_bottom_l)
		plot_surf_stat_map(mesh_l, texture_l, hemi='left', colorbar=False, view='posterior', cmap=colormap,
		                   bg_map=bg_map_l, threshold=threshold, vmin=threshold, figure=fig, axes=ax_surf_bottom_ml)
		plot_surf_stat_map(mesh_r, texture_r, hemi='right', colorbar=False, view='posterior', cmap=colormap,
		                   bg_map=bg_map_r, threshold=threshold, vmin=threshold, figure=fig, axes=ax_surf_bottom_mr)
		plot_surf_stat_map(mesh_r, texture_r, hemi='right', colorbar=False, view='medial', cmap=colormap,
		                   bg_map=bg_map_r, threshold=threshold, vmin=threshold, figure=fig, axes=ax_surf_bottom_r)

		scaling = 85

		ax_surf_top_ml.view_init(90, -90)
		ax_surf_top_ml.set_xlim3d(-scaling, scaling)
		ax_surf_top_ml.set_ylim3d(-scaling, scaling)
		pos_ml = ax_surf_top_ml.get_position()
		ax_surf_top_ml.set_position([pos_ml.bounds[0] + 0.035, pos_ml.bounds[1], pos_ml.bounds[2], pos_ml.bounds[3]])
		ax_surf_top_mr.view_init(90, -90)
		ax_surf_top_mr.set_xlim3d(-scaling, scaling)
		ax_surf_top_mr.set_ylim3d(-scaling, scaling)
		pos_mr = ax_surf_top_mr.get_position()
		ax_surf_top_mr.set_position([pos_mr.bounds[0] - 0.035, pos_mr.bounds[1], pos_mr.bounds[2], pos_mr.bounds[3]])

		ax_surf_bottom_ml.view_init(-90, -90)
		ax_surf_bottom_ml.set_xlim3d(-scaling, scaling)
		ax_surf_bottom_ml.set_ylim3d(-scaling, scaling)
		pos_ml = ax_surf_bottom_ml.get_position()
		ax_surf_bottom_ml.set_position([pos_ml.bounds[0] + 0.035, pos_ml.bounds[1], pos_ml.bounds[2], pos_ml.bounds[3]])
		ax_surf_bottom_mr.view_init(-90, -90)
		ax_surf_bottom_mr.set_xlim3d(-scaling, scaling)
		ax_surf_bottom_mr.set_ylim3d(-scaling, scaling)
		pos_mr = ax_surf_bottom_mr.get_position()
		ax_surf_bottom_mr.set_position([pos_mr.bounds[0] - 0.035, pos_mr.bounds[1], pos_mr.bounds[2], pos_mr.bounds[3]])

		pos_top_l = ax_surf_top_l.get_position()
		ax_surf_top_l.set_position([pos_top_l.bounds[0] + 0.035, pos_top_l.bounds[1], pos_top_l.bounds[2], pos_top_l.bounds[3]])

		pos_top_r = ax_surf_top_r.get_position()
		ax_surf_top_r.set_position(
			[pos_top_r.bounds[0] - 0.035, pos_top_r.bounds[1], pos_top_r.bounds[2], pos_top_r.bounds[3]])

		pos_bottom_l = ax_surf_bottom_l.get_position()
		ax_surf_bottom_l.set_position(
			[pos_bottom_l.bounds[0] + 0.035, pos_bottom_l.bounds[1], pos_bottom_l.bounds[2], pos_bottom_l.bounds[3]])

		pos_bottom_r = ax_surf_bottom_r.get_position()
		ax_surf_bottom_r.set_position(
			[pos_bottom_r.bounds[0] - 0.035, pos_bottom_r.bounds[1], pos_bottom_r.bounds[2], pos_bottom_r.bounds[3]])

		# construct cmap
		norm = mpl.colors.Normalize(vmin=threshold, vmax=max_val)
		cmap_orig = mpl.cm.get_cmap(colormap)

		#fig.subplots_adjust(top=0.9)
		cbar_ax = fig.add_axes([0.2, 0.1, 0.625, 0.02])
		cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap_orig,
		                                norm=norm,
		                                orientation='horizontal')
		cb1.ax.tick_params(labelsize=20)

		# add text
		ax_surf_top_l.text2D(0.05, 0.085, 'Left hemisphere', fontsize=20, color='slategray', ha='center')
		ax_surf_top_r.text2D(-.05, 0.085, 'Right hemisphere', fontsize=20, color='slategray', ha='center')

		#bbox = fig.bbox
		plt.savefig(opj(output_dir, 'contrast-' + con + '.' + output_format), fig=fig, dpi=dpi, format=output_format, bbox_inches='tight')
		plt.close()
		print('Written image: ' + str(opj(output_dir, 'contrast-' + con + '.' + output_format)))


def plot_coord_atlas_bg(con, x, coords, bg):

	conDir = opj(output_dir, 'space-' + 'MNI152NLin2009cAsym', 'contrast-' + con)

	img = load_img(opj(conDir, 'spmT_0001.nii'))
	[img_thr, p] = threshold_stats_img(stat_img=img, alpha=0.05, height_control='fdr', cluster_threshold=20)
	coord_mni = coords
	plot_stat_map(img_thr, bg_img=bg, title='', threshold=threshold, colorbar=True,
	                        vmax=13, vmin=threshold, cut_coords=coord_mni,
	                        cmap='viridis', draw_cross=True, dim=False)

	plt.savefig(opj(output_dir, 'contrast-' + con + '_cluster-' + str(x) + '.' + output_format), dpi=dpi, format=output_format)
	plt.close()

def plot_decode_bg(toDecode, x, coords, bg):

	decDir = opj(decode_dir, toDecode, 'snpm_batch')

	img = load_img(opj(decDir, 'SnPM_filtered_masked.nii'))

	plot_stat_map(img, bg_img=bg, title='', threshold=2.9, colorbar=True,
	              vmax=15, vmin=2.9, cut_coords=coords,
	              cmap='viridis', draw_cross=True, dim=False)

	plt.savefig(opj(decode_dir, 'decode-' + toDecode + '_cluster-' + str(x) + '.' + output_format), dpi=dpi,
	            format=output_format)
	plt.close()

if __name__ == "__main__":
	args = get_args()

	if 'surf' in args.type:
		plotContrastSurf(args.contrasts)

	if 'coord_bg' in args.type:

		bg_img_folder = opj('/home/achtzehnj/data/atlases/100um/Synthesized_FLASH25_in_MNI_v2_200um.nii.gz')
		bg_img = load_img(bg_img_folder)
		bg_img = threshold_img(bg_img, 9.)

		coords = [[24, -84, 44]]
		for x, coord in enumerate(coords):
			plot_coord_atlas_bg(args.contrasts[0], x, coord, bg_img)

	if 'decode_bg' in args.type:

		bg_img_folder = opj('/home/achtzehnj/data/atlases/100um/Synthesized_FLASH25_in_MNI_v2_200um.nii.gz')
		bg_img = load_img(bg_img_folder)
		bg_img = threshold_img(bg_img, 9.)

		coords = [[45, -63, 1], [-48, -72, 1], [0, -84, 28]]
		decode = 'dist-vs-dots'

		for x, coord in enumerate(coords):
			plot_decode_bg(decode, x, coord, bg_img)

# 2nd level
# 15 [[42, -51, -19], [-54, -78, -2], [39, 48, -15], [48, -36, 54], [30, -72, 34], [24, -84, 44]], num > contr
# 13 [[36, 45, -12], [18, -66, -2], [45, -81, 24]] dist > contr
# 11 [[-45, 45, -12], [36, 45, -12], [-54, 27, 1]] time > contr

# decode
# [[9, -87, -9], [48, 45, -9], [-48, 45, -9], [9, -30, -9]] time-vs-lumin
# [[12, -72, -2], [3, -63, 44], [36, 48, -15], [0, -30, -5]] dist-vs-lumin
# [[45, -84, -5], [-48, -87, 1], [27, -69, 34], [0, -72, 54], [-3, -18, 14], [42, 48, -12], [-3, -39, 28]] dots-vs-lumin

# [[45, -66, -2], [-48, -78, 1], [36, -57, 57], [-36, -60, 57], [6, -87, 13]] time-vs-dots
# [[45, -63, 1], [-48, -72, 1], [0, -84, 28]] dist-vs-dots