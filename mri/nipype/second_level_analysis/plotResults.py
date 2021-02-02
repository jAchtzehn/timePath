from os.path import join as opj
from os.path import abspath
from nilearn.plotting import plot_glass_brain, plot_surf_stat_map, plot_stat_map
from nilearn.image import coord_transform
from nilearn import datasets
from nilearn import surface
from scipy import ndimage
from sys import platform
import matplotlib.pyplot as plt
import matplotlib as mpl
import pylab as pl
from mpl_toolkits.mplot3d import axes3d

from nilearn.image import load_img, new_img_like, threshold_img
from nilearn._utils.niimg_conversions import _safe_get_data
import numpy as np

if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/home/achtzehnj/data/timePath/')

# ---- options ----
addMarkers = False
contrast_list = range(2, 3)
colormap = 'viridis'
output_format = 'png'
dpi = 350
plot_threshold = 6.39
fs_average = 'fsaverage5'

contrast_ids = [('%04d' % x) for x in contrast_list]
condition_names = ['time', 'dist', 'lumin', 'dots', 'icon', 'comp']
contrasts = [['2-TDD-avg', 'T', condition_names, [1/3, 1/3, 0, 1/3, 0, 0]]]
output_dir = opj(experiment_dir, 'derivatives', 'output_nipype', '2ndlevel_higher_cluster')
thesis_dir = abspath('/Users/jachtzehn/Documents/Medizin/thesis/figures/results/fmri')
thesis_dir_appendix = abspath('/Users/jachtzehn/Documents/Medizin/thesis/figures/appendix')

# contrasts = [['1-trial-avg', 'T', condition_names, [0.25, 0.25, 0.25, 0.25, 0, 0]],
#                  ['2-TDD-avg', 'T', condition_names, [1/3, 1/3, 0, 1/3, 0, 0]],
#                  ['3-time', 'T', condition_names, [1, 0, 0, 0, 0, 0]],
#                  ['4-dist', 'T', condition_names, [0, 1, 0, 0, 0, 0]],
#                  ['5-lumin', 'T', condition_names, [0, 0, 1, 0, 0, 0]],
#                  ['6-dots', 'T', condition_names, [0, 0, 0, 1, 0, 0]],
#                  ['7-icon', 'T', condition_names, [0, 0, 0, 0, 1, 0]],
#                  ['8-comp', 'T', condition_names, [0, 0, 0, 0, 0, 1]],
#                  ['9-TDD-gt-lumin', 'T', condition_names, [1/3., 1/3., -1, 1/3., 0, 0]],
#                  ['10-lumin-gt-TDD', 'T', condition_names, [-1/3., -1/3., 1, -1/3., 0, 0]],
#                  ['11-time-gt-lumin', 'T', condition_names, [1, 0, -1, 0, 0, 0]],
#                  ['12-lumin-gt-time', 'T', condition_names, [-1, 0, 1, 0, 0, 0]],
#                  ['13-dist-gt-lumin', 'T', condition_names, [0, 1, -1, 0, 0, 0]],
#                  ['14-lumin-gt-dist', 'T', condition_names, [0, -1, 1, 0, 0, 0]],
#                  ['15-dots-gt-lumin', 'T', condition_names, [0, 0, -1, 1, 0, 0]],
#                  ['16-lumin-gt-dots', 'T', condition_names, [0, 0, 1, -1, 0, 0]],
#                  ['17-time-gt-others', 'T', condition_names, [1, -0.5, 0, -0.5, 0, 0]],
#                  ['18-dist-gt-others', 'T', condition_names, [-0.5, 1, 0, -0.5, 0, 0]],
#                  ['19-dots-gt-others', 'T', condition_names, [-0.5, -0.5, 0, 1, 0, 0]],
#                  ['20-time-gt-dist', 'T', condition_names, [1, -1, 0, 0, 0, 0]],
#                  ['21-dist-gt-time', 'T', condition_names, [-1, 1, 0, 0, 0, 0]],
#                  ['22-time-gt-dots', 'T', condition_names, [1, 0, 0, -1, 0, 0]],
#                  ['23-dots-gt-time', 'T', condition_names, [-1, 0, 0, 1, 0, 0]],
#                  ['24-dist-gt-dots', 'T', condition_names, [0, 1, 0, -1, 0, 0]],
#                  ['25-dots-gt-dist', 'T', condition_names, [0, -1, 0, 1, 0, 0]]
#              ]

# mpl.rc('text', usetex=True)
# pl.rcParams['text.latex.preamble'] = [
#     r'\usepackage{tgheros}',    # helvetica font
#     r'\usepackage{sansmath}',   # math-font matching  helvetica
#     r'\sansmath'                # actually tell tex to use it!
#     r'\usepackage{siunitx}',    # micro symbols
#     r'\sisetup{detect-all}',    # force siunitx to use the fonts
# ]
# mpl.rcParams['ytick.labelsize'] = 18
# mpl.rcParams['xtick.labelsize'] = 22
# mpl.rcParams['axes.linewidth'] = 1
# mpl.rcParams['axes.labelsize'] = 28
# mpl.rcParams['axes.labelpad'] = 12
# mpl.rcParams['axes.titlesize'] = 32
# mpl.rcParams['axes.titlepad'] = 20


def plotContrastVolumes():
	for x, con in enumerate(contrast_ids):

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

		conDir = opj(output_dir, 'space-' + 'MNI152NLin2009cAsym', 'contrast-' + con)

		# glass brain
		# find maximum for plot thresholding
		img = load_img(opj(conDir, 'spmT_0001_thr.nii'))
		img_data = _safe_get_data(img, ensure_finite=True)

		max_val = ndimage.maximum(img_data)                 # multiply this value with the percent threshold to get absolute value of threshold
		if max_val > 0:
			min_val = plot_threshold
		else:
			min_val = max_val

		# 3D plots
		fsaverage = datasets.fetch_surf_fsaverage(fs_average)
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
		                   bg_map=bg_map_l, threshold=plot_threshold, figure=fig, axes=ax_surf_top_l)
		plot_surf_stat_map(mesh_l, texture_l, hemi='left', colorbar=False, view='lateral', cmap=colormap,
		                   bg_map=bg_map_l, threshold=plot_threshold, figure=fig, axes=ax_surf_top_ml)
		plot_surf_stat_map(mesh_r, texture_r, hemi='right', colorbar=False, view='lateral', cmap=colormap,
		                   bg_map=bg_map_r, threshold=plot_threshold, figure=fig, axes=ax_surf_top_mr)
		plot_surf_stat_map(mesh_r, texture_r, hemi='right', colorbar=False, view='lateral', cmap=colormap,
		                   bg_map=bg_map_r, threshold=plot_threshold, figure=fig, axes=ax_surf_top_r)
		# bottom
		plot_surf_stat_map(mesh_l, texture_l, hemi='left', colorbar=False, view='medial', cmap=colormap,
		                   bg_map=bg_map_l, threshold=plot_threshold, figure=fig, axes=ax_surf_bottom_l)
		plot_surf_stat_map(mesh_l, texture_l, hemi='left', colorbar=False, view='posterior', cmap=colormap,
		                   bg_map=bg_map_l, threshold=plot_threshold, figure=fig, axes=ax_surf_bottom_ml)
		plot_surf_stat_map(mesh_r, texture_r, hemi='right', colorbar=False, view='posterior', cmap=colormap,
		                   bg_map=bg_map_r, threshold=plot_threshold, figure=fig, axes=ax_surf_bottom_mr)
		plot_surf_stat_map(mesh_r, texture_r, hemi='right', colorbar=False, view='medial', cmap=colormap,
		                   bg_map=bg_map_r, threshold=plot_threshold, figure=fig, axes=ax_surf_bottom_r)

		scaling = 85

		ax_surf_top_ml.view_init(90, -90)
		ax_surf_top_ml.set_xlim3d(-scaling, scaling)
		ax_surf_top_ml.set_ylim3d(-scaling, scaling)
		pos_ml = ax_surf_top_ml.get_position()
		ax_surf_top_ml.set_position([pos_ml.bounds[0] + 0.03, pos_ml.bounds[1], pos_ml.bounds[2], pos_ml.bounds[3]])
		ax_surf_top_mr.view_init(90, -90)
		ax_surf_top_mr.set_xlim3d(-scaling, scaling)
		ax_surf_top_mr.set_ylim3d(-scaling, scaling)
		pos_mr = ax_surf_top_mr.get_position()
		ax_surf_top_mr.set_position([pos_mr.bounds[0] - 0.03, pos_mr.bounds[1], pos_mr.bounds[2], pos_mr.bounds[3]])

		ax_surf_bottom_ml.view_init(-90, -90)
		ax_surf_bottom_ml.set_xlim3d(-scaling, scaling)
		ax_surf_bottom_ml.set_ylim3d(-scaling, scaling)
		pos_ml = ax_surf_bottom_ml.get_position()
		ax_surf_bottom_ml.set_position([pos_ml.bounds[0] + 0.03, pos_ml.bounds[1], pos_ml.bounds[2], pos_ml.bounds[3]])
		ax_surf_bottom_mr.view_init(-90, -90)
		ax_surf_bottom_mr.set_xlim3d(-scaling, scaling)
		ax_surf_bottom_mr.set_ylim3d(-scaling, scaling)
		pos_mr = ax_surf_bottom_mr.get_position()
		ax_surf_bottom_mr.set_position([pos_mr.bounds[0] - 0.03, pos_mr.bounds[1], pos_mr.bounds[2], pos_mr.bounds[3]])

		# construct cmap
		norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
		cmap_orig = mpl.cm.get_cmap(colormap)

		#fig.subplots_adjust(top=0.9)
		cbar_ax = fig.add_axes([0.2, 0.1, 0.625, 0.02])
		cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap_orig,
		                                norm=norm,
		                                orientation='horizontal')
		cb1.ax.tick_params(labelsize=20)

		# add text
		ax_surf_top_l.text2D(0.05, 0.085, 'Left hemisphere', fontsize=20, color='slategray', ha='center')
		ax_surf_top_ml.text2D(0.029, 0.065, 'A', fontsize=20, color='black', ha='center')
		ax_surf_top_ml.text2D(0.029, -0.08, 'P', fontsize=20, color='black', ha='center')

		ax_surf_bottom_ml.text2D(0.029, -0.075, 'A', fontsize=20, color='black', ha='center')
		ax_surf_top_r.text2D(-.05, 0.085, 'Right hemisphere', fontsize=20, color='slategray', ha='center')

		#bbox = fig.bbox
		plt.savefig(opj(output_dir, contrasts[x][0] + '.' + output_format), fig=fig, dpi=dpi, format=output_format, bbox_inches='tight')
		plt.close()
		print('Written image: ' + str(opj(output_dir, contrasts[x][0] + '.' + output_format)))

		fig_gb = plt.figure(figsize=(14, 4))
		img_no_nans = threshold_img(img, 0)
		plot_glass_brain(img_no_nans, display_mode='lyrz', colorbar=True, plot_abs=False,
		                 threshold=plot_threshold, vmax=max_val, cmap=colormap, figure=fig_gb)
		plt.tight_layout()
		plt.savefig(opj(output_dir, contrasts[x][0] + '_gb.' + 'pdf'), fig=fig_gb, dpi=dpi, format='pdf')
		plt.close()


def plotSlice(contrastNr, contrastInfo, threshold, plot_axis, cluster_coords, transform_to_mni):

	conDir = opj(output_dir, 'space-' + 'MNI152NLin2009cAsym', 'contrast-' + str(contrastNr).zfill(4))

	# find maximum for plot thresholding
	img = load_img(opj(conDir, 'spmT_0001_thr.nii'))
	img_data = _safe_get_data(img, ensure_finite=True)
	anat_img = load_img(opj(experiment_dir, 'templates', 'mni_1mm', '1mm_T1.nii.gz'))

	max_val = ndimage.maximum(
		img_data)  # multiply this value with the percent threshold to get absolute value of threshold
	if max_val > 0:
		min_val = threshold
	else:
		min_val = max_val

	fig_anat = plt.figure(figsize=(16, 10))

	# transform coord
	if transform_to_mni:
		cut_coords_mni = ()
		for cluster in cluster_coords:
			new_tuple = coord_transform(0, cluster[1], 0, img.affine)
			cut_coords_mni = cut_coords_mni + (new_tuple[1],)
	else:
		cut_coords_mni = ()
		for cluster in cluster_coords:
			cut_coords_mni = cut_coords_mni + (cluster[1],)

	display = plot_stat_map(img, bg_img=anat_img, title='', dim=1, threshold=threshold,
	                        vmax=max_val, display_mode=plot_axis, cut_coords=cut_coords_mni,
	                        cmap='viridis')

	# add markers
	# display.add_markers(np.array(cluster_coords), marker_color='red', marker_size=20)  # add marker

	plt.tight_layout()
	plt.savefig(opj(output_dir, contrastInfo[contrastNr - 1][0] + '-anat-plot.' + 'pdf'), fig=fig_anat,
	            dpi=300, format='pdf')
	# plt.savefig(opj(thesis_dir, contrastInfo[contrastNr - 1][0] + '-anat-plot.' + 'pdf'), fig=fig_anat,
	#             dpi=300, format='pdf')
	plt.close()
	print('Written image: {}'.format(opj(output_dir, contrastInfo[contrastNr - 1][0] + '-anat-plot.' + output_format)))

	fig_gb = plt.figure(figsize=(14, 4))
	img_no_nans = threshold_img(img, 0)
	plot_glass_brain(img_no_nans, display_mode='lyrz', colorbar=True, plot_abs=False,
	                 threshold=threshold, vmax=max_val, cmap='viridis', figure=fig_gb)
	plt.tight_layout()
	plt.savefig(opj(thesis_dir_appendix, contrastInfo[contrastNr - 1][0] + '_gb' + '.pdf'),
	            fig=fig_gb, dpi=300, format='pdf')
	plt.close()


# plotContrastVolumes()

clusters_to_plot = [[15, [[42, -51, -19], [-54, -78, -2], [39, 48, -15], [30, -72, 34], [48, -36, 54]]],
                    [13, [[36, 45, -12], [36, 44, -12], [36, 46, -12], [18, -66, -2], [18, -65, -2]]],
                    [9, [[39, 45, -12], [51, -75, -2], [-42, 42, -12], [42, -54, -19], [-51, -81, -5]]],
                    [11, [[-45, 45, -12], [-45, 44, -12], [36, 46, -12], [36, 47, -12], [-54, 27, 1]]],
                    [2, [[36, 45, -12], [36, 44, -12], [36, 46, -12], [18, -66, -2], [18, -65, -2]]],
                    [17, [[-39, 57, 4], [30, 42, 24], [30, 48, 24], [6, -84, 28], [-6, -75, 34]]],
                    [18, [[-9, -55, 64], [-9, -56, 64], [-9, -53, 64], [-45, -81, 24], [-45, -83, 24]]],
                    [19, [[48, -66, 1], [48, -64, 1], [-51, -75, 1], [36, -54, 54], [36, -52, 54]]]]


#for cluster in clusters_to_plot:
#	plotSlice(cluster[0], contrasts, plot_threshold, 'y', cluster[1], False)


def plotSliceImg(imgFilename, name, threshold, plot_axis, cluster_coords, transform_to_mni):

	# find maximum for plot thresholding
	img = load_img(imgFilename)
	img_data = _safe_get_data(img, ensure_finite=True)
	anat_img = load_img(opj(experiment_dir, 'templates', 'mni_1mm', '1mm_T1.nii.gz'))

	max_val = ndimage.maximum(
		img_data)  # multiply this value with the percent threshold to get absolute value of threshold
	if max_val > 0:
		min_val = threshold
	else:
		min_val = max_val

	fig_anat = plt.figure(figsize=(16, 10))

	# transform coord
	if transform_to_mni:
		cut_coords_mni = ()
		for cluster in cluster_coords:
			new_tuple = coord_transform(0, cluster[1], 0, img.affine)
			cut_coords_mni = cut_coords_mni + (new_tuple[1],)
	else:
		cut_coords_mni = ()
		for cluster in cluster_coords:
			cut_coords_mni = cut_coords_mni + (cluster[1],)

	display = plot_stat_map(img, bg_img=anat_img, title='', dim=1, threshold=threshold,
	                        vmax=max_val, display_mode=plot_axis, cut_coords=cut_coords_mni,
	                        cmap='viridis')

	# add markers
	# display.add_markers(np.array(cluster_coords), marker_color='red', marker_size=20)  # add marker

	plt.tight_layout()
	plt.savefig(opj(output_dir, name + '-anat-plot.' + 'pdf'), fig=fig_anat,
	            dpi=300, format='pdf')
	plt.close()
	print('Written image: {}'.format(opj(output_dir, name + '-anat-plot.' + output_format)))

	fig_gb = plt.figure(figsize=(14, 4))
	img_no_nans = threshold_img(img, 0)
	plot_glass_brain(img_no_nans, display_mode='lyrz', colorbar=True, plot_abs=False,
	                 threshold=threshold, vmax=max_val, cmap='viridis', figure=fig_gb)
	plt.tight_layout()
	plt.savefig(opj(thesis_dir_appendix, name + '_gb' + '.pdf'),
	            fig=fig_gb, dpi=300, format='pdf')
	plt.close()

fname = opj(experiment_dir, 'conjunction', 'conjAll_thrSPM.nii')
name = 'conj'

plotSliceImg(fname, name, 0, 'y', [[42, 51, -12], [42, 47, -12], [42, 45, -12], [42, 43, -12], [42, 41, -12]], False)

