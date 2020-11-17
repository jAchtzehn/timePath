from os.path import join as opj
from os.path import abspath
from nilearn.plotting import plot_glass_brain, plot_surf_stat_map
from nilearn.image import coord_transform
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

# ---- options ----
addMarkers = False
contrast_list = range(1, 26)
colormap = 'viridis'
output_format = 'pdf'
dpi = 300

contrast_ids = [('%04d' % x) for x in contrast_list]
condition_names = ['time', 'dist', 'lumin', 'dots', 'icon', 'comp']
contrasts = [['1-trial-avg', 'T', condition_names, [0.25, 0.25, 0.25, 0.25, 0, 0]],
                 ['2-TDD-avg', 'T', condition_names, [1/3, 1/3, 0, 1/3, 0, 0]],
                 ['3-time', 'T', condition_names, [1, 0, 0, 0, 0, 0]],
                 ['4-dist', 'T', condition_names, [0, 1, 0, 0, 0, 0]],
                 ['5-lumin', 'T', condition_names, [0, 0, 1, 0, 0, 0]],
                 ['6-dots', 'T', condition_names, [0, 0, 0, 1, 0, 0]],
                 ['7-icon', 'T', condition_names, [0, 0, 0, 0, 1, 0]],
                 ['8-comp', 'T', condition_names, [0, 0, 0, 0, 0, 1]],
                 ['9-TDD>lumin', 'T', condition_names, [1/3., 1/3., -1, 1/3., 0, 0]],
                 ['10-lumin>TDD', 'T', condition_names, [-1/3., -1/3., 1, -1/3., 0, 0]],
                 ['11-time>lumin', 'T', condition_names, [1, 0, -1, 0, 0, 0]],
                 ['12-lumin>time', 'T', condition_names, [-1, 0, 1, 0, 0, 0]],
                 ['13-dist>lumin', 'T', condition_names, [0, 1, -1, 0, 0, 0]],
                 ['14-lumin>dist', 'T', condition_names, [0, -1, 1, 0, 0, 0]],
                 ['15-dots>lumin', 'T', condition_names, [0, 0, -1, 1, 0, 0]],
                 ['16-lumin>dots', 'T', condition_names, [0, 0, 1, -1, 0, 0]],
                 ['17-time>others', 'T', condition_names, [1, -0.5, 0, -0.5, 0, 0]],
                 ['18-dist>others', 'T', condition_names, [-0.5, 1, 0, -0.5, 0, 0]],
                 ['19-dots>others', 'T', condition_names, [-0.5, -0.5, 0, 1, 0, 0]],
                 ['20-time>dist', 'T', condition_names, [1, -1, 0, 0, 0, 0]],
                 ['21-dist>time', 'T', condition_names, [-1, 1, 0, 0, 0, 0]],
                 ['22-time>dots', 'T', condition_names, [1, 0, 0, -1, 0, 0]],
                 ['23-dots>time', 'T', condition_names, [-1, 0, 0, 1, 0, 0]],
                 ['24-dist>dots', 'T', condition_names, [0, 1, 0, -1, 0, 0]],
                 ['25-dots>dist', 'T', condition_names, [0, -1, 0, 1, 0, 0]]
             ]

output_dir = opj(experiment_dir, 'output_nipype', '2ndLevel')


mpl.rc('text', usetex=True)
pl.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',    # helvetica font
    r'\usepackage{sansmath}',   # math-font matching  helvetica
    r'\sansmath'                # actually tell tex to use it!
    r'\usepackage{siunitx}',    # micro symbols
    r'\sisetup{detect-all}',    # force siunitx to use the fonts
]
mpl.rcParams['ytick.labelsize'] = 22
mpl.rcParams['xtick.labelsize'] = 22
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['axes.labelsize'] = 28
mpl.rcParams['axes.labelpad'] = 12
mpl.rcParams['axes.titlesize'] = 32
mpl.rcParams['axes.titlepad'] = 20

for x, con in enumerate(contrast_ids):

	# init figure with gridspec
	fig = plt.figure(figsize=(15, 15), constrained_layout=False)
	gs = fig.add_gridspec(3, 2)
	ax_glass = fig.add_subplot(gs[0, :])   # top row
	ax_surf_top_l = fig.add_subplot(gs[1, 0], projection='3d')
	ax_surf_top_r = fig.add_subplot(gs[1, 1], projection='3d')
	ax_surf_bottom_l = fig.add_subplot(gs[2, 0], projection='3d')
	ax_surf_bottom_r = fig.add_subplot(gs[2, 1], projection='3d')

	conDir = opj(output_dir, 'space-' + 'MNI152NLin2009cAsym', 'contrast-' + con)

	# glass brain
	# find maximum for plot thresholding
	img = load_img(opj(conDir, 'spmT_0001_thr.nii'))
	img_data = _safe_get_data(img, ensure_finite=True)

	max_val = ndimage.maximum(img_data)                 # multiply this value with the percent threshold to get absolute value of threshold
	if max_val > 0:
		min_val = ndimage.minimum(img_data[img_data > 0])
	else:
		min_val = max_val

	img_no_nans = threshold_img(img, 0)

	plot_glass_brain(img_no_nans, display_mode='lyrz', colorbar=False, plot_abs=False,
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

	# construct cmap
	norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
	cmap_orig = mpl.cm.get_cmap(colormap)

	# fig.subplots_adjust(top=0.9)
	cbar_ax = fig.add_axes([0.15, 0.925, 0.7, 0.05])
	cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap_orig,
	                                norm=norm,
	                                orientation='horizontal')

	plt.tight_layout()
	plt.savefig(opj(output_dir, contrasts[x][0] + '.' + output_format), fig=fig, dpi=dpi, format=output_format)
	print('Written image: ' + str(opj(output_dir, contrasts[x][0] + '.' + output_format)))
