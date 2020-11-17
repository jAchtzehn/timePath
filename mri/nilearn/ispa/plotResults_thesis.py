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
from nilearn.plotting import plot_stat_map
from nilearn.image import load_img, new_img_like, threshold_img, math_img
from nilearn._utils.niimg_conversions import _safe_get_data
import numpy as np
from nistats.reporting import compare_niimgs

if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/mnt/work/achtzehnj/data')

nilearn_dir = opj(experiment_dir, 'nilearn')
results_dir = opj(nilearn_dir, 'group_results_ispa_std')
thesis_dir = abspath('/Users/jachtzehn/Documents/Medizin/thesis/figures/results/fmri')

# ---- options ----
colormap = 'viridis'
output_name = '_group_level_mvpa_masked_gb'
output_format = 'pdf'
dpi = 300
conditions_to_decode = [['time', 'dist'], ['time', 'dots'], ['dist', 'dots'], ['time', 'lumin'], ['dist', 'lumin'], ['dots', 'lumin']]
mask_img = True     # mask image with cluster mask?
nrSplits = 24

mpl.rc('text', usetex=True)
pl.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',    # helvetica font
    r'\usepackage{sansmath}',   # math-font matching  helvetica
    r'\sansmath'                # actually tell tex to use it!
    r'\usepackage{siunitx}',    # micro symbols
    r'\sisetup{detect-all}',    # force siunitx to use the fonts
]
mpl.rcParams['ytick.labelsize'] = 18
mpl.rcParams['xtick.labelsize'] = 22
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['axes.labelsize'] = 28
mpl.rcParams['axes.labelpad'] = 12
mpl.rcParams['axes.titlesize'] = 32
mpl.rcParams['axes.titlepad'] = 20


def plotGlassBrains():
	for condition_pair in conditions_to_decode:

		if len(condition_pair) > 1:
			condition_dir = opj(results_dir, condition_pair[0] + '-vs-' + condition_pair[1], 'snpm_batch')
			output_fname = condition_pair[0] + '-vs-' + condition_pair[1] + output_name
		else:
			condition_dir = opj(results_dir, condition_pair[0], 'snpm_batch')
			output_fname = condition_pair[0] + output_name

		img_fname = opj(condition_dir, 'SnPM_filtered.nii')
		mask_fname = opj(condition_dir, 'SnPM_filtered_binary_mask_cluster.nii')
		mask = load_img(mask_fname)

		# 1. plot significant cluster
		fig, ax = plt.subplots(1, 1, figsize=(15, 4), constrained_layout=False)

		if len(condition_pair) > 1:
			split_dir = opj(results_dir, condition_pair[0] + '-vs-' + condition_pair[1])
			output_fname = condition_pair[0] + '-vs-' + condition_pair[1] + output_name + '_accuracies'
		else:
			split_dir = opj(results_dir, condition_pair[0])
			output_fname = condition_pair[0] + output_name + '_accuracies'

		img = mean_img(opj(split_dir, 'ispa_sl_wb_decoding-' + condition_pair[0] + '-vs-' + condition_pair[1] + '_split-*.nii'))

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
		print(max_val)
		print(min_val)

		# construct cmap
		norm = mpl.colors.Normalize(vmin=0.0292, vmax=0.15)
		cmap_orig = mpl.cm.get_cmap(colormap)

		plot_glass_brain(img, display_mode='lyrz', colorbar=False, vmin=0.0292, cmap=cmap_orig, figure=fig, axes=ax)

		# add custom cbar
		fig.subplots_adjust(right=0.8)
		cbar_ax = fig.add_axes([0.95, 0.15, 0.01, 0.7])
		cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap_orig,
		                                norm=norm,
		                                orientation='vertical')
		cbar_ax.yaxis.set_ticks_position('left')
		cbar_ax.yaxis.set_ticks([0.0292, 0.05, 0.075, 0.1, 0.125, 0.15])
		cbar_ax.yaxis.set_ticklabels([r'2.9$\%$', r'5$\%$', r'7.5$\%$', r'10$\%$', r'12.5$\%$', r'15$\%$'])

		#plt.tight_layout()
		plt.savefig(opj(results_dir, output_fname + '.' + output_format), fig=fig, dpi=dpi, format=output_format)
		plt.savefig(opj(thesis_dir, condition_pair[0] + '-vs-' + condition_pair[1] + '.' + output_format), fig=fig, dpi=dpi, format=output_format, bbox_inches='tight')
		print('Written image: ' + str(opj(results_dir, output_fname + '.' + output_format)))
		plt.close()


# plotGlassBrains()

def plotAnatSlices(condition, coordinates, clustername):

	condition_dir = opj(results_dir, condition[0] + '-vs-' + condition[1], 'snpm_batch')
	output_fname = condition[0] + '-vs-' + condition[1] + '-' + clustername

	img = load_img(opj(condition_dir, 'SnPM_filtered_masked_accuracy.nii'))
	anat_img = load_img(opj(experiment_dir, 'templates', 'mni_1mm', '1mm_T1.nii.gz'))

	fig, ax = plt.subplots(1, 1, figsize=(15, 4), constrained_layout=False)

	# find maximum for plot thresholding
	img_data = _safe_get_data(img, ensure_finite=True)

	max_val = ndimage.maximum(
		img_data)  # multiply this value with the percent threshold to get absolute value of threshold
	if max_val > 0:
		min_val = ndimage.minimum(img_data[img_data > 0])
	else:
		min_val = max_val

	img = math_img("img * 100", img=img)

	display = plot_stat_map(img, bg_img=anat_img, title='', draw_cross=True, colorbar=True,
	                        vmax=15, vmin=2.9, threshold=0, display_mode='ortho', cut_coords=coordinates,
	                        cmap='viridis')

	plt.savefig(opj(thesis_dir, output_fname + '.pdf'), fig=fig,
	            dpi=300, format='pdf', bbox_inches='tight')
	print('Written image: {}'.format(opj(thesis_dir, output_fname + '.pdf')))
	plt.close()


#plotAnatSlices(['dots', 'lumin'], (-2, -11, 8), 'thalamus')
#plotAnatSlices(['dots', 'lumin'], (40, 48, -11), 'ifg')
#plotAnatSlices(['time', 'lumin'], (-48, 45, -9), 'time-unique')
#plotAnatSlices(['time', 'lumin'], (9, -87, -9), 'time-unique2')
#plotAnatSlices(['dist', 'lumin'], (9, -63, -5), 'dist-unique')
#plotAnatSlices(['dots', 'lumin'], (42, -63, -15), 'dots-unique')
#plotAnatSlices(['dist', 'dots'], (48, -66, -2), 'fg')
#plotAnatSlices(['time', 'dots'], (-6, -75, 37), 'precuneus')
#plotAnatSlices(['time', 'dots'], (-36, -54, 54), 'ips')
#plotAnatSlices(['time', 'dots'], (-64, -45, 37), 'supramarg')

def compareNiiImgs(img1, img2):
	#masker = NiftiMasker(mask_img=mask_filename, standardize=True, detrend=False)       # create nifti masker (2D array ready) from mask file

	compare_niimgs(img1, img2, masker=mask,
	               ref_label='img1', src_label='img2')

	plt.show()


#compareNiiImgs(opj(results_dir, 'dist-vs-dots', 'snpm_batch', 'SnPM_filtered_masked_accuracy.nii'), opj(results_dir, 'time-vs-dots', 'snpm_batch', 'SnPM_filtered_masked_accuracy.nii'))