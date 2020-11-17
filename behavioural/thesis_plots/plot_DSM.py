from os.path import join as opj
from os.path import abspath
import numpy as np
import os
import pandas as pd
import researchpy as rp
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import itertools
import matplotlib
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib.patches import Polygon
import _pickle as cPickle
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ------------ File I/O ------------
experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
rsa_dir = opj(experiment_dir, 'rsa', 'group_results_space-MNI152NLin2009cAsym')

# ------------ options ------------
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
conditions = ['time', 'dist', 'dots']
condition_pairs = ['time-dist', 'time-dots', 'dist-dots']
rois_plot = ['ips', 'ifg', 'V5', 'cuneus', 'precuneus']
rois_sign_calc = ['ips_lh', 'ips_rh', 'ifg_lh', 'ifg_rh', 'V5_lh', 'V5_rh', 'cuneus_lh', 'cuneus_rh', 'precuneus_lh', 'precuneus_rh']
exclude_outliers = True
outlier_range = 1.5
bp_pos = [0.85, 1.15, 1.85, 2.15, 2.85, 3.15]
bp_width = 0.175

matplotlib.rc('text', usetex=True)
pl.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',    # helvetica font
    r'\usepackage{sansmath}',   # math-font matching  helvetica
    r'\sansmath'                # actually tell tex to use it!
    r'\usepackage{siunitx}',    # micro symbols
    r'\sisetup{detect-all}',    # force siunitx to use the fonts
]
matplotlib.rcParams['ytick.labelsize'] = 18
matplotlib.rcParams['xtick.labelsize'] = 18
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['axes.labelsize'] = 18
matplotlib.rcParams['axes.labelpad'] = 12
matplotlib.rcParams['axes.titlesize'] = 22
matplotlib.rcParams['axes.titlepad'] = 10


pklfile = opj(rsa_dir, 'RDM_space-MNI152NLin2009cAsym_searchlight_mean.pkl')

with open(pklfile, 'rb') as f:
	[rdm_all_subj, rdm_cumm_subj] = cPickle.load(f, encoding="latin1")

for roi in rois_sign_calc:
	rdm_data = rdm_all_subj[roi]

	print('-------- ROI: {} -------'.format(roi))

	rdm_vals = {}
	for condition_pair in condition_pairs:

		rdm_temp = []
		if condition_pair == 'time-dist':
			for i in range(len(subjects)):
				rdm_temp.append(rdm_data[i][2, 0])
			rdm_vals[condition_pair] = rdm_temp

		elif condition_pair == 'dist-dots':
			for i in range(len(subjects)):
				rdm_temp.append(rdm_data[i][0, 1])
			rdm_vals[condition_pair] = rdm_temp

		elif condition_pair == 'time-dots':
			for i in range(len(subjects)):
				rdm_temp.append(rdm_data[i][2, 1])
			rdm_vals[condition_pair] = rdm_temp

	p_vals = []
	for pair in itertools.combinations(condition_pairs, 2):

		[res_w, p_w] = stats.wilcoxon(rdm_vals[pair[0]], rdm_vals[pair[1]], alternative='two-sided')

		mean_diff = np.abs(np.median(rdm_vals[pair[0]]) - np.median(rdm_vals[pair[1]]))

		if p_w < 0.05:
			print('{} vs {}, median_1 = {:.3f}, median_2 = {:.3f}, median_diff = {:.3f}, w = {:.3f}, p = {:.4f} *****'.format(
				pair[0], pair[1], np.mean(rdm_vals[pair[0]]), np.mean(rdm_vals[pair[1]]), mean_diff, res_w, p_w))
		else:
			print('{} vs {}, median_1 = {:.3f}, median_2 = {:.3f}, median_diff = {:.3f}, w = {:.3f}, p = {:.4f}'.format(
				pair[0], pair[1], np.mean(rdm_vals[pair[0]]), np.mean(rdm_vals[pair[1]]), mean_diff, res_w, p_w))

		p_vals.append(p_w)

	result_multitest = multipletests(p_vals, alpha=0.05, method='hs')
	print(result_multitest[0])

# plotting each ROI
for roi in rois_plot:
	fig1, ax = plt.subplots(1, 1, figsize=(15, 8))
	x = range(0, 4)

	ax.set_ylim([0.5, 2])
	bp_data = []
	for condition_pair in condition_pairs:

		for hemi in ['lh', 'rh']:
			rdm_temp = []

			rdm_data = rdm_all_subj[roi + '_' + hemi]

			if condition_pair == 'time-dist':
				for i in range(len(subjects)):
					rdm_temp.append(rdm_data[i][2, 0])
			elif condition_pair == 'dist-dots':
				for i in range(len(subjects)):
					rdm_temp.append(rdm_data[i][0, 1])
				rdm_vals[condition_pair] = rdm_temp

			elif condition_pair == 'time-dots':
				for i in range(len(subjects)):
					rdm_temp.append(rdm_data[i][2, 1])
				rdm_vals[condition_pair] = rdm_temp

			bp_data.append(rdm_temp)

	boxprops = dict(linewidth=0, color='black')
	medianprops = dict(linewidth=2, color='white')
	whiskerprops = dict(linewidth=2, color='slategray')

	bp = ax.boxplot(bp_data, notch=False, vert=True, whis=1.5, positions=bp_pos, widths=bp_width,
	                showfliers=False, zorder=2, showcaps=False, boxprops=boxprops, medianprops=medianprops,
	                whiskerprops=whiskerprops)

	ax.set_xticks([1, 2, 3])
	ax.set_xticklabels(['time-space', 'time-numerosity', 'space-numerosity'])
	ax.set_ylabel('Dissimilarity [a.u.]')

	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)

	for i in range(len(bp_data)):
		x = np.linspace(bp_pos[i], bp_pos[i], len(bp_data[i]))
		# for j in range(len(x)):
		# 	x[j] = x[j] + random.random() * (bp_width/2.5) * random.choice((-1, 1))
		y = bp_data[i]

		if (i % 2) == 0:
			color = 'indianred'
			label = 'Left hemisphere'
		else:
			color = 'GoldenRod'
			label = 'Right hemisphere'

		if i == 0 or i == 1:
			ax.scatter(x, y, alpha=0.65, color=color, zorder=3, linewidths=0, label=label, s=50)
		else:
			ax.scatter(x, y, alpha=0.65, color=color, linewidths=0, zorder=3, s=50)

	# shade bp area
	numBoxes = len(bp_data)
	medians = list(range(numBoxes))
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
		t = med.get_xdata()[0]
		#ax.annotate('{:.3f}'.format(med.get_ydata()[0]), (t, med.get_ydata()[0]), ha='center', fontsize=22, color='black')
		boxPolygon = Polygon(boxCoords, facecolor='slategray')
		ax.add_patch(boxPolygon)

	plt.legend(loc='lower right', ncol=1, fontsize=16)
	plt.tight_layout()
	plt.savefig(opj(rsa_dir, roi + '_dsm.pdf'), format='pdf', dpi=300, fig=fig1)
	plt.close()


# plot matrices


fig, ax = plt.subplots(2, 5, figsize=(35, 12))

plt.tight_layout(5)
fig.subplots_adjust(hspace=.35, wspace=.55)
# fig.suptitle('Mean RDMs for most stable patterns', fontsize=20)

# find max
max_vals = []
for roi in rois_sign_calc:
	max_vals.append(np.max(rdm_cumm_subj[roi] / len(subjects)))

for roi_enum, roi in enumerate(rois_sign_calc):
	rdm_cumm_subj[roi] = rdm_cumm_subj[roi] / len(subjects)  # compute mean

	if (len(rois_sign_calc) % 2) == 0:
		if roi_enum < len(rois_sign_calc) / 2:
			subplot_row = 0
			subplot_col = int(roi_enum)
		else:
			subplot_row = 1
			subplot_col = int(roi_enum - len(rois_sign_calc) / 2)
	else:
		subplot_row = 0
		subplot_col = int(roi_enum)

	# plotting
	rdm_plot = ax[subplot_row, subplot_col].imshow(rdm_cumm_subj[roi], interpolation='nearest')
	div = make_axes_locatable(ax[subplot_row, subplot_col])
	cax = div.append_axes("right", size="5%", pad=0.2)
	cbar = plt.colorbar(rdm_plot, cax=cax)
	rdm_plot.set_clim(vmin=0, vmax=np.max(max_vals))
	ax[subplot_row, subplot_col].set_xticklabels(['space', 'numerosity', 'time'], fontdict=None, minor=False,
	                                             rotation=-45)
	ax[subplot_row, subplot_col].set_xticks(range(len(rdm_cumm_subj[roi])))
	ax[subplot_row, subplot_col].set_yticklabels(['space', 'numerosity', 'time'], fontdict=None, minor=False)
	ax[subplot_row, subplot_col].set_yticks(range(len(rdm_cumm_subj[roi])))
	ax[subplot_row, subplot_col].set_title('{}'.format(roi.replace('_', '-')))

	for i in range(len(rdm_cumm_subj[roi])):
		for j in range(len(rdm_cumm_subj[roi])):
			if i != j:
				ax[subplot_row, subplot_col].text(j, i, '{:.2f}'.format(rdm_cumm_subj[roi][i, j]),
				                                  ha='center', va='center', color='w', fontsize=22)

plt.savefig(opj(rsa_dir, 'matrix' + '_dsm.pdf'), dpi=300)
plt.close()