'''
calculations and plots for magnitude difference check and precision judgement check
'''

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

# ------------ File I/O ------------
experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
behavioral_dir = opj(experiment_dir, 'behavioural')
output_dir = abspath('/Users/jachtzehn/Downloads')
appendix_dir = abspath('/Users/jachtzehn/Documents')

# ------------ options ------------
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
conditions = ['time', 'dist', 'dots', 'lumin']
plot = True
exclude_outliers = True
outlier_range = 1.5
bp_pos = [0.85, 1.15, 1.85, 2.15, 2.85, 3.15, 3.85, 4.15]
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


# load data
pse_data_all = pd.read_csv(opj(behavioral_dir, 'pse_data_all_distinct.tsv'), delimiter='\t')
os.system('clear')

if plot:
	fig1, ax = plt.subplots(1, 1, figsize=(15, 8))
	bp_data = []
# 1. manipulation check if small values are perceived as smaller than large values
for i, condition in enumerate(conditions):
	pse_data = pse_data_all[pse_data_all.condition == condition]

	# plot summary of the data
	print('------ Summary of data for condition {} -------'.format(condition))
	print(rp.summary_cont(pse_data.thresh.groupby(pse_data.magnitude)))
	print(rp.summary_cont(pse_data.width.groupby(pse_data.magnitude)))

	# exclude outliers
	if exclude_outliers:

		# based on threshold
		for magn in ['small', 'large']:
			pse_data_magn = pse_data[(pse_data.magnitude == magn)]

			# thresh
			q25, q75 = np.percentile(pse_data_magn.thresh, 25), np.percentile(pse_data_magn.thresh, 75)
			iqr = q75 - q25
			# print('Condition {}, magnitude large, 25 = {:.3f}, 75 = {:.3f}, IQR = {:.4f}'.format(condition, q25, q75, iqr))

			outlier_mask = np.invert(((q25 - iqr * outlier_range <= pse_data_magn.thresh) & (pse_data_magn.thresh <= q75 + iqr * outlier_range)))
			subject_outliers = pse_data_magn.subject[outlier_mask].values

			pse_data = pse_data[~pse_data.subject.isin(subject_outliers)]
			print('Removed participants bases on threshold: {}, condition: {}, magnitude: {}'.format(subject_outliers, condition, magn))

	# test if data is normally distributed
	[res_ks, p_ks] = stats.normaltest(pse_data.thresh[pse_data.magnitude == 'large'])
	print('Condition {}, magnitude: large,  ks-test result: {:.4f}, p={:.6f}'.format(condition, res_ks, p_ks))
	[res_ks, p_ks] = stats.normaltest(pse_data.thresh[pse_data.magnitude == 'small'])
	print('Condition {}, magnitude: small,  ks-test result: {:.4f}, p={:.6f}'.format(condition, res_ks, p_ks))

	[res_t, p_t] = stats.ttest_ind(pse_data.thresh[pse_data.magnitude == 'large'], pse_data.thresh[pse_data.magnitude == 'small'])
	[res_w, p_w] = stats.wilcoxon(pse_data.thresh[pse_data.magnitude == 'large'], pse_data.thresh[pse_data.magnitude == 'small'], alternative='greater')

	print('Condition {} t-test result: {:.5f}, p={:.6f}'.format(condition, res_t, p_t))
	print('Condition {} w-test result: {:.5f}, p={:.6f}'.format(condition, res_w, p_w))
	print('\n\n')

	if plot:
		bp_data.append(pse_data_all.thresh[(pse_data_all.magnitude == 'large') & (pse_data_all.condition == condition)].values)
		bp_data.append(pse_data_all.thresh[(pse_data_all.magnitude == 'small') & (pse_data_all.condition == condition)].values)

# plot
if plot:
	x = range(0, 5)
	boxprops = dict(linewidth=0, color='black')
	medianprops = dict(linewidth=2, color='white')
	whiskerprops = dict(linewidth=2, color='slategray')

	bp = ax.boxplot(bp_data, notch=False, vert=True, whis=1.5, positions=bp_pos, widths=bp_width,
	                showfliers=False, zorder=2, showcaps=False, boxprops=boxprops, medianprops=medianprops,
	                whiskerprops=whiskerprops)

	ax.set_xticks([1, 2, 3, 4])
	ax.set_xticklabels(['time', 'space', 'numerosity', 'control'])
	ax.set_ylabel('Normalised PSE')

	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)

	for i in range(len(bp_data)):
		x = np.linspace(bp_pos[i], bp_pos[i], len(bp_data[i]))
		# for j in range(len(x)):
		# 	x[j] = x[j] + random.random() * (bp_width/2.5) * random.choice((-1, 1))
		y = bp_data[i]

		if (i % 2) == 0:
			color = 'indianred'
			label = 'Large magnitude'
		else:
			color = 'GoldenRod'
			label = 'Small magnitude'

		if i == 0 or i == 1:
			ax.scatter(x, y, alpha=0.65, color=color, zorder=3, linewidths=0, label=label, s=50)
		else:
			ax.scatter(x, y, alpha=0.65, color=color, linewidths=0, zorder=3, s=50)

	# significance
	for i in range(1, 5):
		ax.annotate('***', (i, 2.8), ha='center', fontsize=22, color='slategray')
		ax.plot([i - 0.15, i + 0.15], [2.8, 2.8], color='slategray', lw=1)
		ax.plot([i - 0.15, i - 0.15], [2.78, 2.8], color='slategray', lw=1)
		ax.plot([i + 0.15, i + 0.15], [2.78, 2.8], color='slategray', lw=1)

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

		boxPolygon = Polygon(boxCoords, facecolor='slategray')
		ax.add_patch(boxPolygon)

	plt.legend(loc='lower left', ncol=1, fontsize=16)
	plt.tight_layout()
	#plt.savefig(opj(output_dir, 'thresholds-magnitudes.pdf'), format='pdf', dpi=300, fig=fig1)
	plt.close()

# -------------------

if plot:
	fig2, ax = plt.subplots(2, 3, figsize=(14, 8), sharex=True, sharey=True)
# 2. check whether precision correlates between conditions
p_values = []
for condition_pair in itertools.combinations(conditions, 2):

	if plot:
		if condition_pair[0] == 'time' and condition_pair[1] == 'dist':
			ax_row = 0
			ax_column = 0
		elif condition_pair[0] == 'time' and condition_pair[1] == 'dots':
			ax_row = 0
			ax_column = 1
		elif condition_pair[0] == 'dist' and condition_pair[1] == 'dots':
			ax_row = 0
			ax_column = 2
		elif condition_pair[0] == 'time' and condition_pair[1] == 'lumin':
			ax_row = 1
			ax_column = 0
		elif condition_pair[0] == 'dist' and condition_pair[1] == 'lumin':
			ax_row = 1
			ax_column = 1
		elif condition_pair[0] == 'dots' and condition_pair[1] == 'lumin':
			ax_row = 1
			ax_column = 2

	pse_data_all = pd.read_csv(opj(behavioral_dir, 'pse_data_all_distinct.tsv'), delimiter='\t')
	pse_data_width = pse_data_all

	if exclude_outliers:

		for c in condition_pair:
			pse_data_c = pse_data_all[pse_data_all.condition == c]

			q25, q75 = np.percentile(pse_data_c.width, 25), np.percentile(pse_data_c.width, 75)
			iqr = q75 - q25

			outlier_mask = np.invert(((q25 - iqr * outlier_range <= pse_data_c.width) & (pse_data_c.width <= q75 + iqr * outlier_range)))
			subject_outliers = pse_data_c.subject[outlier_mask].values

			pse_data_width = pse_data_width[~pse_data_width.subject.isin(subject_outliers)]

			print('Removed participants bases on width: {},'.format(subject_outliers))
		c1 = pse_data_width[pse_data_all.condition == condition_pair[0]]
		c2 = pse_data_width[pse_data_all.condition == condition_pair[1]]

	else:
		c1 = pse_data_all[pse_data_all.condition == condition_pair[0]]
		c2 = pse_data_all[pse_data_all.condition == condition_pair[1]]

	mean_c1 = []
	mean_c2 = []

	for subj in c1.subject.unique():
		mean_c1.append(np.mean([c1.width[(c1.subject == subj) & (c1.magnitude == 'small')].values[0],
		                       c1.width[(c1.subject == subj) & (c1.magnitude == 'large')].values[0]]))
		mean_c2.append(np.mean([c2.width[(c2.subject == subj) & (c2.magnitude == 'small')].values[0],
		                       c2.width[(c2.subject == subj) & (c2.magnitude == 'large')].values[0]]))

	# [res_ks1, p_ks1] = stats.normaltest(mean_c1)
	# [res_ks2, p_ks2] = stats.normaltest(mean_c2)
	# print('Condition {}, ks-test result: {:.4f}, p={:.6f}'.format(condition_pair[0], res_ks1, p_ks1))
	# print('Condition {}, ks-test result: {:.4f}, p={:.6f}'.format(condition_pair[1], res_ks2, p_ks2))

	slope, intercept, r_value, p_value, std_err = stats.linregress(stats.rankdata(mean_c1), stats.rankdata(mean_c2))

	if plot:
		ax[ax_row, ax_column].scatter(stats.rankdata(mean_c1), stats.rankdata(mean_c2), alpha=0.65, color='slategray')

		condition_pair_plot = list(condition_pair)

		for i, condition in enumerate(condition_pair_plot):
			if condition == 'dist':
				condition_pair_plot[i] = 'space'
			if condition == 'dots':
				condition_pair_plot[i] = 'numerosity'
			if condition == 'lumin':
				condition_pair_plot[i] = 'control'
		x = np.linspace(0, np.max(stats.rankdata(mean_c1)), 100)
		ax[ax_row, ax_column].plot(x, intercept + slope * x, 'r', label='fitted line', color='indianred')
		ax[ax_row, ax_column].spines['top'].set_visible(False)
		ax[ax_row, ax_column].spines['right'].set_visible(False)
		ax[ax_row, ax_column].spines['top'].set_visible(False)
		ax[ax_row, ax_column].spines['right'].set_visible(False)
		ax[ax_row, ax_column].set_title('{} \& {}'.format(condition_pair_plot[0], condition_pair_plot[1]))

		ax[ax_row, ax_column].set_xlim([0, 22])
		ax[ax_row, ax_column].set_ylim([0, 22])

		if ax_row == 1:
			ax[ax_row, ax_column].set_xlabel('rank(x)')
		if ax_column == 0:
			ax[ax_row, ax_column].set_ylabel('rank(y)')

	[spear_r, spear_p] = stats.spearmanr(mean_c1, mean_c2)
	if plot:
		ax[ax_row, ax_column].annotate('r={:.2f}, p={:.3f}'.format(spear_r, spear_p), [0.5, 20], fontsize=16, color='black')

	print('Condition pair: {} vs. {}, {:.3f} & {:.4f}, slope: {:.4f}'.format(condition_pair[0], condition_pair[1], spear_r, spear_p, slope))
	p_values.append(spear_p)

if plot:
	plt.tight_layout(h_pad=5, w_pad=-0.5)
	plt.savefig(opj(appendix_dir, 'spearman-corr.pdf'), format='pdf', dpi=300, fig=fig2)
	plt.close()

result_multitest = multipletests(p_values, alpha=0.05, method='hs')
print(result_multitest)
