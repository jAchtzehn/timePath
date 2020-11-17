from os.path import join as opj
from os.path import abspath
import json
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
import pylab as pl
import json
import pandas as pd

# ------------ File I/O ------------
mriqc_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/mriqc')
fmriprep_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/fmriprep')
output_dir = abspath('/Users/jachtzehn/Documents/Medizin/thesis/figures/results/fmri/quality')

# ------------ options ------------
parameters_to_plot = ['tsnr', 'fd_mean', 'fd_perc', 'dvars_std', 'fber', 'efc', 'fwhm_avg', 'aqi', 'aor', 'gsr_x', 'gsr_y']
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
outlier_range = 1.5
bp_width = 0.5
figsize = (24, 10)
boxprops = dict(linewidth=0, color='slategray')
medianprops = dict(linewidth=2, color='white')
whiskerprops = dict(linewidth=2, color='slategray')

plot_fd_subjects = False

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

# read in data of all participants
qdata = []
for subj in range(1, len(subjects) + 1):
	qdata_subj = []
	for ses in range(1, 3):
		for run in range(1, 5):
			# read in functional quality parameter file
			jsonFile = opj(mriqc_dir, 'sub-' + str(subj).zfill(2), 'ses-' + str(ses).zfill(2),
			               'func', 'sub-' + str(subj).zfill(2) + '_ses-' + str(ses).zfill(2)
			               + '_task-class_run-' + str(run).zfill(2) + '_bold' + '.json')
			with open(jsonFile, 'r') as datafile:
				qdata_run = json.load(datafile)

			qdata_subj.append(qdata_run)

	qdata.append(qdata_subj)

# quality parameters for group level
for param in parameters_to_plot:

	# prepare data for boxplot
	boxData = []
	boxData_mean = []
	for subj in range(len(subjects)):
		boxData_runs = []
		for run in range(8):
			boxData_runs.append(qdata[subj][run][param])

		boxData.append(boxData_runs)
		boxData_mean.append(np.mean(boxData_runs))

	# create plot
	fig, ax1 = plt.subplots(figsize=figsize)
	bp = ax1.boxplot(boxData, notch=False, vert=True, whis=1.5, widths=bp_width,
	                showfliers=False, zorder=4, showcaps=False, boxprops=boxprops, medianprops=medianprops,
	                whiskerprops=whiskerprops)
	ax1.spines['top'].set_visible(False)
	ax1.spines['right'].set_visible(False)

	# lim
	maxVal = np.max(boxData)
	minVal = np.min(boxData)
	top = maxVal + 0.05 * maxVal
	bottom = minVal - 0.05 * minVal
	#ax1.set_ylim(bottom, top)
	ax1.set_xlim(0.5, 25.5)
	# grid
	ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
	               alpha=0.5)
	ax1.set_axisbelow(True)

	# label
	ax1.set_xlabel('Participants')

	if param == 'fd_mean':
		ylabel = 'Mean FD [mm]'
	elif param == 'fd_perc':
		ylabel = 'Percentage of FD above 1 mm'
	elif param == 'tsnr':
		ylabel = 'tSNR'
	elif param == 'efc':
		ylabel = 'EFC'
	elif param == 'fber':
		ylabel = 'FBER [a.u.]'
	elif param == 'gsr_x':
		ylabel = 'GSR x-axis'
	elif param == 'gsr_y':
		ylabel = 'GSR y-axis'
	else:
		ylabel = param.replace('_', '-') + ' [a.u.]'

	ax1.set_ylabel(ylabel)

	# scatter of original data
	for i in range(len(boxData)):
		x = np.linspace(i + 1, i + 1, len(boxData[i]))
		y = boxData[i]
		ax1.scatter(x, y, alpha=0.65, color='indianred', zorder=3, linewidths=0, label='', s=50)

	# means and IQRs
	median = np.median(boxData_mean)
	q75, q25 = np.percentile(boxData_mean, 75), np.percentile(boxData_mean, 25)
	iqr = q75 - q25
	upper_boundary = q75 + iqr * outlier_range
	lower_boundary = q25 - iqr * outlier_range
	ax1.plot([0.5, 25.5], [median, median], zorder=2., color='goldenrod')
	# ax1.plot([0.5, 25.5], [upper_boundary, upper_boundary], zorder=1)
	# ax1.plot([0.5, 25.5], [lower_boundary, lower_boundary], zorder=1)
	print('Param: {}, median: {:.4f}'.format(param, median))
	if param == 'fd_perc' or param == 'fd_mean' or param == 'gsr_y':
		lower_boundary = 0

	# add text to upper and lower boundary
	#ax1.annotate('$Q_{75} + IQR \cdot 1.5$', (23, upper_boundary - 0.01 * upper_boundary), fontsize=22, color='slategray')

	ax1.axhspan(lower_boundary, upper_boundary, facecolor='lightgrey', alpha=0.3, zorder=1)
	numBoxes = len(boxData)
	medians = list(range(numBoxes))

	# shade area of bp
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
		ax1.add_patch(boxPolygon)

	fname = param.replace('_', '-')
	plt.savefig(opj(output_dir, fname + '.pdf'), format='pdf', dpi=300)

	plt.close(fig)

if plot_fd_subjects:
	# plot FD for each participant and each run
	for subj in range(1, len(subjects) + 1):
		fig, ax1 = plt.subplots(figsize=figsize)
		for ses in range(1, 3):
			for run in range(1, 5):
				# read in functional quality parameter file
				confounds_filename = opj(fmriprep_dir, 'sub-' + str(subj).zfill(2), 'ses-' + str(ses).zfill(2),
				                         'func', 'sub-' + str(subj).zfill(2) + '_ses-' + str(ses).zfill(2)
				                         + '_task-class_run-' + str(run).zfill(2) + '_bold_confounds' + '.tsv')
				confounds_data = pd.read_csv(confounds_filename, delimiter='\t')
				confounds_data.fillna(0)
				ax1.plot(confounds_data.FramewiseDisplacement)
				ax1.set_ylim([0, 10])

		plt.tight_layout()
		plt.savefig(opj(output_dir, 'sub-' + str(subj).zfill(2) + '-fd.pdf'), format='pdf', dpi=300)
		plt.close(fig)
