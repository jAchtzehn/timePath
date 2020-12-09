from os.path import join as opj
from os.path import abspath
from sys import platform
import json
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/media/sf_data/fMRI/timePath/')

parameters_to_plot = ['tsnr', 'fd_mean', 'fd_perc']

# read in data of all participants
qdata = []
for subj in range(1, 26):
	qdata_subj = []
	for ses in range(1, 3):
		for run in range(1, 5):
			
			# read in functional quality parameter file
			jsonFile = opj(experiment_dir, 'mriqc', 'sub-' + str(subj).zfill(2), 'ses-' + str(ses).zfill(2),
									  'func', 'sub-' + str(subj).zfill(2) + '_ses-' + str(ses).zfill(2)
									  + '_task-class_run-' + str(run).zfill(2) + '_bold' + '.json')
			with open(jsonFile, 'r') as datafile:
				qdata_run = json.load(datafile)
			
			qdata_subj.append(qdata_run)
			
			
	qdata.append(qdata_subj)

for param in parameters_to_plot:
	
	# prepare data for boxplot
	boxData = []
	for subj in range(25):
		boxData_runs = []
		for run in range(8):
			boxData_runs.append(qdata[subj][run][param])
		
		boxData.append(boxData_runs)
	
	# create plot
	fig, ax1 = plt.subplots(figsize=(25, 12))
	bp = ax1.boxplot(boxData, notch=False, vert=True, whis=1.5)
	
	plt.tight_layout(5)
	numBoxes = len(boxData)
	medians = list(range(numBoxes))
	
	# lim
	maxVal = np.max(boxData)
	minVal = np.min(boxData)
	top = maxVal + 0.05 * maxVal
	bottom = minVal - 0.05 * minVal
	ax1.set_ylim(bottom, top)
	
	plt.setp(bp['boxes'], color='black')
	plt.setp(bp['medians'], color='#2ecc71')
	plt.setp(bp['whiskers'], color='black')
	plt.setp(bp['fliers'], color='red', marker='o')

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
	ax1.set_xlabel('Subjects')
	ax1.set_ylabel(param)
	
	# xtick above
	pos = np.arange(numBoxes) + 1
	upperLabels = [str(np.round(s, 2)) for s in medians]
	weights = ['bold', 'semibold']
	for tick, label in zip(range(numBoxes), ax1.get_xticklabels()):
		k = tick % 2
		ax1.text(pos[tick], top, upperLabels[tick],
			 horizontalalignment='center', fontsize=12, weight=weights[k],
		         color=cmap(normalize(medians[tick])))
	
	# save
	plt.savefig(opj(experiment_dir, 'mriqc_plots', param + '.png'))

