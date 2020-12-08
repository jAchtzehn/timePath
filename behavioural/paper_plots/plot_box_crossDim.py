from os.path import join as opj
from os.path import abspath
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

# ------------ options ------------

# update font sizes
matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams['ytick.labelsize'] = 18
matplotlib.rcParams['xtick.labelsize'] = 18

styles = {'time': ['v-', 'SkyBlue', 'full'],
          'dist': ['o:', 'IndianRed', 'full'],
          'dots': ['x--', 'GoldenRod', 'full'],
          'lumin': ['D-.', 'DimGray', 'full']}

# ------------ File I/O ------------
output_filename = 'pse_data_cross_dim_individual_norm'
experiment_dir = abspath('/home/achtzehnj/data/timePath/derivatives')
behavioral_dir = opj(experiment_dir, 'behavioural')
output_dir = abspath('/home/achtzehnj/Desktop')

# ---- code ----
data = pd.read_csv(opj(behavioral_dir, output_filename + '.tsv'), delimiter='\t')

boxData = []
combinations = [['time', 'dist'], ['dist', 'time'], ['time', 'dots'], ['dots', 'time'], ['dots', 'dist'], ['dist', 'dots']]
# aggregate data
for cond_combination in combinations:
	scatterData = data.pse_diff[(data.condition_rel == cond_combination[0]) & (data.condition_irrel == cond_combination[1])].values.tolist()
	sd = np.std(scatterData)
	upper_boundary = np.mean(scatterData) + 3 * sd
	lower_boundary = np.mean(scatterData) - 3 * sd
	scatterData_cleaned = [x for x in scatterData if x < upper_boundary and x > lower_boundary]
	boxData.append(scatterData_cleaned)

fig1, ax = plt.subplots(1, 1, figsize=(12, 9))
bp_pos = [0.7, 0.95, 1.4, 1.65, 2.1, 2.35]
bp_width = 0.175
x = range(0, 4)

boxprops = dict(linewidth=0, color='black')
medianprops = dict(linewidth=2, color='white')
whiskerprops = dict(linewidth=2, color='slategray')

bp = ax.boxplot(boxData, notch=False, vert=True, whis=1.5, positions=bp_pos, widths=bp_width,
                showfliers=False, zorder=2, showcaps=False, boxprops=boxprops, medianprops=medianprops,
                whiskerprops=whiskerprops)

ax.set_xticks(bp_pos)
#ax.set_xticklabels(['time on space', 'space on time', 'time on numerosity', 'numerosity on time', 'numerosity on space', 'space on numerosity'])
ax.set_xticklabels(['', '', '', '', '', ''])
# plt.tick_params(
#     axis='x',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     bottom=False,      # ticks along the bottom edge are off
#     top=False,         # ticks along the top edge are off
#     labelbottom=False) # labels along the bottom edge are off
ax.set_ylabel('PSE difference')
ax.set_ylim([-0.7, 1])
ax.set_xlim([0.5, 2.55])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

for i in range(len(boxData)):
	x = np.linspace(bp_pos[i], bp_pos[i], len(boxData[i]))
	# for j in range(len(x)):
	# 	x[j] = x[j] + random.random() * (bp_width/2.5) * random.choice((-1, 1))
	y = boxData[i]


	ax.scatter(x, y, alpha=0.65, color='indianred', zorder=3, linewidths=0, s=50)


# shade bp area
numBoxes = len(boxData)
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

	annotations = combinations[i]
	for j, condition in enumerate(annotations):
		if condition == 'dots':
			annotations[j] = 'numerosity'
		if condition == 'dist':
			annotations[j] = 'space'

	if annotations[1] == 'numerosity':
		ax.text(bp_pos[i], -0.95, '{}\non {}'.format(annotations[1], annotations[0]), fontsize=18, color='black',
				horizontalalignment='center', rotation=-45, verticalalignment='center')
	else:
		ax.text(bp_pos[i], -0.95, '{} on\n{}'.format(annotations[1], annotations[0]), fontsize=18, color='black',
				horizontalalignment='center', rotation=-45, verticalalignment='center')

	if i == 0:
		ax.text(bp_pos[i], np.max(boxData[i]) + 0.075, '***', fontsize=18, color='black', horizontalalignment='center')
	elif i == 1:
		ax.text(bp_pos[i], np.max(boxData[i]) + 0.075, '**', fontsize=18, color='black', horizontalalignment='center')
	else:
		ax.text(bp_pos[i], np.max(boxData[i]) + 0.075, 'n.s.', fontsize=18, color='black', horizontalalignment='center')

	ax.plot([0, 4], [0, 0], color='gray', lw=0.1, linestyle='--')
plt.tight_layout(6)
plt.savefig(opj(output_dir, 'pse_diff.pdf'), format='pdf', dpi=300, fig=fig1)
plt.close()