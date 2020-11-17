from os.path import join as opj
from os.path import abspath
import pandas as pd
import numpy as np
import matplotlib

import matplotlib.pyplot as plt
import pylab as pl

# ------------ options ------------
subjects = range(1, 26)
conditions = ['time', 'dist', 'dots', 'lumin']
plot_results = True
fillbetween_yerr = False

# plotting options
matplotlib.rc('text', usetex=True)
pl.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',    # helvetica font
    r'\usepackage{sansmath}',   # math-font matching  helvetica
    r'\sansmath'                # actually tell tex to use it!
    r'\usepackage{siunitx}',    # micro symbols
    r'\sisetup{detect-all}',    # force siunitx to use the fonts
]
matplotlib.rcParams['ytick.labelsize'] = 22
matplotlib.rcParams['xtick.labelsize'] = 22
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['axes.labelsize'] = 28
matplotlib.rcParams['axes.labelpad'] = 12
matplotlib.rcParams['axes.titlesize'] = 32
matplotlib.rcParams['axes.titlepad'] = 20

figsize = (16, 10)
tight_pad = 5
lw = 3
ms = 10

styles = {'time': ['v-', 'SkyBlue', 'full'],
          'dist': ['o:', 'IndianRed', 'full'],
          'dots': ['s--', 'GoldenRod', 'full'],
          'lumin': ['D-.', 'DimGray', 'full']}

# ------------ File I/O ------------
experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
behavioral_dir = opj(experiment_dir, 'behavioural')
output_dir = abspath('/Users/jachtzehn/Documents/Medizin/thesis/figures/results/behavioural')

# read in RT data
vas_data_all = pd.read_csv(opj(behavioral_dir, 'vas_data_all_subj.tsv'), delimiter='\t')

if plot_results:

	fig, ax1 = plt.subplots(figsize=(15, 8))
	plt.tight_layout(4)
	for condition in conditions:

		vas = []
		yerr_vas = []

		for run in range(8):
			vas_run = []
			for subj in subjects:
				vas_subj = vas_data_all[(vas_data_all.run == run + 1) & (vas_data_all.subject == subj) & (vas_data_all.trial_type == condition)]
				if not vas_subj.empty:
					vas_run.append(vas_subj.vas.values[0])

			vas.append(np.mean(vas_run))
			yerr_vas.append(1.96 * (np.std(vas_run) / np.sqrt(len(vas_run))))

		x = list(range(1, 9))
		for nr in range(0, 8):
			if condition == 'time':
				x[nr] += .15
			elif condition == 'dist':
				x[nr] -= .15
			elif condition == 'lumin':
				x[nr] += .05
			else:
				x[nr] -= .05

		# rename lumin to control
		if condition == 'lumin':
			label_cond = 'control'
		elif condition == 'dist':
			label_cond = 'space'
		elif condition == 'dots':
			label_cond = 'numerosity'
		else:
			label_cond = condition

		ep = ax1.errorbar(x, vas, yerr=yerr_vas, label=label_cond, alpha=1,
		                  color=styles[condition][1],  fmt=styles[condition][0],
		                  linewidth=lw, ms=ms, fillstyle=styles[condition][2],
		                  capsize=0, capthick=0, elinewidth=0)

		if fillbetween_yerr:
			ax1.fill_between(range(1, 9), vas-yerr_vas, vas+yerr_vas, alpha=0.25)

	ax1.set_xlabel('Runs', fontsize=20)
	ax1.set_ylabel('Difficulty rating', fontsize=20)

	ax1.spines['top'].set_visible(False)
	ax1.spines['right'].set_visible(False)
	ax1.set_ylim([0, 1.])
	ax1.grid(True, axis='y', linestyle='--')
	plt.legend(loc='lower right', ncol=4, fontsize=16)

	plt.savefig(opj(output_dir, 'vas-runs.pdf'), format='pdf', dpi=300)

