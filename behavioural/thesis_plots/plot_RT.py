from os.path import join as opj
from os.path import abspath
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

# from utilFunctions import readSubjData
import researchpy as rp
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.stats import multicomp as mc
import pylab as pl

# ------------ options ------------
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
conditions = ['time', 'dist', 'dots', 'lumin']
plot_results = True
fillbetween_yerr = False
clean_up_data = True            # remove RTs lower than 0.3s and higher than 2s and remove for participants 4, 5, 6 lumin in the first session
boundaries_RT = [0.3, 2.1]

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
RT_data_all = pd.read_csv(opj(behavioral_dir, 'RT_data_all_subj.tsv'), delimiter='\t')

# clean up data
if clean_up_data:
	# delete RTs below 300 ms
	RT_data_all = RT_data_all[(RT_data_all.RT >= boundaries_RT[0]) & (RT_data_all.RT <= boundaries_RT[1])]

	# delete session 1 lumin trials for subjects 4, 5 and 6
	RT_data_all = RT_data_all[np.invert((RT_data_all.subject == 4) & (RT_data_all.session == 1) & (RT_data_all.trial_type == 'lumin'))]
	RT_data_all = RT_data_all[np.invert((RT_data_all.subject == 5) & (RT_data_all.session == 1) & (RT_data_all.trial_type == 'lumin'))]
	RT_data_all = RT_data_all[np.invert((RT_data_all.subject == 6) & (RT_data_all.session == 1) & (RT_data_all.trial_type == 'lumin'))]

	# delete participant 15
	RT_data_all = RT_data_all[RT_data_all.subject != 15]

# 1. plot RTs across all runs for each of the four conditions

# 1.plotting
if plot_results:

	fig, ax1 = plt.subplots(figsize=(15, 8))
	plt.tight_layout(5)
	for condition in conditions:

		mean_rt_per_run = []
		yerr_rt_per_run = []
		for run in range(8):

			mean_rt_per_subj = []
			for subj in subjects:
				rt_run_subj = RT_data_all['RT'][(RT_data_all['trial_type'] == condition) & (RT_data_all['run'] == run + 1) & (RT_data_all['subject'] == subj)]
				if not rt_run_subj.empty:
					mean_rt_per_subj.append(np.mean(rt_run_subj))

			mean_rt_per_run.append(np.mean(mean_rt_per_subj))
			yerr_rt_per_run.append(1.96 * (np.std(mean_rt_per_subj) / np.sqrt(len(mean_rt_per_subj))))

		# jitter x-axis a bit
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

		ep = ax1.errorbar(x, mean_rt_per_run, yerr=yerr_rt_per_run, label=label_cond, alpha=1,
		                  color=styles[condition][1],  fmt=styles[condition][0],
		                  linewidth=lw, ms=ms, fillstyle=styles[condition][2],
		                  capsize=0, capthick=0, elinewidth=0)

		if fillbetween_yerr:
			ax1.fill_between(range(1, 9), mean_rt_per_run - yerr_rt_per_run, mean_rt_per_run + yerr_rt_per_run, alpha=0.25)

	ax1.set_xlabel('Runs', fontsize=20)
	ax1.set_ylabel('Time [s]', fontsize=20)
	ax1.set_ylim([0.9, 1.3])
	ax1.grid(True, axis='y', linestyle='--')
	plt.legend(loc='lower right', ncol=4, fontsize=16)
	ax1.spines['top'].set_visible(False)
	ax1.spines['right'].set_visible(False)

	plt.savefig(opj(output_dir, 'rt-runs.pdf'), format='pdf', dpi=300)


# 2.statistical analysis

# 2.1 ANOVA to test if there is a difference between task types (i.e. difference between the RT of time, dist, dots and lumin?)
#anova_rt_trial = stats.f_oneway(RT_data_all['RT'][RT_data_all['trial_type'] == 'time'], RT_data_all['RT'][RT_data_all['trial_type'] == 'dist'])
#print(anova_rt_trial)

#anova_rt_run = stats.f_oneway(RT_data_all['RT'][RT_data_all['run'] == 1], RT_data_all['RT'][RT_data_all['run'] == 8])
#print(anova_rt_run)

#RT_data_all.boxplot('RT', by='trial_type')
#plt.show()

print('------ Summary of data -------')
print(rp.summary_cont(RT_data_all['RT'].groupby(RT_data_all['trial_type'])))

print('\n\n------- Anova LM -------\n')
anova_ols_rt_trial = smf.ols('RT ~ trial_type', data=RT_data_all).fit()
print(anova_ols_rt_trial.summary())

print('\n\n------- Anova table -------\n')

aov_table = sm.stats.anova_lm(anova_ols_rt_trial, typ=2)
print(aov_table)

print('\n\n------- Multiple Comparisons + Tukeys HSD -------\n')
mc_trial_type = mc.MultiComparison(RT_data_all['RT'], RT_data_all['trial_type'])
mc_results = mc_trial_type.tukeyhsd()
print(mc_results)

