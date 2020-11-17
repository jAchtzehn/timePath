from os.path import join as opj
from os.path import abspath
from sys import platform
import _pickle as cPickle
import pandas as pd

# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/media/sf_data/fMRI/timePath/')

nilearn_dir = opj(experiment_dir, 'nilearn')

# ------------ options ------------
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]                       # subjects to plot
data_suffix = 'solver-' + 'linRegr'
file_suffix = 'solver-' + 'linRegr'
space = 'MNI152NLin2009cAsym'

pd_subjID = []
pd_roi = []
pd_cond1 = []
pd_cond2 = []
pd_cv_accuracy = []
pd_p = []

for subj in subjects:

	subj_str = 'sub-' + str(subj).zfill(2)
	
	with open(opj(nilearn_dir, subj_str, 'space-' + space, 'results', subj_str + '_' + data_suffix + '_classData.pkl'), 'rb') as f:
			classData = cPickle.load(f)
	
	roi_names = sorted(classData.keys())                    # names of the ROIs
	conditions = sorted(classData[roi_names[0]].keys(), key=lambda v: ('lumin' in v, v))
	
	data_cv_score = {'roi': roi_names}
	data_p = {'roi': roi_names}

	for cond in conditions:
		data_cv_score[cond] = []
		data_p[cond] = []

		score_index = cond.find('-')
		cond1 = cond[:score_index]
		cond2 = cond[score_index + 1:]
		
		for roi in roi_names:
			data_cv_score[cond].append(classData[roi][cond][0])
			data_p[cond].append(classData[roi][cond][3])

			pd_subjID.append(subj)
			pd_roi.append(roi)
			pd_cond1.append(cond1)
			pd_cond2.append(cond2)
			pd_cv_accuracy.append(classData[roi][cond][0])
			pd_p.append(classData[roi][cond][3])

	df_cv_score = pd.DataFrame.from_dict(data_cv_score)
	df_p = pd.DataFrame.from_dict(data_p)
	
	df_cv_score = df_cv_score.reindex(sorted(df_cv_score.columns, key=lambda v: ('roi' not in v, v)), axis=1)
	df_p = df_p.reindex(sorted(df_p.columns, key=lambda v: ('roi' not in v, v)), axis=1)
	
	df_cv_score.to_csv(opj(nilearn_dir, subj_str, 'space-' + space, 'results', subj_str + '_' + data_suffix + '_cv_scores.tsv'), sep='\t', index=False)
	df_p.to_csv(opj(nilearn_dir, subj_str, 'space-' + space, 'results', subj_str + '_' + data_suffix + '_p.tsv'), sep='\t', index=False)

df = pd.DataFrame({'subject': pd_subjID,
                   'roi': pd_roi,
                   'condition_1': pd_cond1,
                   'condition_2': pd_cond2,
                   'accuracy': pd_cv_accuracy,
                   'p': pd_p})
df.to_csv(opj(nilearn_dir, 'group_results_space-' + space, 'decoding_accuracies_' + file_suffix + '.csv'), sep='\t',
          columns=['subject', 'roi', 'condition_1', 'condition_2', 'accuracy', 'p'], index=False)

