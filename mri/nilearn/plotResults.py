from nilearn.plotting import plot_stat_map
from nilearn.image import coord_transform
import numpy as np
import os
from os.path import join as opj
from os.path import abspath
from sys import platform
import _pickle as cPickle
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import ttest_1samp
from statsmodels.stats.multitest import multipletests
from scipy import ndimage
from tqdm import tqdm
from util_functions import create_html_report


# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/media/sf_data/fMRI/timePath/')

nilearn_dir = opj(experiment_dir, 'nilearn')

# ------------ options ------------
subjects = range(1, 26)                        # subjects to plot
# [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]       # all except 15
width = 0.2                                # width of the bars
matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams['ytick.labelsize'] = 20
matplotlib.rcParams['xtick.labelsize'] = 20

data_suffix = 'score-balanced-accuracy'
file_suffix = 'score-balanced-accuracy'

subject_level = False
group_level = True

plot_weights = False
plot_cv = True
write_html = False
correct_multiple_comparisons = True

output_space = ['T1w', 'MNI152NLin2009cAsym']                                    # MNI152NLin2009cAsym
condition_replacements = {'time': 'time', 'dots': 'numerosity', 'dist': 'space', 'lumin': 'control'}


def autolabel(data, data_std, p_data, axis,  xpos='center'):
	"""
    Attach a text label above each bar in *cv_score*, displaying its height.

    *xpos* indicates which side to place the text w.r.t. the center of
    the bar. It can be one of the following {'center', 'right', 'left'}.
    """
	
	xpos = xpos.lower()  # normalize the case of the parameter
	ha = {'center': 'center', 'right': 'left', 'left': 'right'}
	offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

	for cv_score, cv_std, p_value in zip(data, data_std, p_data):
		
		height = cv_score.get_height()
		axis.text(cv_score.get_x() + cv_score.get_width()*offset[xpos], 1.01*height,
		          '{:.2f}'.format(height), ha=ha[xpos], va='bottom', fontsize=20)
		
		additional_voffset = 0
		if cv_std <= 0.1:
			additional_voffset = 0.025
		
		if p_value <= 0.001:
			axis.text(cv_score.get_x() + cv_score.get_width()/2, height + cv_std + additional_voffset,
			          '***', ha='center', va='bottom', fontsize=20)
		elif p_value <= 0.01:
			axis.text(cv_score.get_x() + cv_score.get_width()/2, height + cv_std + additional_voffset,
			          '**', ha='center', va='bottom', fontsize=20)
		elif p_value <= 0.05:
			axis.text(cv_score.get_x() + cv_score.get_width()/2, height + cv_std + additional_voffset,
			          '*', ha='center', va='bottom', fontsize=20)
		else:
			axis.text(cv_score.get_x() + cv_score.get_width()/2, (height + 0.01) + cv_std + additional_voffset,
			          'ns', ha='center', va='bottom', fontsize=20)
	
	
for space in output_space:
	os.system('clear')
	
	print('\nPlotting for subjects: %s, space: %s, file_suffix: %s' % (subjects, space, file_suffix))
	
	# make group results
	if not os.path.exists(opj(nilearn_dir, 'group_results_space-' + space)):
		os.mkdir(opj(nilearn_dir, 'group_results_space-' + space))
	
	if subject_level:
		
		print('\nSubject level\n')
		sbar = tqdm(subjects, leave=False)
		for subj in sbar:
			# -- preparations
			subj_str = 'sub-' + str(subj).zfill(2)      # create the subject string (just for easier handling)
			os.system('mkdir -p %s' % opj(nilearn_dir, subj_str,  'space-' + space, 'results'))        # make sure the results dir exists
			
			with open(opj(nilearn_dir, subj_str, 'space-' + space, 'results', subj_str + '_' + data_suffix + '_classData.pkl'), 'rb') as f:
				classData = cPickle.load(f)
	
			roi_names = sorted(classData.keys())                    # names of the ROIs
			conditions = sorted(classData[roi_names[0]].keys(), key=lambda v: ('lumin' in v, v))
			
			# -- CV accuracy plots
			if plot_cv:
				sbar.set_postfix({'Image file': opj(subj_str, 'space-' + space, 'results', subj_str + '_' + file_suffix + '_cv_scores.png')})
				sbar.set_description('Processing cv scores for subject %s' % str(subj))
				# plotting
				fig, ax = plt.subplots(len(conditions)/2, 2, figsize=(35, 20))
				plt.tight_layout(5)
				fig.subplots_adjust(hspace=0.5)
				
				for cond_enum, condition in enumerate(conditions):
					
					if cond_enum < len(conditions)/2:
						subplot_row = cond_enum
						subplot_col = 0
					else:
						subplot_row = cond_enum - len(conditions)/2
						subplot_col = 1
					
					roi_xticklabels = []                        # list for xtick labels of rois
					lh_data, lh_data_std, lh_p = [], [], []
					rh_data, rh_data_std, rh_p = [], [], []
					
					for roi in roi_names:
						
						roi_xticklabel = roi[0:roi.find('_')]
						if roi_xticklabel not in roi_xticklabels:
							roi_xticklabels.append(roi[0:roi.find('_')])
						
						if 'lh' in roi:
							lh_data.append(classData[roi][condition][0])        # cv score
							lh_data_std.append(classData[roi][condition][1]/2)  # std of cv score
							lh_p.append(classData[roi][condition][3])           # p-value
						else:
							rh_data.append(classData[roi][condition][0])        # cv score
							rh_data_std.append(classData[roi][condition][1]/2)  # std of cv score
							rh_p.append(classData[roi][condition][3])           # p-value
							
					len_samples = classData[roi_names[0]][condition][5]
					
					x_loc = np.arange(len(lh_data))
					lh = ax[subplot_row, subplot_col].bar(x_loc - width / 2, lh_data, width, yerr=lh_data_std, color='SkyBlue', label='lh')
					rh = ax[subplot_row, subplot_col].bar(x_loc + width / 2, rh_data, width, yerr=rh_data_std, color='IndianRed', label='rh')
					ax[subplot_row, subplot_col].set_ylim([0.4, 1.0])
					ax[subplot_row, subplot_col].set_title((condition + ', n = [%s, %s]' % (len_samples[0], len_samples[1])), fontsize=18)
					ax[subplot_row, subplot_col].set_xticks(x_loc)
					ax[subplot_row, subplot_col].set_xticklabels(roi_xticklabels)
					ax[subplot_row, subplot_col].axhline(y=0.5, linewidth=0.5, color='grey', linestyle='--')
					
					autolabel(lh, lh_data_std, lh_p, ax[subplot_row, subplot_col], "left")
					autolabel(rh, rh_data_std, rh_p, ax[subplot_row, subplot_col], "right")
					
				fig.legend(handles=[lh, rh], labels=['left hemisphere', 'right hemisphere'], loc="upper left", borderaxespad=0.1, fontsize=18)
				fig.suptitle('Cross validation accuracy for subject %s; %s' % (subj, file_suffix), fontsize=24)
				
				plt.savefig(opj(nilearn_dir, subj_str, 'space-' + space, 'results', subj_str + '_' + file_suffix + '_cv_scores.png'), dpi=300)
				plt.close()
			
			# -- weight plots
			if plot_weights:
				sbar.set_description('Processing svm weights for subject %s' % str(subj))
	
				cbar = tqdm(conditions, leave=False)
				for condition in cbar:
	
					cbar.set_description('Processing conditions %s' % condition)
					sbar.set_postfix({'Image file': opj(subj_str + '_' + 'condition-' + condition + '_' + file_suffix + '_weights' + '.png')})
	
					fig, ax = plt.subplots(len(roi_names)/2, 2, figsize=(25, 20))
					plt.tight_layout(3)
					fig.subplots_adjust(hspace=0.025)
					fig.subplots_adjust(wspace=0.0125)
	
					row_comp = 0
					rbar = tqdm(roi_names, leave=True)
					for roi_enum, roi in enumerate(rbar):
	
						rbar.set_description('Processing roi %s' % str(roi))
	
						if (roi_enum % 2) == 0:
							subplot_row = roi_enum - row_comp
							subplot_col = 0
							row_comp += 1
						elif roi_enum == 1:
							row_comp = 1
							subplot_row = roi_enum - row_comp
							subplot_col = 1
						else:
							subplot_row = roi_enum - row_comp
							subplot_col = 1
	
						coef_ = classData[roi][condition][4][0]
						coef_squared = np.square(coef_)
						coef_img = classData[roi][condition][4][1].inverse_transform(coef_squared)
						
						# find maximum
						max_data = coef_img.get_data()
						max_pos = ndimage.maximum_position(max_data)
						max_pos_mni = coord_transform(max_pos[0], max_pos[1], max_pos[2], coef_img.affine)
						
						if space == 'T1w':
							anat_data_name = opj(experiment_dir, 'fmriprep', subj_str, 'anat', subj_str + '_T1w_preproc.nii.gz')    # T1 image
						else:
							anat_data_name = opj(experiment_dir, 'fmriprep', subj_str, 'anat', subj_str + '_T1w_space-MNI152NLin2009cAsym_preproc.nii.gz')    # MNI image
							
						#figtitle = ('%s - %s' % (str(conditions_of_interest[0]), str(conditions_of_interest[1])))
						plot_stat_map(coef_img, bg_img=anat_data_name, display_mode='ortho', cut_coords=(max_pos_mni[0], max_pos_mni[1], max_pos_mni[2]),
						              draw_cross=True, cmap='viridis', title=roi, axes=ax[subplot_row, subplot_col])
	
					fig.suptitle('subject: %s, condition: %s' % (subj, condition), fontsize=24)
					plt.savefig(opj(nilearn_dir, subj_str, 'space-' + space, 'results', subj_str + '_' + 'condition-' + condition + '_' + file_suffix + '_weights' + '.png'), dpi=300)
					plt.close()
	
	if write_html and plot_cv:
		# summarize individual plots into a html file
		tqdm.write('Summarizing cv plots into html file...')
		image_list = []
		for subj in subjects:
			subj_str = 'sub-' + str(subj).zfill(2)      # create the subject string (just for easier handling)
			imagePath = opj(nilearn_dir, subj_str, 'space-' + space, 'results', subj_str + '_' + file_suffix + '_cv_scores.png')
			image_list.append(imagePath)
		
		create_html_report(image_list, opj(nilearn_dir, 'group_results_space-' + space, 'group_space-' + space + '_' + file_suffix + '_cv_scores' + '.html'), 'CV scores')
		
	if write_html and plot_weights:
		# summarize individual plots into a html file
		tqdm.write('Summarizing weight plots into html file...')
		image_files = {}
		for subj in subjects:
			subj_str = 'sub-' + str(subj).zfill(2)      # create the subject string (just for easier handling)
	
			with open(opj(nilearn_dir, subj_str, 'space-' + space, 'results', subj_str + '_' + data_suffix + '_classData.pkl'), 'rb') as f:
				classData = cPickle.load(f)
			
			roi_names = sorted(classData.keys())                    # names of the ROIs
			conditions = sorted(classData[roi_names[0]].keys(), key=lambda v: ('lumin' in v, v))
			
			for condition in conditions:
				if condition not in image_files.keys():
					image_files[condition] = [opj(nilearn_dir, subj_str, 'space-' + space, 'results', subj_str + '_' + 'condition-' + condition + '_' + file_suffix + '_weights' + '.png')]
				else:
					image_files[condition].append(opj(nilearn_dir, subj_str, 'space-' + space, 'results', subj_str + '_' + 'condition-' + condition + '_' + file_suffix + '_weights' + '.png'))
		
		for condition in image_files.keys():
			create_html_report(image_files[condition], opj(nilearn_dir, 'group_results_space-' + space, 'group_' + 'condition-' + condition + '_' + file_suffix + '_weights' + '.html'),
			                   'SWM weights for condition: %s' % condition)
			tqdm.write('Created html report: %s' % opj(nilearn_dir, 'group_results_space-' + space, 'group_' + 'condition-' + condition + '_' + file_suffix + '_weights' + '.html'))
			
	if group_level:
		
		print('\nGroup level\n')
		
		# load data of all subjects
		classData = {}
		for subj in subjects:
			subj_str = 'sub-' + str(subj).zfill(2)
			with open(opj(nilearn_dir, subj_str, 'space-' + space, 'results', subj_str + '_' + data_suffix + '_classData.pkl'), 'rb') as f:
				classData_subj = cPickle.load(f)
			
			classData[subj_str] = classData_subj
		
		subj_strs = sorted(classData.keys())
		roi_names = sorted(classData[subj_strs[0]].keys())
		conditions = sorted(classData[subj_strs[0]][roi_names[0]].keys(), key=lambda v: ('lumin' in v, v))
		
		# plotting
		fig, ax = plt.subplots(len(conditions)/2, 2, figsize=(30, 20))
		plt.tight_layout(5)
		fig.subplots_adjust(hspace=0.2, wspace=0.15)
		
		cbar = tqdm(conditions, desc='Conditions', leave=True)
		for cond_enum, condition in enumerate(cbar):
			
			cbar.set_description('Processing condition %s' % str(condition))
			if cond_enum < len(conditions)/2:
				subplot_row = cond_enum
				subplot_col = 0
			else:
				subplot_row = cond_enum - len(conditions)/2
				subplot_col = 1
			
			roi_xticklabels = []                        # list for xtick labels of rois
			
			subjs_mean_data_lh, subjs_ci_data_lh, subjs_p_data_lh = [], [], []
			subjs_mean_data_rh, subjs_ci_data_rh, subjs_p_data_rh = [], [], []
			
			for roi in roi_names:
				roi_xticklabel = roi[0:roi.find('_')]
				if roi_xticklabel not in roi_xticklabels:
					roi_xticklabels.append(roi[0:roi.find('_')])
			
				# now go through all subjects and accumulate the data
				subjs_data_lh, subjs_data_rh = [], []
				for subj in subjects:
					subj_str = 'sub-' + str(subj).zfill(2)
					
					# load all classification data into a list for mean calculation
					if 'lh' in roi:
						subjs_data_lh.append(classData[subj_str][roi][condition][0])
					else:
						subjs_data_rh.append(classData[subj_str][roi][condition][0])
				
				if 'lh' in roi:
					subjs_mean_data_lh.append(np.mean(subjs_data_lh))
					subjs_ci_data_lh.append(1.96 * (np.std(subjs_data_lh) / np.sqrt(len(subjs_data_lh))))
					stats, p = ttest_1samp(subjs_data_lh, 0.5)
					subjs_p_data_lh.append(p)
					
				elif 'rh' in roi:
					subjs_mean_data_rh.append(np.mean(subjs_data_rh))
					subjs_ci_data_rh.append(1.96 * (np.std(subjs_data_rh) / np.sqrt(len(subjs_data_rh))))
					stats, p = ttest_1samp(subjs_data_rh, 0.5)
					subjs_p_data_rh.append(p)
					
				else:
					print('--> WARNING, script has not been programmed for ROIS that are not rh or lh!!!')
			
			# if multiple comparisons is
			if correct_multiple_comparisons:
				print('done')
			# get the correct figure title
			dash_idx = condition.find('-')
			cond1 = condition[:dash_idx]
			cond2 = condition[dash_idx + 1:]
			
			# switch around lumin to always be in second place
			if 'lumin' in cond1:
				cond1 = cond2
				cond2 = 'lumin'
				
			subplot_title = condition_replacements[cond1] + ' vs. ' + condition_replacements[cond2]
			
			x_loc = np.arange(len(subjs_mean_data_lh))
			lh = ax[subplot_row, subplot_col].bar(x_loc - width / 2, subjs_mean_data_lh, width, yerr=subjs_ci_data_lh, color='SkyBlue', label='lh')
			rh = ax[subplot_row, subplot_col].bar(x_loc + width / 2, subjs_mean_data_rh, width, yerr=subjs_ci_data_rh, color='IndianRed', label='rh')
			ax[subplot_row, subplot_col].set_ylim([0.45, .9])
			ax[subplot_row, subplot_col].set_title(subplot_title, fontsize=24)
			ax[subplot_row, subplot_col].set_ylabel('Classification accuracy', fontsize=20)
			ax[subplot_row, subplot_col].set_xticks(x_loc)
			ax[subplot_row, subplot_col].set_xticklabels(roi_xticklabels)
			ax[subplot_row, subplot_col].axhline(y=0.5, linewidth=1.5, color='IndianRed', linestyle='--')
			ax[subplot_row, subplot_col].axhline(y=0.6, linewidth=.5, color='grey', linestyle='--')
			ax[subplot_row, subplot_col].axhline(y=0.7, linewidth=.5, color='grey', linestyle='--')
			ax[subplot_row, subplot_col].axhline(y=0.8, linewidth=.5, color='grey', linestyle='--')
			ax[subplot_row, subplot_col].axhline(y=0.9, linewidth=.5, color='grey', linestyle='--')
			
			autolabel(lh, subjs_ci_data_lh, subjs_p_data_lh, ax[subplot_row, subplot_col], "left")
			autolabel(rh, subjs_ci_data_rh, subjs_p_data_rh, ax[subplot_row, subplot_col], "right")
			
		fig.legend(handles=[lh, rh], labels=['left hemisphere', 'right hemisphere'], loc="lower center", ncol=2, borderaxespad=0.1, fontsize=24)
		#fig.suptitle('Group level CV accuracy', fontsize=24)
		
		plt.savefig(opj(nilearn_dir, 'group_results_space-' + space, 'group_level_' + file_suffix + '_cv_scores.png'), dpi=300)
		tqdm.write('\tImage file: %s' % opj('group_level_' + file_suffix + '_cv_scores.png'))
	