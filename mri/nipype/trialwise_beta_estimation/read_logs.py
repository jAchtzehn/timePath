"""
Gets the onsets and durations of the desired events (defined by input "conditions") and creates a
bunch variable

(c) Johannes Achtzehn, 06.05.2019
"""

def get_classification_onsets(subj_id, base_dir, ses_id, run_id, trial_id, addRegressors, conditions):
	
	from nipype.interfaces.base import Bunch
	from os.path import join as opj
	import pandas as pd
	
	subject_info = []           # final list containing the subject info per run
	onsets = []                 # list containing sublists for each condition
	durations = []              # list containing sublists for each condition
	
	# read in events.tsv file
	df_all = pd.read_csv(opj(base_dir, 'bids', subj_id, 'ses-0' + str(ses_id), 'func', subj_id + '_ses-0' + str(ses_id) + '_task-class_run-0' + str(run_id) + '_events.tsv'),
                 delimiter='\t')
	
	# for the current trial, delete all the other events in the event file except time, dist, lumin and dots to get the correct onset
	df_trial = df_all[(df_all.trial_nr == trial_id) & ((df_all.trial_type != 'icon') & (df_all.trial_type != 'resp') &
	                                                   (df_all.trial_type != 'comp') & (df_all.trial_type != 'warning'))]
	trial_idx = df_all[(df_all.trial_nr == trial_id) & ((df_all.trial_type != 'icon') & (df_all.trial_type != 'resp') &
	                                                    (df_all.trial_type != 'comp') & (df_all.trial_type != 'warning'))].index.values
	
	# add current trial in a single list for this trial
	onsets.append(df_trial['onset'].tolist())
	durations.append(df_trial['duration'].tolist())
	
	# insert information into trial_info, this is later used to write out current onset and dur in a csv file
	trial_info = pd.DataFrame({'trial_type': df_all.trial_type[trial_idx].values,
	                           'onset': df_trial['onset'].values,
	                           'duration': df_trial['duration'].values})
	
	# now go through the other conditions and add these as additional regressors
	conditions.remove('cross')              # first remove 'cross', as we do not want to model zero trials
	
	regr_info = pd.DataFrame({'trial_type': [],
	                          'onset': [],
	                          'duration': []})
	
	for condition in conditions[1:]:
			idx = df_all[df_all['trial_type'] == condition].index.tolist()
			
			# remove the current trial from the indices
			if trial_idx in idx:
				idx.remove(trial_idx)
			
			trialTypes_cond = df_all['trial_type'].iloc[idx].values.tolist()
			onsets_cond = df_all['onset'].iloc[idx].values.tolist()
			durations_cond = df_all['duration'].iloc[idx].values.tolist()
			
			# check if the current run includes a thumb press, if not remove thumb from the conditions
			if not onsets_cond == []:
				onsets.append(onsets_cond)
				durations.append(durations_cond)
				# add info to regr info
				regr_info = regr_info.append(pd.DataFrame({'trial_type': trialTypes_cond, 'onset': onsets_cond, 'duration':durations_cond}),
				                             ignore_index=True)
				
			else:
				conditions.remove(condition)
	
	regr_info = regr_info.sort_values(by='onset')   # sort by onset
	regr_info = regr_info.reset_index(drop=True)    # update index
	
	if not addRegressors:
		subject_info.insert(run_id,
		                    Bunch(conditions=conditions,
		                          onsets=onsets,
		                          durations= durations))
	else:
		# read additional regressors for run
		confoundFile = opj(base_dir, 'fmriprep', subj_id, 'ses-0' + str(ses_id), 'func', subj_id + '_ses-0' + str(ses_id) +  '_task-class_run-0' + str(run_id) + '_bold_confounds_cut.tsv')
		
		regressors = pd.read_csv(confoundFile, delimiter='\t', header=0)
		
		subject_info.insert(run_id,
		                    Bunch(conditions=conditions,
		                          onsets=onsets,
		                          durations=durations,
		                          regressors=[list(regressors.FramewiseDisplacement.fillna(0)),
		                                      list(regressors.aCompCor00),
		                                      list(regressors.aCompCor01),
		                                      list(regressors.aCompCor02),
		                                      list(regressors.aCompCor03),
		                                      list(regressors.aCompCor04),
		                                      list(regressors.aCompCor05),
		                                      ],
		                          regressor_names = ['FramewiseDisplacement', 'aCompCor00',
		                                       'aCompCor01', 'aCompCor02', 'aCompCor03',
		                                       'aCompCor04', 'aCompCor05']))
	
	return subject_info, trial_info, regr_info
	
	
	
	
	