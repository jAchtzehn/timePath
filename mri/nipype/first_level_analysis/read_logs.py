def get_classification_onsets(subj_id, base_dir, addRegressors, conditions):

	from nipype.interfaces.base import Bunch
	from os.path import join as opj
	import pandas as pd
	
	subject_info = []           # final list containing the subject info per run
	
	# go through each session
	for ses in range(1, 3):
		# go through each run and create subject info bunch
		for runID in range(4):
			
			# read in events.tsv file
			df = pd.read_csv(opj(base_dir, 'bids', subj_id, 'ses-0' + str(ses), 'func', subj_id + '_ses-0' + str(ses) + '_task-class_run-0' + str(runID + 1) + '_events.tsv'),
		                 delimiter='\t')
			
			onsets = []         # list containing sublists for each run
			durations = []      # list containing sublists for each run
			
			# for each condition, read out onsets and duration
			for condition in conditions:
				
				idx = df[df['trial_type'] == condition].index
				onsets_cond = df['onset'].iloc[idx].values.tolist()
				durations_cond = df['duration'].iloc[idx].values.tolist()
				
				onsets.append(onsets_cond)
				durations.append(durations_cond)
				
			# get corrected runID for subject_info --> b/c in subject_info the runIDs will run from 0 - 7, in the files they are from 1 - 4 with 2 sessions
			if ses == 2:
				runID_subjInfo = runID + 4
			else:
				runID_subjInfo = runID
				
			if not addRegressors:
				subject_info.insert(runID_subjInfo,
				                    Bunch(conditions=conditions,
				                          onsets=onsets,
				                          durations= durations))
			else:
				# read additional regressors for run
				confoundFile = opj(base_dir, 'fmriprep', subj_id, 'ses-0' + str(ses), 'func', subj_id + '_ses-0' + str(ses) +  '_task-class_run-0' + str(runID + 1) + '_bold_confounds_cut.tsv')
				
				regressors = pd.read_csv(confoundFile, delimiter='\t', header=0)
				
				subject_info.insert(runID_subjInfo,
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
				
	return subject_info
	
	
	
	
	