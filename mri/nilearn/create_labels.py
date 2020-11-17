"""
Script that creates a .tsv file with labels used for decoding. Also creates class_chunk, which pairs 2 runs into one for
cross-valdation
"""

from os.path import join as opj
from os.path import abspath
from sys import platform
import pandas as pd
import os
import glob


# ------------ options ------------
subjects = range(1, 26)
ses_list = range(1, 3)      # sessions that should be included
run_list = range(1, 5)      # runs that should be included
output_space = 'MNI152NLin2009cAsym'        # MNI152NLin2009cAsym

condition_names = ['time', 'dist', 'lumin', 'dots', 'cross']       # which conditions should be included in the

# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
	nilearn_dir = opj(experiment_dir, 'nilearn')
else:
	experiment_dir = abspath('/media/sf_data/fMRI/timePath/')
	nilearn_dir = opj(experiment_dir, 'nilearn')
	
for subj in subjects:
	
	os.system('mkdir -p %s' % opj(experiment_dir, 'nilearn', 'sub-' + str(subj).zfill(2), 'space-' + output_space, 'betas'))
		
	trial_types, trial_duration, trial_distance, trial_num, trial_speed = [], [], [], [], []         # create list for durations
	runs = []                   # create list for runs
	sessions = []               # create corresponding list for session numbers
	class_chunk = []            # create list for classification chunks (combining 2 runs into one)
	subject_list = []
	sub_str = 'sub-' + str(subj).zfill(2)   # for readability
	
	for ses in ses_list:
		
		# read in the log file for trial_distance information
		log_file = pd.read_csv(glob.glob(opj(experiment_dir, 'behavioural', 'logfiles', 'tp2_' + str(subj).zfill(2)
		                                     + '_' + str(ses).zfill(2) + '_*.txt'))[0],
		                       sep=' ', header=None)

		for run in run_list:
			# read in event file
			ses_str = 'ses-' + str(ses).zfill(2)    # for readability
			run_str = 'run-' + str(run).zfill(2)    # for readability
			        
			event_file = pd.read_csv(opj(experiment_dir, 'bids', sub_str, ses_str, 'func',
			                             sub_str + '_' + ses_str + '_task-class' + '_' + run_str + '_events.tsv'),
			                         delimiter='\t')
			
			events = event_file.trial_type
			durations = event_file.duration
			
			# create another column "class_chunks" to combine two runs into one for even trial distribution
			if ses == 1:
				if run == 1 or run == 2:
					chunk = 1
				elif run == 3 or run == 4:
					chunk = 2
			else:
				if run == 1 or run == 2:
					chunk = 3
				elif run == 3 or run == 4:
					chunk = 4

			# store temporarily in a list (could be added to the data frame later on)
			for event in range(len(events)):
				if events.iloc[event] in condition_names:       # first check if the trial should be included
					
					trial = event_file['trial_nr'].iloc[event]
					
					# remove for every odd run the last trial and for every even run the last two trials, also remove the first two trials
					if (run % 2 == 0 and trial == 65) or (trial == 1 or trial == 2):
						#print('Removed trial %s of run %s' %(trial, run))
						pass
					else:
						if ses == 1:
							runs.append(run)
							if run == 1 or run == 2:
								class_chunk.append(1)
							elif run == 3 or run == 4:
								class_chunk.append(2)
						else:
							runs.append(run + 4)
							if run == 1 or run == 2:
								class_chunk.append(3)
							elif run == 3 or run == 4:
								class_chunk.append(4)
					
						sessions.append(ses)
						trial_types.append(events.iloc[event])
						
						# assess duration and create trial duration
						if durations.iloc[event] <= 4:
							trial_duration.append('low')
						else:
							trial_duration.append('high')

						# append subject
						subject_list.append(subj)

						# access distance and create trial distance
						for line in log_file.values:
							event_type = line[2]
							ses_nr = line[4]
							run_nr = line[5]
							trial_nr = line[6]
							speed = float(line[9])/float(line[8])
							
							if event_type == 'movementOnset' and int(ses_nr) == ses and int(run_nr) == run and int(trial_nr) == trial:
								if float(line[9]) > 19:
									trial_distance.append('high')
								else:
									trial_distance.append('low')
								if float(line[10]) > 70:
									trial_num.append('high')
								else:
									trial_num.append('low')
								if speed > 7:
									trial_speed.append('high')
								elif 3 < speed and speed < 6:
									trial_speed.append('mid')
								elif speed < 4:
									trial_speed.append('low')

							elif event_type == 'crossOnset' and int(ses_nr) == ses and int(run_nr) == run and int(trial_nr) == trial:
								trial_distance.append('na')
								trial_num.append('na')
								trial_speed.append('na')

	# create a dataframe for output
	trialData = pd.DataFrame({'subject': subject_list,
	                          'trial': trial_types,
	                          'duration': trial_duration,
	                          'distance': trial_distance,
	                          'nrdots': trial_num,
	                          'speed': trial_speed,
	                          'run': runs,
	                          'session': sessions,
	                          'class_chunk': class_chunk})
		
	# save subject specific label files in nilearn data folder
	trialData.to_csv(opj(nilearn_dir, sub_str, 'space-' + output_space, 'betas', sub_str + '_merged_events.tsv'),
	                 sep="\t", columns=trialData.columns, index=False)
	print('Written file %s' % (opj(experiment_dir, 'nilearn', sub_str, 'space-' + output_space, 'betas', sub_str + '_merged_events.tsv')))