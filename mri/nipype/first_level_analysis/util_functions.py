def createSubjectlist(subjects):

	subjList = []
	
	for subj in subjects:
		subjList.append('sub-' + str(subj).zfill(2))
	
	return subjList


def extractMotionParameters(experimentPath, subj_list):

	import pandas as pd
	import os
	from os.path import join as opj

	for subj in subj_list:
		
		for ses in range(1, 3):
			dataFolder = opj(experimentPath, 'fmriprep', subj, 'ses-0' + str(ses), 'func')

			confoundFiles = [x for x in os.listdir(dataFolder) if '.tsv' in x]

			for conFile in confoundFiles:
				motionFileName = conFile[:-4] + '_mc_params.par'
				motionParameters = pd.read_csv(opj(dataFolder, conFile), delimiter='\t')
				motionParameters.to_csv(path_or_buf=opj(dataFolder, motionFileName), sep='\t', columns=['X','Y','Z', 'RotX', 'RotY', 'RotZ'], header=False, index=False)
		 
def generateOutlierFiles(experimentPath, subj_list, dummyTRs):
	
	import os
	import pandas as pd
	from os.path import join as opj
	
	for subj in subj_list:
		
		for ses in range(1, 3):
			output_dir = opj(experimentPath, 'fmriprep', 'sub-' + str(subj).zfill(2) , 'ses-' + str(ses).zfill(2), 'func')
			
			for runID in range(1, 5):
					
				trialsToExclude = range(dummyTRs)        # exclude first TRs
				
				# create a file for each run in the according folder
				filename = 'sub-' + str(subj).zfill(2) + '_ses-' + str(ses).zfill(2) + '_task-class_run-' + str(runID).zfill(2) + '_bold_outliers.txt'
				df = pd.DataFrame(trialsToExclude)
				df.to_csv(path_or_buf=opj(output_dir, filename), sep='\t', header=False, index=False)
				

def getVolNr(behav_dir, subj_id):
	import pandas as pd
	from os.path import join as opj
	
	logfile = pd.read_csv(opj(behav_dir, 'tp2_stana_Nscans.csv'), delimiter=';')    # read in file
	nrVol = logfile['volumes'].iloc[logfile[logfile['VP'] == int(subj_id[4:])].index.tolist()].tolist()                  # get number of volumes for current subj
	return nrVol
