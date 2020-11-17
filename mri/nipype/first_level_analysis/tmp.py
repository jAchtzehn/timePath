from os.path import join as opj
from os.path import abspath


behav_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/behavioural')  # where is the onset and duration information stored?


import pandas as pd
from os.path import join as opj

subj_id = 'sub-04'

logfile = pd.read_csv(opj(behav_dir, 'tp2_stana_Nscans.csv'), delimiter=';')    # read in file
nrVol = logfile['volumes'].iloc[logfile[logfile['VP'] == int(subj_id[4:])].index.tolist()].tolist()                  # find all indices of the current

print('done')