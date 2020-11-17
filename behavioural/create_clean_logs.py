"""
Script that creates a .tsv file with labels used for decoding. Also creates class_chunk, which pairs 2 runs into one for
cross-valdation
"""

from os.path import join as opj
from os.path import abspath
import pandas as pd
import os
import glob

# ------------ options ------------
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]

ses_list = range(1, 3)  # sessions that should be included
run_list = range(1, 5)  # runs that should be included
trial_list = range(1, 67)
exclude_events = ['TriggerMRI', 'startOfBlock', 'showVAS', 'getVASresponse', 'pauseInstrOnset']
write_data = True

# ------------ File I/O ------------
experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
behavioral_dir = opj(experiment_dir, 'behavioural')

data_all_subj_rt = pd.DataFrame()
data_all_subj_vas = pd.DataFrame()

total_nr_trials = 0
too_slow_rt_nr_trials = 0
made_resp_forgot_nr_trials = 0
rt_below_thresh = 0

for subj in subjects:

    sub_str = 'sub-' + str(subj).zfill(2)  # for readability

    os.system(
        'mkdir -p %s' % opj(behavioral_dir, sub_str))

    trial_nr, trial_type, trial_type_nr, std_time, std_dist, std_num, std_lumin = [], [], [], [], [], [], []  # create list for values
    trial_time, trial_dist, trial_num = [], [], []
    compVal, RT, resp = [], [], []
    trial_type_vas, VASresponse, run_vas, sessions_vas = [], [], [], []

    subject_nr = []
    runs = []  # create list for runs
    sessions = []  # create corresponding list for session numbers

    for ses in ses_list:

        # read in the log file for trial_distance information
        log_file = pd.read_csv(glob.glob(opj(experiment_dir, 'behavioural', 'logfiles', 'tp2_' + str(subj).zfill(2)
                                             + '_' + str(ses).zfill(2) + '_*.txt'))[0],
                               sep=' ', header=None)

        for run in run_list:

            # get each trial info
            for trial in trial_list:
                total_nr_trials += 1
                # go through each line and find start and end indices of current trial
                idc_current_trial = []

                for idx, line in enumerate(log_file.values):

                    if line[2] not in exclude_events:
                        if int(line[4]) == ses and int(line[5]) == run and int(line[6]) == trial:
                            idc_current_trial.append(idx)

                log_file_current_trial = log_file[idc_current_trial[0]:idc_current_trial[len(idc_current_trial) - 1] + 1]   # data for the current trial

                if 'instrRespondFasterOffset' in log_file_current_trial.loc[:, 2].values:
                    too_slow_rt_nr_trials += 1
                if 'discrimResponse' in log_file_current_trial.loc[:,2].values and 'thumb_button' in log_file_current_trial.loc[:, 2].values:
                    made_resp_forgot_nr_trials += 1

                # save only when response was made
                if 'discrimResponse' in log_file_current_trial.loc[:, 2].values and 'thumb_button' not in log_file_current_trial.loc[:, 2].values:
                    trial_data = log_file_current_trial[log_file_current_trial.loc[:, 2].str.match('discrimResponse')].values.transpose()

                    trial_type_nr.append(int(trial_data[7]))

                    if int(trial_data[7]) == 0:
                        trial_type.append('time')
                    elif int(trial_data[7]) == 1:
                        trial_type.append('dist')
                    elif int(trial_data[7]) == 2:
                        trial_type.append('dots')
                    elif int(trial_data[7]) == 3:
                        trial_type.append('lumin')

                    if ses == 1:
                        runs.append(run)
                    else:
                        runs.append(run + 4)

                    subject_nr.append(subj)
                    trial_nr.append(trial)

                    sessions.append(ses)
                    std_time.append(float(trial_data[8]))
                    std_dist.append(float(trial_data[9]))
                    std_num.append(float(trial_data[10]))
                    std_lumin.append(float(trial_data[11]))
                    trial_time.append(round(float(trial_data[12]), 3))
                    trial_dist.append(round(float(trial_data[13]), 3))
                    trial_num.append(round(float(trial_data[14]), 3))
                    compVal.append(float(trial_data[15]))
                    RT.append(round(float(trial_data[16]), 5))
                    # response
                    if int(trial_data[17]) == 1:
                        resp.append(0)
                    else:
                        resp.append(1)

                    if float(trial_data[16]) < 0.3:
                        rt_below_thresh += 1

            # get VAS difficulty rating
            idc_current_run = []
            for idx, line in enumerate(log_file.values):
                if int(line[4]) == int(ses) and int(line[5]) == int(run):
                    idc_current_run.append(idx)

            log_file_current_run = log_file[idc_current_run[0]:idc_current_run[len(idc_current_run) - 1]]  # data for the current trial

            for line in log_file_current_run.values:
                if line[2] == 'getVASresponse':
                    VASresponse.append(round(float(line[17]), 2))

                    if ses == 1:
                        run_vas.append(run)
                    else:
                        run_vas.append(run + 4)

                    sessions_vas.append(ses)
                    if int(line[7]) == 0:
                        trial_type_vas.append('time')
                    elif int(line[7]) == 1:
                        trial_type_vas.append('dist')
                    elif int(line[7]) == 2:
                        trial_type_vas.append('dots')
                    elif int(line[7]) == 3:
                        trial_type_vas.append('lumin')

    # create a dataframe for trial output
    trialData = pd.DataFrame({'subject': subject_nr, 'trial_nr': trial_nr, 'trial_type': trial_type, 'trial_type_nr': trial_type_nr, 'run': runs, 'session': sessions,
                              'std_time': std_time, 'std_dist': std_dist, 'std_dots': std_num, 'std_lumin': std_lumin,
                              'trial_time': trial_time, 'trial_dist': trial_dist, 'trial_num': trial_num,
                              'compVal': compVal, 'RT': RT, 'resp': resp
                              })

    # create a dataframe for VAS output
    vasData = pd.DataFrame({'subject': [subj] * len(VASresponse), 'run': run_vas, 'session': sessions_vas,
                            'trial_type': trial_type_vas, 'vas': VASresponse})

    # summarize into one dataFrame
    data_all_subj_rt = data_all_subj_rt.append(trialData, ignore_index=True)
    data_all_subj_vas = data_all_subj_vas.append(vasData, ignore_index=True)

    if write_data:
        # save subject specific label files in nilearn data folder
        trialData.to_csv(opj(behavioral_dir, sub_str, sub_str + '_merged_events.tsv'),
                         sep="\t", columns=trialData.columns,
                         index=False)

        # write VAS data
        vasData.to_csv(opj(behavioral_dir, sub_str, sub_str + '_vas_responses.tsv'),
                       sep="\t", columns=vasData.columns, index=False)

        print('Written file %s' % (
            opj(behavioral_dir, sub_str, sub_str + '_merged_events.tsv')))

        print('Written file %s' % (
            opj(behavioral_dir, sub_str, sub_str + '_vas_responses.tsv')))

if write_data:
    # save summary file
    data_all_subj_rt.to_csv(opj(behavioral_dir, 'RT_data_all_subj.tsv'), sep='\t', columns=data_all_subj_rt.columns,
                            index=False)
    data_all_subj_vas.to_csv(opj(behavioral_dir, 'VAS_data_all_subj.tsv'), sep='\t', columns=data_all_subj_vas.columns,
                             index=False)
    print('Written file %s' % (
        opj(behavioral_dir, 'RT_data_all_subj.tsv')))

    print('Written file %s' % (
        opj(behavioral_dir, 'VAS_data_all_subj.tsv')))

# print('Total nr of trials: %s, too slow: %s(%s), forgot cue: %s(%s), below thresh:%s(%s)' %
#       (total_nr_trials, too_slow_rt_nr_trials, too_slow_rt_nr_trials/total_nr_trials, made_resp_forgot_nr_trials, made_resp_forgot_nr_trials/total_nr_trials,
#        rt_below_thresh, rt_below_thresh/total_nr_trials))
