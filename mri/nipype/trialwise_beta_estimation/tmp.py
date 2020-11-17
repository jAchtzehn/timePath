from os.path import join as opj
from os.path import abspath
import pandas as pd
import os

behav_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/behavioural/onset_and_duration_files_with_warning')  # where is the onset and duration information stored?
target_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/behavioural/onset_and_duration_files')
for sub in range(1, 26):

    sub_str = 'sub-' + str(sub).zfill(2)

    for ses in range(1, 3):
        for run in range(1, 5):
            warn_ons_fname = opj(behav_dir, sub_str, 'ses-' + str(ses).zfill(2),
                                 'func', 'run-' + str(run).zfill(2), 'onsdur_stana', 'warning')

            f_ons = open(warn_ons_fname + '.ons')
            filedata_ons = f_ons.read().splitlines()
            filedata_ons = [float(x) for x in filedata_ons if x != '']
            f_ons.close()

            print('subj: {}, ses: {}, run: {}, onsets: {}'.format(str(sub), str(ses), str(run), str(filedata_ons)))

            # os.system('cp {} {}'.format(warn_ons_fname + '.ons',
            #                             opj(target_dir, sub_str, 'ses-' + str(ses).zfill(2),
            #                                 'run-' + str(run).zfill(2),
            #                                 'onsdur_stana')))
            #
            # os.system('cp {} {}'.format(warn_ons_fname + '.dur',
            #                             opj(target_dir, sub_str, 'ses-' + str(ses).zfill(2),
            #                                 'run-' + str(run).zfill(2),
            #                                 'onsdur_stana')))