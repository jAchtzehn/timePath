from os.path import join as opj
from os.path import abspath
from sys import platform
import os
from scipy.spatial.distance import squareform
from mvpa2.suite import fmri_dataset, zscore, mean_group_sample, remove_nonfinite_features
from mvpa2.measures import rsa
from mvpa2.measures.searchlight import sphere_searchlight
from mvpa2.datasets.mri import map2nifti
from mvpa2.base.learner import ChainLearner
from mvpa2.mappers.shape import TransposeMapper
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from nilearn.plotting import plot_stat_map
from nilearn.image.resampling import coord_transform
from nilearn.image import smooth_img, load_img, new_img_like
from utilFunctions import create_html_report
import cPickle
import itertools
from sklearn.preprocessing import MinMaxScaler

# ------------ File I/O ------------
if platform == 'darwin':
    experiment_dir = abspath('/Volumes/Seagate/Backups/16102020/data/timePath')
else:
    experiment_dir = abspath('/home/achtzehnj/data/timePath/')

nilearn_dir = opj(experiment_dir, 'derivatives', 'nilearn')
rsa_dir = opj(experiment_dir, 'derivatives', 'rsa')
behav_dir = opj(experiment_dir, 'derivatives', 'behavioural')

os.system('mkdir -p %s' % rsa_dir)                 # create the mvpa2 folder for data storage

# ------------ options ------------
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]     # subjects to plot
# [4, 6, 7, 9, 12, 14, 17, 19, 23, 24]                                                                  # 10 best subjects (TSNR >= 50)
# [1, 2, 3, 5, 8, 10, 11, 13, 15, 16, 18, 20, 21, 22, 25]                                               # the rest
# [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]               # all except 15
conditions = ['time', 'dist', 'dots']
file_suffix = 'cdist'
sl_radius = 3
spaces = ['MNI152NLin2009cAsym']

os.system('clear')                              # clean up terminal
# if file_suffix is empty, add _
if file_suffix == '':
    file_suffix = '_'
else:
    file_suffix = '_' + file_suffix + '_'

minmaxscaler = MinMaxScaler(feature_range=(1, 2))
for space in spaces:
    print('Calculating space: %s' % space)

    for subj in subjects:
        subj_str = 'sub-' + str(subj).zfill(2)      # create the subject string (just for easier handling)

        print('Subject: %s' % subj_str)
        result_dir = opj(rsa_dir, subj_str, 'space-' + space, 'results')
        os.system('mkdir -p %s' % result_dir)                 # create the mvpa2 folder for data storage

        # -- data loading
        data_path = opj(nilearn_dir, subj_str, 'space-' + space, 'betas', subj_str + '_betas_merged.nii.gz')
        behav_path = opj(nilearn_dir, subj_str, 'space-' + space, 'betas', subj_str + '_merged_events.tsv')

        mask_path = opj(nilearn_dir, 'group_masks', 'space-' + space, 'group_mask_wb_binarized.nii.gz')  # mask filename
        # -- calculate participants' crossDim matrix
        pse_diff_data = pd.read_csv(opj(behav_dir, 'pse_data_cross_dim_individual_norm.tsv'), delimiter='\t')
        crossDim = {}
        for condition_pair in itertools.permutations(conditions, 2):
            crossDim['{}_vs_{}'.format(condition_pair[0], condition_pair[1])] = np.mean(np.abs(pse_diff_data['pse_diff'][(((pse_diff_data.condition_rel == condition_pair[0]) &
                                                             (pse_diff_data.condition_irrel == condition_pair[1])) |
                                                              ((pse_diff_data.condition_rel == condition_pair[1]) &
                                                               (pse_diff_data.condition_irrel == condition_pair[0]))) &
                                                             (pse_diff_data.subject == subj)].values))

        crossDim_matrix = -1 * np.array([crossDim['dist_vs_dots'], crossDim['time_vs_dist'], crossDim['time_vs_dots']])
        minmaxscaler.fit(crossDim_matrix.reshape(-1, 1))
        crossDim_matrix_norm = minmaxscaler.transform(crossDim_matrix.reshape(-1, 1))

        # read in behav data
        behav_data = pd.read_csv(behav_path, delimiter='\t')
        targets = behav_data['trial'].tolist()                  # labels are called targets in mvpa
        chunks = behav_data['class_chunk'].tolist()             # chunks are runpairs of 2 runs
        runs = behav_data['run'].tolist()

        # create dataset
        fmri_data = fmri_dataset(data_path, mask=mask_path, targets=targets, chunks=chunks)

        fmri_data.sa['runs'] = runs
        fmri_data = remove_nonfinite_features(fmri_data)

        # -- preprocessing
        zscore(fmri_data, chunks_attr='chunks')                     # zscoring the data to remove mean
        fmri_data = fmri_data[fmri_data.sa.targets != 'cross']      # remove cross from the conditions (targets)

        # remove cross condition
        for cond in ['time', 'dist', 'dots', 'lumin']:
            if cond not in conditions:
                fmri_data = fmri_data[fmri_data.sa.targets != cond]      # remove cross from the conditions (targets)

        # -- computation
        # first compute RDMs
        mtgs = mean_group_sample(['targets'])       # create the mean of all trials
        mtds = mtgs(fmri_data)

        dsm = rsa.CDist(auto_train=True, force_train=True)
        #dsm = rsa.PDist(square=False, center_data=True)

        sl = sphere_searchlight(dsm, radius=sl_radius, nproc=46)
        slres = sl(mtds)

        rdm_data_new_order = np.zeros(shape=(3, len(slres.samples[0, :])))
        rdm_data_new_order[0, :] = slres.samples[2, :]  # time vs. dist
        rdm_data_new_order[1, :] = slres.samples[5, :]  # time vs. dots
        rdm_data_new_order[2, :] = slres.samples[1, :]  # dots vs dist

        rdm_data_new_order = np.nan_to_num(rdm_data_new_order)

        # now compute RDM consistency by cross-validation with chunks
        mtcgs = mean_group_sample(['targets', 'chunks'])
        mtcds = mtcgs(fmri_data)
        dscm = rsa.PDistConsistency(center_data=True)

        sl_stab = sphere_searchlight(dscm, radius=sl_radius, nproc=46)
        slres_stab = sl_stab(mtcds)
        mean_consistency = np.mean(slres_stab, axis=0)

        # now compute similarity with participants' behavioral matrix
        target_matrix = np.array([crossDim_matrix_norm[0], crossDim_matrix_norm[1], crossDim_matrix_norm[2]]).reshape(-1, )
        t = slres.samples[:, mean_consistency.argmax()]
        tdsm = rsa.PDistTargetSimilarity(target_matrix)
        sl_tdsm = sphere_searchlight(ChainLearner([tdsm, TransposeMapper()]), radius=sl_radius, nproc=46)
        slres_tdsm = sl_tdsm(mtds)

        #target_data_idx = slres_tdsm.samples[1, :] > 0.05
        #slres_tdsm.samples[0, target_data_idx] = 0

        niimg_targ = map2nifti(fmri_data, np.nan_to_num(slres_tdsm.samples[0, :]))
        niimg_targ.to_filename(opj(result_dir,
                              subj_str + '_space-' + space + file_suffix + 'wb_targ_sim_values.nii.gz'))

        niimg_res = map2nifti(fmri_data, rdm_data_new_order)
        niimg_res.to_filename(opj(result_dir,
                              subj_str + '_space-' + space + file_suffix + 'wb_rdm_values.nii.gz'))

        niimg_cons = map2nifti(fmri_data, mean_consistency)
        niimg_cons.to_filename(opj(result_dir,
                              subj_str + '_space-' + space + file_suffix + 'wb_cons_values.nii.gz'))

        niiimg_cons_tmp = load_img(opj(result_dir,
                              subj_str + '_space-' + space + file_suffix + 'wb_cons_values.nii.gz'))
        niiimg_cons_smooth = smooth_img(niiimg_cons_tmp, 6)
        niiimg_cons_smooth.to_filename(opj(result_dir,
                              subj_str + '_space-' + space + file_suffix + 'wb_cons_smoothed_values.nii.gz'))