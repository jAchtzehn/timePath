from os.path import join as opj
from os.path import abspath
from sys import platform
import os
from mvpa2.suite import fmri_dataset, zscore, mean_group_sample, remove_nonfinite_features
from mvpa2.measures import rsa
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
from utilFunctions import create_html_report

# ------------ File I/O ------------
if platform == 'darwin':
    experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
    experiment_dir = abspath('/media/sf_data/fMRI/timePath/')

nilearn_dir = opj(experiment_dir, 'nilearn')
rsa_dir = opj(experiment_dir, 'rsa')

os.system('mkdir -p %s' % rsa_dir)                 # create the mvpa2 folder for data storage

# ------------ options ------------
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]     # subjects to plot
# [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]  # all except 15
# [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]  # all
conditions = ['time', 'dist', 'dots']
masks = ['ips_lh', 'ips_rh', 'ifg_lh', 'ifg_rh', 'hc_lh', 'hc_rh', 'V5_rh', 'V5_lh', 'insula_rh', 'insula_lh']
# ['ips_lh', 'ips_rh', 'ifg_lh', 'ifg_rh', 'hc_lh', 'hc_rh', 'V5_rh', 'V5_lh', 'insula_rh', 'insula_lh']
plot_html = False
plot_group_level = True
spaces = ['MNI152NLin2009cAsym']
file_suffix = 'centered'
os.system('clear')                              # clean up terminal

for space in spaces:
    
    group_result_dir = opj(rsa_dir, 'group_results_space-' + space)
    os.system('mkdir -p %s' % group_result_dir)
    
    rdm_images_list = []
    # init dict for saving the RDMs for all subjects
    rdms_cummulative_allSubj = {}
    rdms_stab_allSubj = {}
    rdms_allSubj = {}

    for roi in masks:
        rdms_cummulative_allSubj[roi] = np.zeros((len(conditions), len(conditions)))
        rdms_stab_allSubj[roi] = []
        rdms_allSubj[roi] = []

    for subj in subjects:
        subj_str = 'sub-' + str(subj).zfill(2)      # create the subject string (just for easier handling)
        
        print('Subject: %s' % subj_str)
        result_dir = opj(rsa_dir, subj_str, 'space-' + space, 'results')
        os.system('mkdir -p %s' % result_dir)                 # create the mvpa2 folder for data storage
       
        fig, ax = plt.subplots(2, len(masks)/2, figsize=(18, 8))
        plt.tight_layout(5)
        fig.suptitle('Subject %s' % str(subj), fontsize=20)
        fig.subplots_adjust(hspace=0.35)
        
        for roi_enum, roi in enumerate(masks):
            
            if roi_enum < len(masks)/2:
                subplot_row = 0
                subplot_col = roi_enum
            else:
                subplot_row = 1
                subplot_col = roi_enum - len(masks)/2
     
            # -- data loading
            data_path = opj(nilearn_dir, subj_str, 'space-' + space, 'betas', subj_str + '_betas_merged.nii.gz')
            if roi == 'V5_rh' or roi == 'V5_lh':
                mask_path = opj(nilearn_dir, 'group_masks', 'space-' + space,
                                'group_mask_' + roi + '_binarized.nii.gz')  # mask filename
            else:
                mask_path = opj(nilearn_dir, subj_str, 'space-' + space, 'masks',
                                subj_str + '_' + roi + '_mask_binarized.nii.gz')  # mask filename

            behav_path = opj(nilearn_dir, subj_str, 'space-' + space, 'betas', subj_str + '_merged_events.tsv')
            
            # read in behav data
            behav_data = pd.read_csv(behav_path, delimiter='\t')
            targets = behav_data['trial'].tolist()                  # labels are called targets in mvpa
            chunks = behav_data['class_chunk'].tolist()             # chunks are runpairs of 2 runs
            runs = behav_data['run'].tolist()
            
            # create dataset
            fmri_data = fmri_dataset(data_path, mask=mask_path, targets=targets, chunks=chunks)
            fmri_data.sa['runs'] = runs
            fmri_data = remove_nonfinite_features(fmri_data)

            # # -- preprocessing
            zscore(fmri_data)                                           # zscoring the data to remove mean
            fmri_data = fmri_data[fmri_data.sa.targets != 'cross']      # remove cross from the conditions (targets)

            for cond in ['time', 'dist', 'dots', 'lumin']:
                if cond not in conditions:
                    fmri_data = fmri_data[fmri_data.sa.targets != cond]      # remove cross from the conditions (targets)

            # -- computation
            # first compute RDMs for all trials
            mtgs = mean_group_sample(['targets'])       # create the mean of all trials
            mtds = mtgs(fmri_data)
            dsm = rsa.PDist(pairwise_metric='correlation', center_data=True)
            res = dsm(mtds)

            # # second compute CV RDMs for all trials
            # mtgs_cv = mean_group_sample(['targets', 'chunks'])  # create the mean of all trials
            # mtds_cv = mtgs_cv(fmri_data)
            # dsm_cv = rsa.PDist(square=True, pairwise_metric='mahalanobis')
            # res_cv = dsm_cv(mtds_cv)
            #
            # rdm_cv = np.zeros([3, 3])
            # for i in list(range(1, 11, 3)):
            #     m = res_cv.samples[i-1:i+2, i-1:i+2]
            #     rdm_cv += res_cv.samples[i-1:i+2, i-1:i+2]
            # rdm_cv = rdm_cv / 4

            # third compute stability of RDMs
            # pattern stability across chunks analysis
            mtgs_stab = mean_group_sample(['targets', 'chunks'])
            mtds_stab = mtgs_stab(fmri_data)
            dsm_stab = rsa.PDistConsistency(center_data=True)
            res_stab = dsm_stab(mtds_stab)
            
            # plotting
            rdm_plot = ax[subplot_row, subplot_col].imshow(res, interpolation='nearest')
            div = make_axes_locatable(ax[subplot_row, subplot_col])
            cax = div.append_axes("right", size="5%", pad=0.2)
            cbar = plt.colorbar(rdm_plot, cax=cax)
            rdm_plot.set_clim(vmin=0, vmax=2)
            ax[subplot_row, subplot_col].set_xticklabels(list(mtds.sa.targets), fontdict=None, minor=False, rotation=-45)
            ax[subplot_row, subplot_col].set_xticks(range(len(res)))
            ax[subplot_row, subplot_col].set_yticklabels(list(mtds.sa.targets), fontdict=None, minor=False)
            ax[subplot_row, subplot_col].set_yticks(range(len(res)))
            ax[subplot_row, subplot_col].set_title(roi + '; RDM stability: %s'
                                                   % '{:.2f}'.format(np.mean(res_stab.samples)))
            
            # add text of value to each field
            for i in range(len(res)):
                for j in range(len(res)):
                    if i != j:
                        ax[subplot_row, subplot_col].text(j, i, '{:.2f}'.format(res.samples[i, j]), ha='center', va='center', color='w')
        
            # save RDM for group level
            rdms_cummulative_allSubj[roi] = rdms_cummulative_allSubj[roi] + res.samples
            rdms_stab_allSubj[roi].append(np.mean(res_stab.samples))
            rdms_cummulative_allSubj['targets'] = mtds.sa.targets

            rdms_allSubj[roi].append(res.samples)
            
        plt.savefig(opj(result_dir, subj_str + '_space-' + space + '_RDM_roi_' + file_suffix + '.png'), dpi=300)
        rdm_images_list.append(opj(result_dir, subj_str + '_space-' + space + '_RDM_roi_' + file_suffix + '.png'))
        plt.close()
    
    # plot individual RDMs
    if plot_html:
        create_html_report(rdm_images_list, opj(group_result_dir, 'RDM_space-' + space + '_ ' + file_suffix + '_roi.html'), title='RDMs for all subjects')
    
    if plot_group_level:

        # write CSV
        # write .csv file with results
        pd_subjId, pd_roi, pd_cond1, pd_cond2, pd_dm, pd_cons = [], [], [], [], [], []

        for s, subjId in enumerate(subjects):
            for roi in masks:
                for i, cond_comb in enumerate([['dots', 'dist'], ['dist', 'time'], ['dots', 'time']]):
                    pd_subjId.append(subjId)
                    pd_roi.append(roi)
                    pd_cond1.append(cond_comb[0])
                    pd_cond2.append(cond_comb[1])

                    if i == 0:
                        pd_dm.append(rdms_allSubj[roi][s][0, 1])
                    elif i == 1:
                        pd_dm.append(rdms_allSubj[roi][s][0, 2])
                    else:
                        pd_dm.append(rdms_allSubj[roi][s][1, 2])

                    pd_cons.append(rdms_stab_allSubj[roi][s])

        df = pd.DataFrame({'subject': pd_subjId,
                           'roi': pd_roi,
                           'condition_1': pd_cond1,
                           'condition_2': pd_cond2,
                           'dissimilarity': pd_dm,
                           'consistency': pd_cons})

        df.to_csv(opj(group_result_dir, 'RDM_space-MNI152NLin2009cAsym_roi_' + file_suffix + '.csv'), sep='\t',
                  columns=['subject', 'roi', 'condition_1', 'condition_2', 'dissimilarity', 'consistency'],
                  index=False)

        # -- plot group level RDM
        fig, ax = plt.subplots(2, len(masks)/2, figsize=(18, 8))
        plt.tight_layout(5)
        fig.suptitle('Mean RDMs', fontsize=20)
        fig.subplots_adjust(hspace=0.35)
        
        # find max
        max_vals = []
        for roi in masks:
            max_vals.append(np.max(rdms_cummulative_allSubj[roi] / len(subjects)))
            
        for roi_enum, roi in enumerate(masks):
            rdms_cummulative_allSubj[roi] = rdms_cummulative_allSubj[roi] / len(subjects)               # compute mean
            
            if roi_enum < len(masks)/2:
                subplot_row = 0
                subplot_col = roi_enum
            else:
                subplot_row = 1
                subplot_col = roi_enum - len(masks)/2
        
            # plotting
            rdm_plot = ax[subplot_row, subplot_col].imshow(rdms_cummulative_allSubj[roi], interpolation='nearest')
            div = make_axes_locatable(ax[subplot_row, subplot_col])
            cax = div.append_axes("right", size="5%", pad=0.2)
            cbar = plt.colorbar(rdm_plot, cax=cax)
            rdm_plot.set_clim(vmin=0, vmax=np.max(max_vals))
            ax[subplot_row, subplot_col].set_xticklabels(list(rdms_cummulative_allSubj['targets']), fontdict=None, minor=False, rotation=-45)
            ax[subplot_row, subplot_col].set_xticks(range(len(rdms_cummulative_allSubj[roi])))
            ax[subplot_row, subplot_col].set_yticklabels(list(rdms_cummulative_allSubj['targets']), fontdict=None, minor=False)
            ax[subplot_row, subplot_col].set_yticks(range(len(rdms_cummulative_allSubj[roi])))
            ax[subplot_row, subplot_col].set_title(roi + '; RDM stability: %s +/- %s'
                                                   % ('{:.2f}'.format(np.mean(rdms_stab_allSubj[roi])),
                                                      '{:.2f}'.format(np.std(rdms_stab_allSubj[roi]))))
            
            for i in range(len(rdms_cummulative_allSubj[roi])):
                for j in range(len(rdms_cummulative_allSubj[roi])):
                    if i != j:
                        ax[subplot_row, subplot_col].text(j, i, '{:.2f}'.format(rdms_cummulative_allSubj[roi][i, j]), ha='center', va='center', color='w')
            
        plt.savefig(opj(group_result_dir, 'RDM_space-' + space + '_roi_' + file_suffix + '.png'), dpi=300)
        plt.close()
