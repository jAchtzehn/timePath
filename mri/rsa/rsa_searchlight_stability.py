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
from utilFunctions import create_html_report
import cPickle

# ------------ File I/O ------------
if platform == 'darwin':
    experiment_dir = abspath('/Volumes/Seagate/Backups/16102020/data/timePath')
else:
    experiment_dir = abspath('/home/achtzehnj/data/timePath/')

nilearn_dir = opj(experiment_dir, 'derivatives', 'nilearn')
rsa_dir = opj(experiment_dir, 'derivatives', 'rsa')

os.system('mkdir -p %s' % rsa_dir)                 # create the mvpa2 folder for data storage

# ------------ options ------------
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]     # subjects to plot
#subjects = [1, 2]     # subjects to plot
# [4, 6, 7, 9, 12, 14, 17, 19, 23, 24]                                                                  # 10 best subjects (TSNR >= 50)
# [1, 2, 3, 5, 8, 10, 11, 13, 15, 16, 18, 20, 21, 22, 25]                                               # the rest
# [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]               # all except 15
conditions = ['time', 'dist', 'dots']
masks = ['rois', 'hc_lh', 'hc_rh', 'V5_lh']
# ['ips_lh', 'ips_rh', 'ifg_lh', 'ifg_rh', 'V5_lh', 'V5_rh', 'hc_lh', 'hc_rh', 'insula_lh', 'insula_rh']
# ['ips_lh', 'ips_rh', 'ifg_lh', 'ifg_rh', 'V5_lh', 'V5_rh']
file_suffix = ''
plot_html = False
group_level = True
save_figs = True
sl_radius = 6
spaces = ['MNI152NLin2009cAsym']

os.system('clear')                              # clean up terminal
# if file_suffix is empty, add _
if file_suffix == '':
    file_suffix = '_'

for space in spaces:
    print('Calculating space: %s' % space)
    
    # create arrays for group level
    rdms_cummulative_allSubj = {}
    rdms_allSubj = {}
    max_cons_allSubj = {}
    mean_cons_allSubj = {}

    for roi in masks:
        rdms_cummulative_allSubj[roi] = np.zeros((len(conditions), len(conditions)))
        rdms_allSubj[roi] = []
        max_cons_allSubj[roi] = []
        mean_cons_allSubj[roi] = []

    # make group results dir
    group_result_dir = opj(rsa_dir, 'group_results_space-' + space)
    os.system('mkdir -p %s' % group_result_dir)
    
    for roi in masks:
        print('Calculating ROI: %s' % roi)
        
        rdm_location_img_list = []
        rdm_subjects = []
        for subj in subjects:
            subj_str = 'sub-' + str(subj).zfill(2)      # create the subject string (just for easier handling)
            
            print('Subject: %s' % subj_str)
            result_dir = opj(rsa_dir, subj_str, 'space-' + space, 'results')
            os.system('mkdir -p %s' % result_dir)                 # create the mvpa2 folder for data storage
            
            # -- data loading
            data_path = opj(nilearn_dir, subj_str, 'space-' + space, 'betas', subj_str + '_betas_merged.nii.gz')
            behav_path = opj(nilearn_dir, subj_str, 'space-' + space, 'betas', subj_str + '_merged_events.tsv')
            if roi == 'V5_rh' or roi == 'V5_lh' or roi == 'rois':
                mask_path = opj(nilearn_dir, 'group_masks', 'space-' + space,
                                    'group_mask_' + roi + '_binarized.nii.gz')  # mask filename
            else:
                mask_path = opj(nilearn_dir, subj_str, 'space-' + space, 'masks',
                                    subj_str + '_' + roi + '_mask_binarized.nii.gz')  # mask filename

            anat_path = opj(experiment_dir, 'derivatives', 'templates', 'mni_1mm', '1mm_T1.nii.gz')
            
            # read in behav data
            behav_data = pd.read_csv(behav_path, delimiter='\t')
            targets = behav_data['trial'].tolist()                  # labels are called targets in mvpa
            chunks = behav_data['class_chunk'].tolist()             # chunks are runpairs of 2 runs
            runs = behav_data['run'].tolist()
            
            # create dataset
            if roi != 'brain':
                fmri_data = fmri_dataset(data_path, mask=mask_path, targets=targets, chunks=chunks)
            else:
                fmri_data = fmri_dataset(data_path, targets=targets, chunks=chunks)
            
            fmri_data.sa['runs'] = runs
            fmri_data = remove_nonfinite_features(fmri_data)
            
            # -- preprocessing
            zscore(fmri_data, chunks_attr='chunks')                     # zscoring the data to remove mean
            fmri_data = fmri_data[fmri_data.sa.targets != 'cross']      # remove cross from the conditions (targets)
            
            for cond in ['time', 'dist', 'dots', 'lumin']:
                if cond not in conditions:
                    fmri_data = fmri_data[fmri_data.sa.targets != cond]      # remove cross from the conditions (targets)
                    
            # -- computation
            # first compute RDMs
            mtgs = mean_group_sample(['targets'])       # create the mean of all trials
            mtds = mtgs(fmri_data)
            dsm = rsa.PDist(square=False, center_data=True)

            slres = dsm(mtds)

            # sl = sphere_searchlight(dsm, radius=sl_radius, nproc=24)
            # slres = sl(mtds)
            #
            # # now compute RDM consistency by cross-validation with chunks
            # mtcgs = mean_group_sample(['targets', 'chunks'])
            # mtcds = mtcgs(fmri_data)
            # dscm = rsa.PDistConsistency(center_data=True)
            #
            # sl_stab = sphere_searchlight(dscm, radius=sl_radius, nproc=24)
            # slres_stab = sl_stab(mtcds)
            # mean_consistency = np.mean(slres_stab, axis=0)
            #
            # # compare all the other voxels to the most stable pattern and
            # # find out if they are alike
            # tdsm = rsa.PDistTargetSimilarity(slres.samples[:, mean_consistency.argmax()])
            # sl_tdsm = sphere_searchlight(ChainLearner([tdsm, TransposeMapper()]), radius=sl_radius, nproc=8)
            # slres_tdsm = sl_tdsm(mtds)
            
            # # -- create .nii.gz
            # niimg = map2nifti(fmri_data, mean_consistency)
            # niimg.to_filename(opj(result_dir, subj_str + '_space-' + space + '_mask-' + roi + file_suffix + 'RDM_location_searchlight.nii.gz'))

            # save individual results
            with open(opj(result_dir, subj_str + '_space-' + space + '_mask-' + roi + file_suffix + 'roiData.pkl'), 'wb+') as f:
                cPickle.dump([slres, mtds.sa.targets], f)
                f.close()

            # -- plotting
            fig, ax = plt.subplots(1, 3, figsize=(20, 5))
            fig.suptitle('Subject %s, mask: %s' % (subj, roi))
            plt.tight_layout(5)
            fig.subplots_adjust(wspace=0.05)
            
            # # plot location
            # coords = tuple(coord_transform(mtds.fa.voxel_indices[mean_consistency.argmax()][0],
            #                                mtds.fa.voxel_indices[mean_consistency.argmax()][1],
            #                                mtds.fa.voxel_indices[mean_consistency.argmax()][2],
            #                                fmri_data.a.imgaffine))
            # plot_stat_map(opj(result_dir, subj_str + '_space-' + space + '_mask-' + roi + file_suffix + 'RDM_location_searchlight.nii.gz'), bg_img=anat_path,
            #               display_mode='ortho', draw_cross=True,
            #               cmap='viridis', symmetric_cbar=False,
            #               cut_coords=coords, threshold=0, axes=ax[0])
            # ax[0].set_title('Location of most stable pattern')
            
            # plot RDM
            # reformat so that the matrix has the order time-space-numerosity
            rdm_data_ns = np.zeros(shape=3)
            rdm_max_cons = slres.samples
            rdm_data_ns[0] = rdm_max_cons[1]
            rdm_data_ns[1] = rdm_max_cons[2]
            rdm_data_ns[2] = rdm_max_cons[0]

            #rdm_data = squareform(slres.samples[:, mean_consistency.argmax()])
            rdm_data = squareform(rdm_data_ns)

            # rdm_plot = ax[1].imshow(rdm_data, interpolation='nearest')
            # div = make_axes_locatable(ax[1])
            # cax = div.append_axes("right", size="5%", pad=0.2)
            # cbar = plt.colorbar(rdm_plot, cax=cax)
            # rdm_plot.set_clim(vmin=0, vmax=2)
            #
            # mlabels = ['time', 'space', 'numerosity']
            # ax[1].set_xticklabels(mlabels, fontdict=None, fontsize=18, minor=False, rotation=-45)
            # ax[1].set_xticks(range(len(mlabels)))
            # ax[1].set_yticklabels(mlabels, fontdict=None, fontsize=18, minor=False)
            # ax[1].set_yticks(range(len(mlabels)))
            # ax[1].set_title('Most stable dissimilarity pattern: {:.2f}'.format(np.mean(mean_consistency[mean_consistency.argmax()])))
            #
            # # add text of value to each field
            # for i in range(len(rdm_data)):
            #     for j in range(len(rdm_data)):
            #         if i != j:
            #             ax[1].text(j, i, '{:.2f}'.format(rdm_data[i, j]), ha='center', va='center', color='w', fontsize=18)
            
            rdms_cummulative_allSubj[roi] = rdms_cummulative_allSubj[roi] + rdm_data        # save current RDM for group level
            #max_cons_allSubj[roi].append(mean_consistency[mean_consistency.argmax()])    # save current rdm consistency for group level
            #mean_cons_allSubj[roi].append(np.mean(mean_consistency))
            rdms_cummulative_allSubj['targets'] = mtds.sa.targets               # save condition names for group level
            rdms_allSubj[roi].append(rdm_data)
            #
            # # plot correlation of other voxels to most stable pattern
            # ax[2].hist(slres_tdsm.samples[0], normed=False, bins=75, color='SkyBlue')
            # ax[2].set_ylabel('Number of voxels')
            # ax[2].set_xlabel('Target similarity structure correlation')
            # ax[2].set_title('Correlation of other searchlight voxels with most stable pattern')
            #
            # if save_figs:
            #     plt.savefig(opj(result_dir, subj_str + '_space-' + space + '_mask-' + roi + file_suffix + 'RDM_location_searchlight_new.png'), dpi=300)
            
            plt.close()
            
            # rdm_location_img_list.append(opj(result_dir, subj_str + '_space-' + space + '_mask-' + roi + file_suffix + 'RDM_location_searchlight.png'))

        if plot_html:
            create_html_report(rdm_location_img_list,
                               opj(group_result_dir, 'RDM_space-' + space + '_mask-' + roi + file_suffix + 'searchlight_location.html'), title='Most stable RDM and location for mask: %s' % roi)

    if group_level:

        # save
        with open(opj(group_result_dir, 'RDM_space-' + space + file_suffix + 'searchlight_mean.pkl'), 'wb+') as f:
            cPickle.dump([rdms_allSubj, rdms_cummulative_allSubj], f)
            f.close()

        print('Group level...')

        # write .csv file with results
        pd_subjId, pd_roi, pd_cond1, pd_cond2, pd_dm, pd_max_cons, pd_mean_cons = [], [], [], [], [], [], []
        for s, subjId in enumerate(subjects):
            for roi in masks:
                for i, cond_comb in enumerate([['dots', 'dist'], ['dist', 'time'], ['dots', 'time']]):
                    pd_subjId.append(subjId)
                    pd_roi.append(roi)
                    pd_cond1.append(cond_comb[0])
                    pd_cond2.append(cond_comb[1])

                    if i == 0:
                        pd_dm.append(rdms_allSubj[roi][s][2, 1])
                    elif i == 1:
                        pd_dm.append(rdms_allSubj[roi][s][1, 0])
                    else:
                        pd_dm.append(rdms_allSubj[roi][s][2, 0])

                    #pd_max_cons.append(max_cons_allSubj[roi][s])
                    #pd_mean_cons.append(mean_cons_allSubj[roi][s])

        df = pd.DataFrame({'subject': pd_subjId,
                           'roi': pd_roi,
                           'condition_1': pd_cond1,
                           'condition_2': pd_cond2,
                           'dissimilarity': pd_dm})

        df.to_csv(opj(group_result_dir, 'RDM_space-MNI152NLin2009cAsym_searchlight_mean.csv'), sep='\t',
                  columns=['subject', 'roi', 'condition_1', 'condition_2', 'dissimilarity'],
                  index=False)

        if (len(masks) % 2) == 0:
            fig, ax = plt.subplots(2, len(masks)/2, figsize=(20, 15))
        else:
            fig, ax = plt.subplots(1, len(masks), figsize=(20, 5), squeeze=False)

        fig.subplots_adjust(hspace=.5, wspace=.6)
        #fig.suptitle('Mean RDMs for most stable patterns', fontsize=20)
        
        # find max
        max_vals = []
        for roi in masks:
            max_vals.append(np.max(rdms_cummulative_allSubj[roi] / len(subjects)))
        
        for roi_enum, roi in enumerate(masks):
            rdms_cummulative_allSubj[roi] = rdms_cummulative_allSubj[roi] / len(subjects)               # compute mean
            
            if (len(masks) % 2) == 0:
                if roi_enum < len(masks)/2:
                    subplot_row = 0
                    subplot_col = roi_enum
                else:
                    subplot_row = 1
                    subplot_col = roi_enum - len(masks)/2
            else:
                subplot_row = 0
                subplot_col = roi_enum
        
            # plotting
            rdm_plot = ax[subplot_row, subplot_col].imshow(rdms_cummulative_allSubj[roi], interpolation='nearest')
            div = make_axes_locatable(ax[subplot_row, subplot_col])
            cax = div.append_axes("right", size="5%", pad=0.2)
            cbar = plt.colorbar(rdm_plot, cax=cax)
            cbar.ax.tick_params(labelsize=16)
            rdm_plot.set_clim(vmin=0, vmax=np.max(max_vals))

            mlabels = ['time', 'space', 'numerosity']
            ax[subplot_row, subplot_col].set_xticklabels(mlabels, fontdict=None, minor=False, rotation=-45, fontsize=18)
            ax[subplot_row, subplot_col].set_xticks(range(len(mlabels)))
            ax[subplot_row, subplot_col].set_yticklabels(mlabels, fontdict=None, minor=False, fontsize=18)
            ax[subplot_row, subplot_col].set_yticks(range(len(mlabels)))
            ax[subplot_row, subplot_col].set_title('{}'.format(roi, np.mean(max_cons_allSubj[roi])))
            
            for i in range(len(rdms_cummulative_allSubj[roi])):
                for j in range(len(rdms_cummulative_allSubj[roi])):
                    if i != j:
                        ax[subplot_row, subplot_col].text(j, i, '{:.2f}'.format(rdms_cummulative_allSubj[roi][i, j]), ha='center', va='center', color='black', fontsize=18)
            
        plt.savefig(opj(group_result_dir, 'RDM_space-' + space + file_suffix + 'searchlight_mean.pdf'), format='pdf', dpi=300)
        plt.close()