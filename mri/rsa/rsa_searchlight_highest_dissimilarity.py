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
from nilearn.image import coord_transform
from utilFunctions import create_html_report

# ------------ File I/O ------------
if platform == 'darwin':
    experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
    experiment_dir = abspath('/mnt/work/achtzehnj/data/')

nilearn_dir = opj(experiment_dir, 'nilearn')
rsa_dir = opj(experiment_dir, 'rsa')

os.system('mkdir -p %s' % rsa_dir)                 # create the mvpa2 folder for data storage

# ------------ options ------------
subjects = [1]     # subjects to calculate
# [4, 6, 7, 9, 12, 14, 17, 19, 23, 24]                                                      # 10 best subjects (TSNR >= 50)
# [1, 2, 3, 5, 8, 10, 11, 13, 15, 16, 18, 20, 21, 22, 25]                                   # the rest
# [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]  # all except 15
conditions = ['time', 'dist', 'dots']
file_suffix = 'wb_similarity'
plot_html = False
group_level = False
save_figs = False
sl_radius = 3
space = 'MNI152NLin2009cAsym'

plot_results = False

os.system('clear')                              # clean up terminal
# if file_suffix is empty, add _
if file_suffix == '':
    file_suffix = '_'

# make group results dir
group_result_dir = opj(rsa_dir, 'group_results_space-' + space)
os.system('mkdir -p %s' % group_result_dir)

rdm_location_img_list = []
for subj in subjects:
    subj_str = 'sub-' + str(subj).zfill(2)      # create the subject string (just for easier handling)

    print('Calculating subject: %s' % subj_str)
    result_dir = opj(rsa_dir, subj_str, 'space-' + space, 'results')
    os.system('mkdir -p %s' % result_dir)                 # create the mvpa2 folder for data storage

    # -- data loading
    data_path = opj(nilearn_dir, subj_str, 'space-' + space, 'betas', subj_str + '_betas_merged.nii.gz')
    behav_path = opj(nilearn_dir, subj_str, 'space-' + space, 'betas', subj_str + '_merged_events.tsv')
    mask_path = opj(nilearn_dir, 'group_masks', 'space-MNI152NLin2009cAsym', 'group_mask_' + 'wb_binarized.nii.gz')
    anat_path = opj(experiment_dir, 'fmriprep', subj_str, 'anat', subj_str + '_T1w_preproc.nii.gz')

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

    for cond in ['time', 'dist', 'dots', 'lumin']:
        if cond not in conditions:
            fmri_data = fmri_data[fmri_data.sa.targets != cond]      # remove cross unwanted conditions from the conditions (targets)

    # -- computation
    # first compute RDMs for all trials
    mtgs = mean_group_sample(['targets'])       # create the mean of all trials
    mtds = mtgs(fmri_data)
    dsm = rsa.PDist(square=False)

    sl = sphere_searchlight(dsm, radius=sl_radius, nproc=40)
    slres = sl(mtds)

    # now compute RDM consistency by cross-validation with chunks
    mtcgs = mean_group_sample(['targets', 'chunks'])
    mtcds = mtcgs(fmri_data)
    dscm = rsa.PDistConsistency()

    sl_cons = sphere_searchlight(dscm, radius=sl_radius, nproc=40)
    slres_cons = sl_cons(mtcds)
    mean_consistency = np.mean(slres_cons, axis=0)

    # -- create .nii.gz
    combinations = ['dots-dist', 'dist-time', 'dots-time']
    for i in range(3):
        niimg = map2nifti(fmri_data, slres.samples[i, :])
        niimg.to_filename(opj(result_dir, subj_str + '_space-' + space + '_' + file_suffix + '_' + combinations[i] + '_searchlight.nii.gz'))

        #niimg = map2nifti(fmri_data, mean_consistency)
        #niimg.to_filename(opj(result_dir, subj_str + '_space-' + space + '_' + file_suffix + '_' + combinations[i] + 'consistency_searchlight.nii.gz'))


    if plot_results:
        # -- plotting
        if roi == 'all':
            fig, ax = plt.subplots(1, 2, figsize=(15, 5))
        else:
            fig, ax = plt.subplots(1, 3, figsize=(20, 5))

        fig.suptitle('Subject %s, mask: %s' % (subj, roi))
        plt.tight_layout(5)
        fig.subplots_adjust(wspace=0.05)

        # plot location
        coords = tuple(coord_transform(mtds.fa.voxel_indices[mean_consistency.argmax()][0],
                                       mtds.fa.voxel_indices[mean_consistency.argmax()][1],
                                       mtds.fa.voxel_indices[mean_consistency.argmax()][2],
                                       fmri_data.a.imgaffine))
        plot_stat_map(opj(result_dir, subj_str + '_space-' + space + '_mask-' + roi + file_suffix + 'RDM_location_searchlight.nii.gz'), bg_img=anat_path,
                      display_mode='ortho', draw_cross=True,
                      cmap='viridis', symmetric_cbar=False,
                      cut_coords=coords, threshold=0, axes=ax[0])
        ax[0].set_title('Location of most stable pattern')

        maxim = mean_consistency.argmax()
        dat = slres.samples[:, mean_consistency.argmax()]
        # plot RDM
        rdm_data = squareform(slres.samples[:, mean_consistency.argmax()])
        rdm_plot = ax[1].imshow(rdm_data, interpolation='nearest')
        div = make_axes_locatable(ax[1])
        cax = div.append_axes("right", size="5%", pad=0.2)
        cbar = plt.colorbar(rdm_plot, cax=cax)
        rdm_plot.set_clim(vmin=0, vmax=2)

        ax[1].set_xticklabels(list(mtds.sa.targets), fontdict=None, minor=False, rotation=-45)
        ax[1].set_xticks(range(len(mtds.sa.targets)))
        ax[1].set_yticklabels(list(mtds.sa.targets), fontdict=None, minor=False)
        ax[1].set_yticks(range(len(mtds.sa.targets)))
        ax[1].set_title('Most stable dissimilarity pattern: {:.2f}'.format(np.mean(mean_consistency[mean_consistency.argmax()])))

        # add text of value to each field
        for i in range(len(rdm_data)):
            for j in range(len(rdm_data)):
                if i != j:
                    ax[1].text(j, i, '{:.2f}'.format(rdm_data[i, j]), ha='center', va='center', color='w')

        rdms_allSubj[roi] = rdms_allSubj[roi] + rdm_data        # save current RDM for group level
        rdms_allSubj['targets'] = mtds.sa.targets               # save condition names for group level

        # plot correlation of other voxels to most stable pattern
        if roi != 'all':
            ax[2].hist(slres_tdsm.samples[0], normed=False, bins=75, color='SkyBlue')
            ax[2].set_ylabel('Number of voxels')
            ax[2].set_xlabel('Target similarity structure correlation')
            ax[2].set_title('Correlation of other searchlight voxels with most stable pattern')

        if save_figs:
            plt.savefig(opj(result_dir, subj_str + '_space-' + space + '_mask-' + roi + file_suffix + 'RDM_location_searchlight.png'), dpi=300)

        plt.show()
        # plt.close()

        rdm_location_img_list.append(opj(result_dir, subj_str + '_space-' + space + '_mask-' + roi + file_suffix + 'RDM_location_searchlight.png'))

if plot_html:
    create_html_report(rdm_location_img_list,
                       opj(group_result_dir, 'RDM_space-' + space + '_mask-' + roi + file_suffix + 'searchlight_location.html'), title='Most stable RDM and location for mask: %s' % roi)

if group_level:

    if (len(masks) % 2) == 0:
        fig, ax = plt.subplots(2, len(masks)/2, figsize=(18, 8))
    else:
        fig, ax = plt.subplots(1, len(masks), figsize=(20, 5), squeeze=False)

    plt.tight_layout(5)
    fig.suptitle('Mean RDMs for most stable patterns', fontsize=20)

    # find max
    max_vals = []
    for roi in masks:
        max_vals.append(np.max(rdms_allSubj[roi] / len(subjects)))

    for roi_enum, roi in enumerate(masks):
        rdms_allSubj[roi] = rdms_allSubj[roi] / len(subjects)               # compute mean

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
        rdm_plot = ax[subplot_row, subplot_col].imshow(rdms_allSubj[roi], interpolation='nearest')
        div = make_axes_locatable(ax[subplot_row, subplot_col])
        cax = div.append_axes("right", size="5%", pad=0.2)
        cbar = plt.colorbar(rdm_plot, cax=cax)
        rdm_plot.set_clim(vmin=0, vmax=np.max(max_vals))
        ax[subplot_row, subplot_col].set_xticklabels(list(rdms_allSubj['targets']), fontdict=None, minor=False, rotation=-45)
        ax[subplot_row, subplot_col].set_xticks(range(len(rdms_allSubj[roi])))
        ax[subplot_row, subplot_col].set_yticklabels(list(rdms_allSubj['targets']), fontdict=None, minor=False)
        ax[subplot_row, subplot_col].set_yticks(range(len(rdms_allSubj[roi])))
        ax[subplot_row, subplot_col].set_title('Mean most stable patterns %s' % roi)

        for i in range(len(rdms_allSubj[roi])):
            for j in range(len(rdms_allSubj[roi])):
                if i != j:
                    ax[subplot_row, subplot_col].text(j, i, '{:.2f}'.format(rdms_allSubj[roi][i, j]), ha='center', va='center', color='w')

    plt.savefig(opj(group_result_dir, 'RDM_space-' + space + file_suffix + 'searchlight_mean.png'), dpi=300)
    plt.close()