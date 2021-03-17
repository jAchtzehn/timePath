import os
from os.path import join as opj
from os.path import abspath
from sys import platform
import pandas as pd
import itertools
from nilearn.image import load_img, index_img, mean_img, new_img_like, threshold_img, coord_transform
from nilearn.glm import threshold_stats_img
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from nilearn.input_data import NiftiMasker
import scipy.stats as stats
from tqdm import tqdm
import matplotlib.pyplot as plt
from nilearn.plotting import plot_glass_brain
import seaborn as sns
from scipy import ndimage
from scipy.spatial.distance import squareform


# ------------ File I/O ------------
#experiment_dir = abspath('/home/achtzehnj/data/timePath/')
experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath')

rsa_dir = opj(experiment_dir, 'derivatives', 'rsa')
behav_file = opj(experiment_dir, 'derivatives', 'behavioural', 'pse_data_cross_dim_individual_norm.tsv')
behav_data = pd.read_csv(behav_file, delimiter='\t')
mask_img = load_img(opj(experiment_dir, 'derivatives', 'rsa', 'group_mask_wb_binarized.nii.gz'))

# ------------ options ------------
subjects = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
conditions = ['time', 'dist', 'dots']

masker = NiftiMasker(mask_img=mask_img, standardize=False, detrend=False,
                     memory=abspath('/Users/jachtzehn/data/fMRI/nilearn_cache'), memory_level=1)

for x, condition_pair in enumerate([['time', 'dist'], ['dist', 'time']]):

    condition_str = 'rel-' + condition_pair[0] + '_irrel-' + condition_pair[1]

    # first plot glassbrain plots

    corr_img = load_img(opj(rsa_dir, 'corrImg_' + condition_str + '_cdist_p.nii.gz'))
    corr_values = masker.fit_transform(corr_img)

    corr_mask_img = load_img(opj(rsa_dir,  'corrImg_' + condition_str + '_cdist_p_cluster_mask.nii'))

    corr_img_masked = threshold_img(corr_img, mask_img=corr_mask_img, threshold=-1)

    plot_glass_brain(corr_img_masked, colorbar=True, cmap='viridis', vmin=-1, vmax=1, plot_abs=False)

    plt.savefig(opj(rsa_dir, 'corrImg_' + condition_str + '_glassBrain_thresholded.pdf'), dpi=400, format='pdf')
    plt.close()

    # plot RDM of maximum correlation voxel
    max_vx_idx = np.abs(corr_values).argmax()
    max_vx_loc = ndimage.maximum_position(corr_img_masked.get_data())
    max_mni_loc = coord_transform(max_vx_loc[0], max_vx_loc[1], max_vx_loc[2], corr_img_masked.affine)
    min_vx_loc = ndimage.minimum_position(corr_img_masked.get_data())
    min_mni_loc = coord_transform(min_vx_loc[0], min_vx_loc[1], min_vx_loc[2], corr_img_masked.affine)

    print(max_mni_loc)
    print(min_mni_loc)

    # correlate
    crossDim = behav_data['pse_diff'][(behav_data.condition_rel == condition_pair[0]) &
                                      (behav_data.condition_irrel == condition_pair[1])].values

    rdms_mean = np.zeros(shape=(3, 24))
    rdms_crossDim = np.zeros(shape=(1, 23))
    for i, sub in enumerate(subjects):
        sub_str = 'sub-' + str(sub).zfill(2)

        rdm_img = load_img(opj(rsa_dir, sub_str, 'space-MNI152NLin2009cAsym', 'results', sub_str + '_space-MNI152NLin2009cAsym_cdist_wb_rdm_values.nii.gz'))

        if x == 1:
            rdm_value = rdm_img.get_data()[max_vx_loc[0], max_vx_loc[1], max_vx_loc[2]]
        else:
            rdm_value = rdm_img.get_data()[min_vx_loc[0], min_vx_loc[1], min_vx_loc[2]]

        rdms_crossDim[0, i] = rdm_value[0]
        rdms_mean[:, i] = rdm_value

    [vx_corr_r, vx_corr_p] = stats.spearmanr(crossDim, rdms_crossDim.reshape(-1))
    print('{}, {}'.format(vx_corr_r, vx_corr_p))

    # plot RDM
    # plotting
    mean_rdm = np.mean(rdms_mean, axis=1)

    mean_rdm_sq = squareform(mean_rdm)

    # fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    #
    # rdm_plot = ax.imshow(mean_rdm_sq, interpolation='nearest')
    # div = make_axes_locatable(ax)
    # cax = div.append_axes("right", size="5%", pad=0.2)
    # cbar = plt.colorbar(rdm_plot, cax=cax)
    # cbar.ax.tick_params(labelsize=16)
    # rdm_plot.set_clim(vmin=0, vmax=np.max(1.6))
    #
    # mlabels = ['time', 'space', 'numerosity']
    # ax.set_xticklabels(mlabels, fontdict=None, minor=False, rotation=-45, fontsize=18)
    # ax.set_xticks(range(len(mlabels)))
    # ax.set_yticklabels(mlabels, fontdict=None, minor=False, fontsize=18)
    # ax.set_yticks(range(len(mlabels)))
    #
    # for i in range(len(mean_rdm_sq)):
    #     for j in range(len(mean_rdm_sq)):
    #         if i != j:
    #             ax.text(j, i, '{:.2f}'.format(mean_rdm_sq[i, j]),
    #                                               ha='center', va='center', color='black', fontsize=18)
    #
    # plt.savefig(opj(rsa_dir, 'corrImg_rel-{}_irrel-{}_corr_plot_max_rdm.pdf'.format(condition_pair[0], condition_pair[1])),
    #                format='pdf')
    # plt.close()


    # sns.set_theme(style="darkgrid")
    #
    # sns.regplot(x=rdms_crossDim.reshape(-1), y=crossDim)
    # plt.xlabel('RDM')
    # plt.ylabel('CrossDim')
    # plt.title('rel-{}_irrel-{}'.format(condition_pair[0], condition_pair[1]))
    # plt.savefig(opj(rsa_dir, 'corrImg_rel-{}_irrel-{}_corr_plot_max.pdf'.format(condition_pair[0], condition_pair[1])),
    #             format='pdf')
    # plt.close()