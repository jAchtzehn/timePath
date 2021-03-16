from nilearn.plotting import plot_stat_map, plot_glass_brain
from os.path import join as opj
from os.path import abspath
import matplotlib.pyplot as plt
from nilearn.image import load_img, math_img, threshold_img, resample_to_img, new_img_like
from nilearn._utils.niimg import _safe_get_data

import numpy as np

atlas_img_path = abspath('/Users/jachtzehn/data/fMRI/atlases/100um/Synthesized_FLASH25_in_MNI_v2_200um.nii.gz')
atlas_img = load_img(atlas_img_path)
atlas_img = threshold_img(atlas_img, 9.)



stat_img_fname = abspath('/Users/jachtzehn/data/fMRI/timePath/derivatives/rsa/corrImg_rel-time_irrel-dist_cdist_p.nii.gz')
stat_img = load_img(stat_img_fname)
stat_img__ = threshold_img(stat_img, threshold=-1, mask_img='/Users/jachtzehn/data/fMRI/timePath/derivatives/rsa/corrImg_rel-time_irrel-dist_cdist_p_cluster_mask.nii')

plot_stat_map(stat_img__, black_bg=True, cut_coords=[-39, -63, 1.2], draw_cross=True, bg_img=atlas_img,
              cmap='viridis', vmin=-1, vmax=1, colorbar=True, dim=0)

plt.savefig('/Users/jachtzehn/Desktop/corrImg_rel-time_irrel-dist_200um_thresholded.pdf', dpi=500, format='pdf')
plt.close()


stat_img_fname = abspath('/Users/jachtzehn/data/fMRI/timePath/derivatives/rsa/corrImg_rel-dist_irrel-time_cdist_p.nii.gz')
stat_img = load_img(stat_img_fname)
stat_img__ = threshold_img(stat_img, threshold=-1, mask_img='/Users/jachtzehn/data/fMRI/timePath/derivatives/rsa/corrImg_rel-dist_irrel-time_cdist_p_cluster_mask.nii')

plot_stat_map(stat_img__, black_bg=True, cut_coords=[-27, -78, 44], draw_cross=True, bg_img=atlas_img,
              cmap='viridis', vmin=-1, vmax=1, colorbar=True, dim=0)

plt.savefig('/Users/jachtzehn/Desktop/corrImg_rel-dist_irrel-time_200um_thresholded.pdf', dpi=500, format='pdf')
plt.close()




# for fname in ['dotsvslumin', 'distvslumin', 'timevslumin', 'conjAll']:
#
#     stat_img_folder = abspath('/Users/jachtzehn/Documents/DZNE/timePath/paper/figures/fig3/')
#     stat_img_fname = opj(stat_img_folder, fname + '_thrSPM.nii')
#
#     plot_glass_brain(stat_img_fname, black_bg=False, cmap='plasma', vmin=5.1, threshold=0, vmax=11.28, colorbar=False)
#     plt.savefig(opj(stat_img_folder, fname + '_gb.png'), dpi=500, format='png')
#     plt.close()