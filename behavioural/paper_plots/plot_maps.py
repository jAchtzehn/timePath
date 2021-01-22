from nilearn.plotting import plot_stat_map, plot_glass_brain
from os.path import join as opj
from os.path import abspath
import matplotlib.pyplot as plt
from nilearn.image import load_img, math_img, threshold_img, resample_to_img, new_img_like
from nilearn._utils.niimg import _safe_get_data

import numpy as np

stat_img_fname = abspath('/Users/jachtzehn/Documents/DZNE/timePath/paper/figures/fig3/conjAll_thrSPM.nii')
atlas_img = abspath('/Users/jachtzehn/data/fMRI/atlases/100um/Synthesized_FLASH25_in_MNI_v2_200um.nii.gz')

stat_img = load_img(stat_img_fname)
stat_img_ = _safe_get_data(stat_img, ensure_finite=True)

stat_img__ = new_img_like(stat_img, data=stat_img_)
#stat_img_resampled = resample_to_img(stat_img__, atlas_img_500)
#stat_img_resampled.to_filename('/home/achtzehnj/data/timePath/conjAll_thrSPM_resampled.nii')

plot_stat_map(stat_img, black_bg=False, cut_coords=[42, 45, -13], draw_cross=False, bg_img=atlas_img,
              cmap='plasma', vmin=5.1, threshold=0, vmax=11.28, colorbar=True, dim=0)

plt.savefig('/Users/jachtzehn/Documents/DZNE/timePath/paper/figures/fig3/conj_white_bg.pdf', dpi=500, format='pdf')
plt.close()

# for fname in ['dotsvslumin', 'distvslumin', 'timevslumin', 'conjAll']:
#
#     stat_img_folder = abspath('/Users/jachtzehn/Documents/DZNE/timePath/paper/figures/fig3/')
#     stat_img_fname = opj(stat_img_folder, fname + '_thrSPM.nii')
#
#     plot_glass_brain(stat_img_fname, black_bg=False, cmap='plasma', vmin=5.1, threshold=0, vmax=11.28, colorbar=False)
#     plt.savefig(opj(stat_img_folder, fname + '_gb.png'), dpi=500, format='png')
#     plt.close()