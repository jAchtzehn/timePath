from nilearn.image import mean_img, threshold_img, math_img
from os.path import join as opj
from os.path import abspath
import os


datapath = abspath('/home/achtzehnj/data/timePath/derivatives/nilearn/')
subjects = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]

mask_files = []
for subj in subjects:
    sub_str = 'sub-' + str(subj).zfill(2)
    sub_folder = opj(datapath, sub_str, 'space-MNI152NLin2009cAsym', 'masks')
    mask_files.append(opj(sub_folder, sub_str + '_ifg_mask_binarized.nii.gz'))

avg_image = mean_img(mask_files)
avg_image_thr = threshold_img(avg_image, 0.5)
avg_image_thr_bin = math_img('img >= 0.35', img=avg_image_thr)
avg_image.to_filename('/home/achtzehnj/data/timePath/avg_ifg.nii')
avg_image_thr.to_filename('/home/achtzehnj/data/timePath/avg_ifg_thresholded.nii')
avg_image_thr_bin.to_filename('/home/achtzehnj/data/timePath/avg_ifg_binarized.nii')