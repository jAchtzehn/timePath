
from os.path import join as opj
from os.path import abspath
from sys import platform
import os
from nilearn.image import load_img, index_img, mean_img, threshold_img


# ------------ File I/O ------------
if platform == 'darwin':
    experiment_dir = abspath('/Volumes/Seagate/Backups/16102020/data/timePath')
else:
    experiment_dir = abspath('/home/achtzehnj/data/timePath/')

nilearn_dir = opj(experiment_dir, 'derivatives', 'nilearn')
rsa_dir = opj(experiment_dir, 'derivatives', 'rsa')


conImgs = {}
conFolder = abspath('/Users/jachtzehn/Documents/DZNE/timePath/paper/figures/fig3/')

for contrast in ['dotsvslumin', 'distvslumin', 'timevslumin']:

	fslcmd = 'fslmaths ' + opj(conFolder, contrast + '_thrSPM.nii') + ' -nan ' + opj(conFolder, contrast + '_thrSPM_nan_zero.nii')
	conImgs[contrast] = opj(conFolder, contrast + '_thrSPM_nan_zero.nii')
	os.system(fslcmd)
	print('Running {}'.format(fslcmd))


fslcmd = 'fslmaths ' + conImgs['dotsvslumin'] + ' -add ' + conImgs['distvslumin'] + ' -add ' + conImgs['timevslumin'] + ' /Users/jachtzehn/Desktop/mask.nii'
print('Running {}'.format(fslcmd))
os.system(fslcmd)

fslcmd = 'mri_binarize --i /Users/jachtzehn/Desktop/mask.nii.gz ' + '--o /Users/jachtzehn/Desktop/mask_thr.nii.gz --min 0.01'

os.system(fslcmd)