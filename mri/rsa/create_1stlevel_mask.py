
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
for contrast in [11, 13, 15]:

	conFolder = '/Volumes/Seagate/Backups/16102020/data/timePath/derivatives/output_nipype/2ndlevel/space-MNI152NLin2009cAsym/contrast-00' + str(contrast)

	fslcmd = 'fslmaths ' + opj(conFolder, 'spmT_0001_thr.nii') + ' -nan ' + opj(conFolder, 'spmT_0001_thr_nan_zero.nii')
	conImgs[str(contrast)] = opj(conFolder, 'spmT_0001_thr_nan_zero.nii')
	os.system(fslcmd)
	print('Running {}'.format(fslcmd))


fslcmd = 'fslmaths ' + conImgs['11'] + ' -add ' + conImgs['13'] + ' -add ' + conImgs['15'] + ' /Users/jachtzehn/Desktop/mask.nii'
print('Running {}'.format(fslcmd))
os.system(fslcmd)

fslcmd = 'mri_binarize --i /Users/jachtzehn/Desktop/mask.nii.gz ' + '--o /Users/jachtzehn/Desktop/mask_thr.nii.gz --min 0.01'

os.system(fslcmd)