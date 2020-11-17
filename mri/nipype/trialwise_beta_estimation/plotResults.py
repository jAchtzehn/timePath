import os
from os.path import join as opj
from os.path import abspath
from nilearn.plotting import plot_stat_map
from nilearn.plotting import plot_glass_brain
from sys import platform
from util_functions import createSubjectlist

import matplotlib.pyplot as plt
from nilearn.image import smooth_img as smooth

if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/media/sf_data/fMRI/timePath/')

subjects = range(1, 26)
contrasts = [ ['T', 'average'], ['T', 'time'], ['T', 'dist'], ['T', 'lumin'], ['T', 'dots'],
              ['T', 'time_vs_lumin'], ['T', 'dist_vs_lumin'], ['T', 'dots_vs_lumin'],
              ['F', 'activation'], ['F', 'differences']]


# ---- options ----
plot_threshold = 4.8
plotMNI = True
plotT1 = False

output_dir = opj(experiment_dir, 'output_nipype', '1stlevel')
subj_list = createSubjectlist(subjects)     # create a list with subjects in the format of ['sub-01', 'sub-02', ...]

# plot MNI
if plotMNI:
	
	for subj in subj_list:
		subjDir = opj(output_dir, subj, 'space-' + 'MNI152NLin2009cAsym')
	
		for i, contr in enumerate(contrasts):
			
			plot_glass_brain(opj(subjDir, 'spm' + contr[0] + '_' + str(i + 1).zfill(4) + '.nii'), title=subj, display_mode='lyrz', output_file=opj(subjDir, contr[1] + '.png'),
								   colorbar=True, plot_abs=False, threshold=plot_threshold, cmap='viridis')
	
			print('Written image: ' + str(opj(subjDir, contr[1] + '.png')))
	
	# write a html file for each contrast
	for contr in contrasts:
		f = open(opj(output_dir, contr[1] + '_space-' + 'MNI152NLin2009cAsym' + '.html'), 'w')
	
		f.write('<html>\n')
		f.write('<h1 style="font-family:helvetica;">')
		f.write('<head><title>Contrast: ' + contr[1] + ', ' + contr[0] + '</title></head>\n')
		f.write('<body><p><font size="14">Contrast: ' + contr[1] + ', ' + contr[0] + '</font></p></body>\n')
		f.write('<hr>')
	
		for subj in subj_list:
			subjDir = opj(output_dir, subj, 'space-' + 'MNI152NLin2009cAsym')
			imagePath = opj(subjDir, contr[1] + '.png')
	
			image = open(imagePath, 'rb').read().encode('base64').replace('\n','')
			f.write('<img src="data:image/png;base64,{0}">'.format(image))
			f.write('\n')
		
		f.write('</h1>')
		f.write('</html>\n')
		f.close()

# plot T1w
if plotT1:
	
	for subj in subj_list:

		anat_img = opj(experiment_dir, 'fmriprep', subj, 'anat', subj + '_T1w_preproc.nii.gz')
		subjDir = opj(output_dir, subj, 'space-' + 'T1w')

		for i, contr in enumerate(contrasts):
			plot_x = plot_stat_map(opj(subjDir, 'spm' + contr[0] + '_' + str(i + 1).zfill(4) + '.nii'), output_file=opj(subjDir, contr[1] + '_x.png'), title='',
								bg_img=anat_img, threshold=plot_threshold, display_mode='x', cut_coords=7, dim=-1, cmap='plasma', draw_cross=True)
			plot_y = plot_stat_map(opj(subjDir, 'spm' + contr[0] + '_' + str(i + 1).zfill(4) + '.nii'), output_file=opj(subjDir, contr[1] + '_y.png'), title='',
								bg_img=anat_img, threshold=plot_threshold, display_mode='y', cut_coords=7, dim=-1, cmap='plasma', draw_cross=True)
			plot_z = plot_stat_map(opj(subjDir, 'spm' + contr[0] + '_' + str(i + 1).zfill(4) + '.nii'), output_file=opj(subjDir, contr[1] + '_z.png'), title='',
								bg_img=anat_img, threshold=plot_threshold, display_mode='z', cut_coords=7, dim=-1, cmap='plasma', draw_cross=True)

		print('Written xyz images: ' + str(subjDir))
	
	# write a html file for each contrast
	for contr in contrasts:
		f = open(opj(output_dir, contr[1] + '_space-' + 'T1w' + '.html'), 'w')
		
		f.write('<html>\n')
		f.write('<head><title>Contrast: ' + contr[1] + ', ' + contr[0] + '</title></head>\n')
		f.write('<body><p><font size="14">Contrast: ' + contr[1] + ', ' + contr[0] + '</font></p></body>\n')
		f.write('<hr>')
		
		for subj in subj_list:
			subjDir = opj(output_dir, subj, 'space-' + 'T1w')
			
			f.write('<body><p><font size="14">Subject: ' + subj + '</font></p></body>\n')
			f.write('<hr>')
			
			for cut_coords in ['x', 'y', 'z']:
				imagePath = opj(subjDir, contr[1] + '_' + cut_coords + '.png')
				
				image = open(imagePath, 'rb').read().encode('base64').replace('\n', '')
				f.write('<img src="data:image/png;base64,{0}">'.format(image))
				f.write('\n')
		
		f.write('</html>\n')
		f.close()
