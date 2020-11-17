import os
from os.path import join as opj
from os.path import abspath
from sys import platform

import base64
import datetime
import numpy as np

from nipype import config, logging
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.interfaces.spm import Threshold

from nilearn.plotting import plot_glass_brain
from nilearn.image import coord_transform
from scipy import ndimage
from nilearn.image import load_img
import matplotlib as mpl

if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/mnt/work/achtzehnj/data')

working_dir = 'working_nipype'                                          # working directory for nipype
output_dir = 'output_nipype'                                            # where should the relevant output be stored?
results_dir = '1stLevel'

print('Deleting work folder %s...' % opj(experiment_dir, working_dir))
os.system('rm -rf %s'%opj(experiment_dir, working_dir))                 # del working dir for clean s

log_dir_timestamp = datetime.datetime.now().strftime("%y-%m-%d-%H-%M-%S")
log_dir = opj(experiment_dir, 'logs_nipype', results_dir + '_' + log_dir_timestamp)
os.makedirs(log_dir)    # make output dir for current workflow

condition_names = ['time', 'dist', 'lumin', 'dots', 'icon', 'comp']

contrast_list = [['1-trial-avg', 'T', condition_names, [0.25, 0.25, 0.25, 0.25, 0, 0]],
                 ['2-TDD-avg', 'T', condition_names, [1/3, 1/3, 0, 1/3, 0, 0]],
                 ['3-time', 'T', condition_names, [1, 0, 0, 0, 0, 0]],
                 ['4-dist', 'T', condition_names, [0, 1, 0, 0, 0, 0]],
                 ['5-lumin', 'T', condition_names, [0, 0, 1, 0, 0, 0]],
                 ['6-dots', 'T', condition_names, [0, 0, 0, 1, 0, 0]],
                 ['7-icon', 'T', condition_names, [0, 0, 0, 0, 1, 0]],
                 ['8-comp', 'T', condition_names, [0, 0, 0, 0, 0, 1]],
                 ['9-TDD>lumin', 'T', condition_names, [1/3., 1/3., -1, 1/3., 0, 0]],
                 ['10-lumin>TDD', 'T', condition_names, [-1/3., -1/3., 1, -1/3., 0, 0]],
                 ['11-time>lumin', 'T', condition_names, [1, 0, -1, 0, 0, 0]],
                 ['12-lumin>time', 'T', condition_names, [-1, 0, 1, 0, 0, 0]],
                 ['13-dist>lumin', 'T', condition_names, [0, 1, -1, 0, 0, 0]],
                 ['14-lumin>dist', 'T', condition_names, [0, -1, 1, 0, 0, 0]],
                 ['15-dots>lumin', 'T', condition_names, [0, 0, -1, 1, 0, 0]],
                 ['16-lumin>dots', 'T', condition_names, [0, 0, 1, -1, 0, 0]],
                 ['17-time>others', 'T', condition_names, [1, -0.5, 0, -0.5, 0, 0]],
                 ['18-dist>others', 'T', condition_names, [-0.5, 1, 0, -0.5, 0, 0]],
                 ['19-dots>others', 'T', condition_names, [-0.5, -0.5, 0, 1, 0, 0]],
                 ['20-time>dist', 'T', condition_names, [1, -1, 0, 0, 0, 0]],
                 ['21-dist>time', 'T', condition_names, [-1, 1, 0, 0, 0, 0]],
                 ['22-time>dots', 'T', condition_names, [1, 0, 0, -1, 0, 0]],
                 ['23-dots>time', 'T', condition_names, [-1, 0, 0, 1, 0, 0]],
                 ['24-dist>dots', 'T', condition_names, [0, 1, 0, -1, 0, 0]],
                 ['25-dots>dist', 'T', condition_names, [0, -1, 0, 1, 0, 0]],
                 ]

contrasts_ids = np.add(list(range(len(contrast_list))), 1).tolist()

# ---- options ----
nr_procs = 4
subjects = range(1, 2)
cluster_size = 0
fwe_correction = True
fwe_p = 0.05
topo_fdr = False
run_thresh = True
output_space = ['MNI152NLin2009cAsym']                  # T1w and/or MNI152NLin2009cAsym are possible
stop_on_first_crash = 'True'


config_dict = {'execution': {'crashdump_dir': log_dir,
                             'poll_sleep_duration': 2,
                             'crashfile_format': 'txt'} ,
               'logging': {'util_level': 'DEBUG',
                         'log_directory':  log_dir,
                         'log_to_file': 'true',
                         'log_size': 95400000},
               'monitoring': {'enabled': 'true',
                              'summary_file': log_dir}}
config.update_config(config_dict)
logging.update_logging(config)

n_threshold = Node(Threshold(extent_threshold=cluster_size, use_fwe_correction=fwe_correction,
                             height_threshold=fwe_p,  height_threshold_type='p-value', use_topo_fdr=topo_fdr),
                   name='n_threshold')

#contrast_list = [str(x).zfill(4) for x in contrasts_ids]

n_infosource = Node(IdentityInterface(fields=['subj_id', 'output_space', 'contrast_id']),
					   name="n_infosource")
n_infosource.iterables = [('subj_id', subjects),
						  ('output_space', output_space),
						  ('contrast_id', contrasts_ids)]

templates = {'mat': opj(output_dir, results_dir, 'sub-' + '{subj_id:02d}', 'space-{output_space}', 'SPM.mat'),
			 't_img': opj(output_dir, results_dir, 'sub-' + '{subj_id:02d}'.zfill(2), 'space-{output_space}',
			              'spmT_' + '{contrast_id:04d}'.zfill(4) + '.nii')}

n_selectfiles = Node(SelectFiles(templates,
							   base_directory=experiment_dir,
							   sort_filelist=True),
				   name="n_selectfiles")

# datasink
n_datasink = Node(DataSink(base_directory=experiment_dir, container=output_dir),
					   name="n_datasink")

substitutions = []
subjFolders = [('_contrast_id_%s_output_space_%s_subj_id_%s' % (i, o, sub), 'sub-%02d/space-%s' % (sub, o))
               for i in contrasts_ids
			   for o in output_space
			   for sub in subjects]
substitutions.extend(subjFolders)
n_datasink.inputs.substitutions = substitutions

# ------------ Connect workflows ------------
w_l1analysis = Workflow(base_dir=opj(experiment_dir, working_dir), name='w_threshold')
w_l1analysis.config['execution'] = {'stop_on_first_crash': stop_on_first_crash}

w_l1analysis.connect([(n_infosource, n_selectfiles, [('subj_id', 'subj_id'),
													 ('output_space', 'output_space'),
													 ('contrast_id', 'contrast_id')]),
                      (n_selectfiles, n_threshold, [('mat', 'spm_mat_file'),
                                                     ('t_img', 'stat_image')]),
                      (n_infosource, n_threshold, [('contrast_id', 'contrast_index')]),
                      (n_threshold, n_datasink, [('thresholded_map', results_dir + '.@T')])
					  ])

if run_thresh:
	w_l1analysis.run('MultiProc', plugin_args={'n_procs': nr_procs})
	os.system('rm -rf %s'%opj(experiment_dir, working_dir))                 # del working dir after workflow is finished

# plot png images
for subj_id in subjects:
	subj = 'sub-' + str(subj_id).zfill(2)

	subjDir = opj(experiment_dir, output_dir, results_dir, subj, 'space-' + 'MNI152NLin2009cAsym')

	for i, contr in enumerate(contrast_list):

		img = load_img(opj(subjDir, 'spm' + contr[1] + '_' + str(i + 1).zfill(4) + '_thr.nii'))
		img_data = img.get_data()
		max_val = ndimage.maximum(img_data)
		min_val = ndimage.maximum(img_data)

		display = plot_glass_brain(img, title=subj, display_mode='lyrz',
		                           colorbar=True, threshold=0, symmetric_cbar=True, vmin=0, vmax=max_val, cmap='plasma')

		display.savefig(opj(subjDir, contr[0] + '.png'), dpi=200)
		display.close()
		print('Written image: ' + str(opj(subjDir, contr[0] + '.png')))

# write a html file for each contrast
# for contr in contrasts:
# 	f = open(opj(experiment_dir, output_dir, results_dir, contr[0] + '_space-' + 'MNI152NLin2009cAsym' + '.html'), 'w')
#
# 	f.write('<html>\n')
# 	f.write('<h1 style="font-family:helvetica;">')
# 	f.write('<head><title>Contrast: ' + contr[0] + ', ' + contr[1] + '</title></head>\n')
# 	f.write('<body><p><font size="14">Contrast: ' + contr[0] + ', ' + contr[1] + '</font></p></body>\n')
# 	f.write('<hr>')
#
# 	for subj in subj_list:
# 		subjDir = opj(output_dir, subj, 'space-' + 'MNI152NLin2009cAsym')
# 		imagePath = opj(subjDir, contr[0] + '.png')
#
# 		image = base64.b64encode(open(imagePath, 'rb').read()).decode().replace('\n', '')
# 		f.write('<img src="data:image/png;base64,{0}">'.format(image))
# 		f.write('\n')
#
# 	f.write('</h1>')
# 	f.write('</html>\n')
# 	f.close()
