import os
from os.path import join as opj
from os.path import abspath
from sys import platform
import datetime

from nipype import config, logging
from nipype.interfaces.spm import (OneSampleTTestDesign, EstimateModel,
                                   EstimateContrast, Threshold)
from nipype.interfaces.utility import IdentityInterface
from nipype.pipeline.engine import Workflow, Node
from nipype.interfaces.io import SelectFiles, DataSink

from util_functions import createSubjectlist

# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
else:
	experiment_dir = abspath('/mnt/work/achtzehnj/data')

working_dir = 'working_nipype'          # working directory for nipype
output_dir = 'output_nipype'            # where should the relevant output be stored?
first_level_dir = '1stlevel'            # parent dir of 1st level results
results_dir = '2ndlevel_higher_cluster' # output of second level results

print('Deleting work folder %s...' % opj(experiment_dir, working_dir))
os.system('rm -rf %s' % opj(experiment_dir, working_dir))  # del working dir for clean start

log_dir_timestamp = datetime.datetime.now().strftime("%y-%m-%d-%H-%M-%S")
log_dir = opj(experiment_dir, 'logs_nipype', results_dir + '_' + log_dir_timestamp)

# create log folders for logging if they have been deleted
os.makedirs(log_dir)  # make output dir for current workflow

# check if output_dir already exists and if so if it should be deleted
try:
	if results_dir in os.listdir(opj(experiment_dir, output_dir)):
		print('Output dir %s already exists' % opj(experiment_dir, output_dir, results_dir))
		rmOutput = input('Should it be deleted before proceeding? (y/n): ')
		answerGiven = False
		
		while not answerGiven:
			if rmOutput == 'y':
				os.system('rm -rf %s' % opj(experiment_dir, output_dir, results_dir))
				print('Deleted folder %s' % opj(experiment_dir, output_dir, results_dir))
				answerGiven = True
			elif rmOutput == 'n':
				print('Not deleting output folder %s' % opj(experiment_dir, output_dir, results_dir))
				answerGiven = True
			else:
				print('Please enter (y/n)')
				answerGiven = False
except OSError:
	print('Output dir %s does not exist, proceeding...' % opj(experiment_dir, output_dir))

# ------------ Options ------------
output_space = ['MNI152NLin2009cAsym']  # T1w and/or MNI152NLin2009cAsym are possible
stop_on_first_crash = 'True'
nr_procs = 4
contrasts = range(2, 3)

cluster_size = 60
fwe_correction = True
fwe_p = 0.05
topo_fdr = False

# ------------ nipype config ------------
config_dict = {'execution': {'crashdump_dir': log_dir, 'poll_sleep_duration': 2, 'crashfile_format': 'txt'},
               'logging': {'util_level': 'DEBUG',
                           'log_directory': log_dir,
                           'log_to_file': 'true',
                           'log_size': 95400000},
               'monitoring': {'enabled': 'true',
                              'summary_file': log_dir}}
config.update_config(config_dict)
logging.update_logging(config)

# ------------ contrasts ------------
contrast_ids = [('%04d' % x) for x in contrasts]
second_level_contrast = ['Group', 'T', ['mean'], [1]]

# ------------ Node definition ------------
# One Sample T-Test Design - creates one sample T-Test Design
n_onesamplettestdes = Node(OneSampleTTestDesign(),
                           name="n_onesampttestdes")

# EstimateModel - estimate the parameters of the model
n_level2estimate = Node(EstimateModel(estimation_method={'Classical': 1}),
                        name="n_level2estimate")

# EstimateContrast - estimates simple group contrast
n_level2conestimate = Node(EstimateContrast(group_contrast=True),
                           name="n_level2conestimate")
n_level2conestimate.inputs.contrasts = [second_level_contrast]

# Threshold - thresholds contrasts
n_level2thresh = Node(Threshold(contrast_index=1,
                                use_topo_fdr=topo_fdr,
                                use_fwe_correction=fwe_correction,
                                extent_threshold=cluster_size,
                                height_threshold=fwe_p,
                                height_threshold_type='p-value',
                                extent_fdr_p_threshold=0.05),
                      name="level2thresh")

# --- IO stream

# infosource
n_infosource = Node(IdentityInterface(fields=['output_space', 'contrast_id']),
                    name="n_infosource")
n_infosource.iterables = [('output_space', output_space),
                          ('contrast_id', contrast_ids)]

templates = {'cons': opj('output_nipype', first_level_dir, 'sub-*',
                         'space-{output_space}', 'con_{contrast_id}.nii')}

n_selectfiles = Node(SelectFiles(templates,
                                 base_directory=experiment_dir,
                                 sort_filelist=True),
                     name="n_selectfiles")

# datasink
n_datasink = Node(DataSink(base_directory=experiment_dir, container=output_dir),
                  name="n_datasink")

substitutions = []
subjFolders = [('_contrast_id_%s_output_space_%s' % (con.zfill(4), o), 'space-%s/contrast-%s' % (o, con))
               for o in output_space
               for con in contrast_ids]

substitutions.extend(subjFolders)
n_datasink.inputs.substitutions = substitutions

# ------------ Connect workflows ------------
w_l2analysis = Workflow(base_dir=opj(experiment_dir, working_dir), name='w_l2analysis')
w_l2analysis.config['execution'] = {'stop_on_first_crash': stop_on_first_crash}

w_l2analysis.connect([(n_infosource, n_selectfiles, [('output_space', 'output_space'),
                                                     ('contrast_id', 'contrast_id')]),
                      (n_selectfiles, n_onesamplettestdes, [('cons', 'in_files')]),
                      (n_onesamplettestdes, n_level2estimate, [('spm_mat_file', 'spm_mat_file')]),
                      (n_level2estimate, n_level2conestimate, [('spm_mat_file', 'spm_mat_file'),
                                                               ('beta_images', 'beta_images'),
                                                               ('residual_image', 'residual_image')]),
                      (n_level2conestimate, n_level2thresh, [('spm_mat_file', 'spm_mat_file'),
                                                             ('spmT_images', 'stat_image')]),
                      (n_level2conestimate, n_datasink, [('spm_mat_file', results_dir + '.@spm_mat'),
                                                         ('spmT_images', results_dir + '.@T'),
                                                         ('con_images', results_dir + '.@con')]),
                      (n_level2thresh, n_datasink, [('thresholded_map', results_dir + '.@threshold')])
                      ])

print('\nRunning analysis output spaces: %s, contrasts: %s \n' % (str(output_space), str(contrast_ids)))

w_l2analysis.write_graph(graph2use='flat', format='pdf')
w_l2analysis.run('MultiProc', plugin_args={'n_procs': nr_procs})

# os.system('cp %s %s' % (opj(experiment_dir, working_dir, w_l2analysis.name, 'graph_detailed.png'),
#                         opj(experiment_dir, output_dir, results_dir)))
# os.system('cp %s %s' % (
# opj(experiment_dir, working_dir, w_l2analysis.name, 'graph.png'), opj(experiment_dir, output_dir, results_dir)))
#os.system('rm -rf %s' % opj(experiment_dir, working_dir))  # del working dir after workflow is finished
