from read_logs import get_classification_onsets as getonsets
from util_functions import createSubjectlist, getVolNr
import os
from os.path import join as opj
from os.path import abspath
from sys import platform
import datetime

from nipype import config, logging
from nipype.utils.profiler import log_nodes_cb, get_system_total_memory_gb
from nipype.utils.draw_gantt_chart import generate_gantt_chart

from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.algorithms.misc import Gunzip
from nipype.algorithms.modelgen import SpecifySPMModel
from nipype.interfaces.io import SelectFiles, DataSink
import nipype.interfaces.spm as spm
from nipype.interfaces.fsl import ExtractROI

# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
	behav_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/behavioural')  # where is the onset and duration information stored?
else:
	experiment_dir = abspath('/mnt/work/achtzehnj/data')
	behav_dir = abspath('/mnt/work/achtzehnj/data/behavioural')  # where is the onset and duration information stored?

working_dir = 'working_nipype'                                          # working directory for nipype
output_dir = 'output_nipype'                                            # where should the relevant output be stored?
results_dir = '1stlevel_tmp'

print('Deleting work folder %s...' % opj(experiment_dir, working_dir))
os.system('rm -rf %s'%opj(experiment_dir, working_dir))                 # del working dir for clean start

log_dir_timestamp = datetime.datetime.now().strftime("%y-%m-%d-%H-%M-%S")
log_dir = opj(experiment_dir, 'logs_nipype', results_dir + '_' + log_dir_timestamp)

# create log folders for logging if they have been deleted
os.makedirs(log_dir)    # make output dir for current workflow

# check if output_dir already exists and if so if it should be deleted
try:
	if results_dir in os.listdir(opj(experiment_dir, output_dir)):
		print('Output dir %s already exists'%opj(experiment_dir, output_dir, results_dir))
		rmOutput = input('Should it be deleted before proceeding? (y/n): ')
		answerGiven = False
	
		while not answerGiven:
			if rmOutput == 'y':
				os.system('rm -rf %s'%opj(experiment_dir, output_dir, results_dir))
				print('Deleted folder %s'%opj(experiment_dir, output_dir, results_dir))
				answerGiven = True
			elif rmOutput == 'n':
				print('Not deleting output folder %s'%opj(experiment_dir, output_dir, results_dir))
				answerGiven = True
			else:
				print('Please enter (y/n)')
				answerGiven = False
except OSError:
	print('Output dir %s does not exist, proceeding...'%opj(experiment_dir, output_dir))

# ------------ Options ------------
settings = {'n_proc': 15,
            'perc_gb': 0.9,
            'dummyTRs': 0,
            'subjects': range(1, 2),
            'fwhm_size': 6,
            'TR': 2,
            'hpf_cutoff': 128,
            'additionalRegressors': True,
            'output_space': ['MNI152NLin2009cAsym'],
            'stop_on_first_crash': 'False',
            'compress_mat_files': False}

# ------------ nipype config ------------
config_dict = {'execution': {'crashdump_dir': log_dir,
                             'poll_sleep_duration': 2,
                             'crashfile_format': 'txt',
                             'profile_runtime': 'true'} ,
               'logging': {'utils_level': 'INFO',
                         'log_directory':  log_dir,
                         'log_to_file': 'true',
                         'log_size': 95400000},
               'monitoring': {'enabled': 'true',
                              'summary_file': log_dir}}
# file logging
config.update_config(config_dict)
logging.update_logging(config)
config.enable_resource_monitor()

# ------------ contrasts ------------
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
# ------------ Preparation -------------
subj_list = createSubjectlist(settings['subjects'])                         # create a list with subjects in the format of ['sub-01', 'sub-02', ...]

# ------------ Node definition ------------

print('\n\n -- Review of settings -- \n')
for key in settings.keys():
	print('Setting {}: \t{}'.format(key, settings[key]))

answerGiven = False
rmOutput = input('Proceed? (y/n): ')

while not answerGiven:
	if rmOutput == 'y':
		answerGiven = True
	elif rmOutput == 'n':
		exit()
	else:
		print('Please enter (y/n)')
		answerGiven = False

# --- own functions
n_getOnsets = Node(Function(input_names=['subj_id', 'base_dir', 'addRegressors', 'conditions'], output_names=['subj_info'], function=getonsets),
					  name='n_getOnsets')
n_getOnsets.inputs.base_dir = experiment_dir
n_getOnsets.inputs.addRegressors = settings['additionalRegressors']
n_getOnsets.inputs.conditions = condition_names

n_getVolNr = Node(Function(input_names=['behav_dir', 'subj_id'], output_names=['nrVol'], function=getVolNr),
                  name='n_getVolNr')
n_getVolNr.inputs.behav_dir = behav_dir

# --- utility
n_gunzip_func = MapNode(Gunzip(), name="n_gunzip_func", iterfield=['in_file'])
n_gunzip_mask = Node(Gunzip(), name="n_gunzip_mask")

n_extractroi = MapNode(ExtractROI(), name='n_extractroi', iterfield=['in_file', 't_size'])
n_extractroi.inputs.output_type = 'NIFTI'
n_extractroi.inputs.t_min = 0

# preproc
n_smooth = Node(spm.Smooth(fwhm=settings['fwhm_size']),
			  name='n_smooth')
n_smooth.interface.mem_gb = 0.5

# --- model
# generate SPM-specific model
n_modelspec = Node(SpecifySPMModel(concatenate_runs=False, input_units='secs', output_units='secs', time_repetition=settings['TR'], high_pass_filter_cutoff=settings['hpf_cutoff']),
				 name='n_modelspec')
n_modelspec.inputs.parameter_source = 'FSL'         # fmriprep uses FSL's mcflirt algorithm for head motion estimation

# level1Design - Generates an SPM design matrix
n_level1design = Node(spm.Level1Design(bases={'hrf': {'derivs': [0, 0]}}, timing_units='secs', interscan_interval=settings['TR'], model_serial_correlations='FAST',
                                       microtime_resolution=36, microtime_onset=1),
					name='n_level1design')
n_level1design.interface.mem_gb = 2

# estimate model
n_level1estimate = Node(spm.EstimateModel(estimation_method={'Classical': 1}),
						   name='n_level1estimate')
n_level1estimate.interface.mem_gb = 6
n_level1estimate.interface.n_procs = 8

# estimate contrast
n_level1conest = Node(spm.EstimateContrast(),
						 name='n_level1conest')

# --- IO stream
# infosource
n_infosource = Node(IdentityInterface(fields=['subj_id', 'output_space', 'contrasts'], contrasts=contrast_list),
					   name="n_infosource")
n_infosource.iterables = [('subj_id', subj_list),
						  ('output_space', settings['output_space'])]

if settings['output_space'] == 'T1w':
	mask_image_name = '{subj_id}_T1w_brainmask.nii.gz'
else:
	mask_image_name = '{subj_id}_T1w_space-{output_space}_brainmask.nii.gz'
	
templates = {'func': opj('fmriprep', '{subj_id}', 'ses-0*', 'func', '{subj_id}_ses-0*_task-class_run-0*_bold_space-{output_space}_preproc.nii.gz'),
			 'mc_param': opj('fmriprep', '{subj_id}', 'ses-0*', 'func', '{subj_id}_ses-0*_task-class_run-0*_bold_confounds_mc_params_cut.par'),
             'mask_file': opj('fmriprep', '{subj_id}', 'anat', mask_image_name)}

n_selectfiles = Node(SelectFiles(templates,
							   base_directory=experiment_dir,
							   sort_filelist=True),
				   name="n_selectfiles")

# datasink
n_datasink = Node(DataSink(base_directory=experiment_dir, container=output_dir),
					   name="n_datasink")

substitutions = []
subjFolders = [('_output_space_%s_subj_id_%s' % (o, sub), '%s/space-%s' % (sub, o))
			   for o in settings['output_space']
			   for sub in subj_list]
substitutions.extend(subjFolders)
n_datasink.inputs.substitutions = substitutions

# ------------ Connect workflows ------------
w_l1analysis = Workflow(base_dir=opj(experiment_dir, working_dir), name='w_l1analysis')
w_l1analysis.config['execution'] = {'stop_on_first_crash': settings['stop_on_first_crash']}

w_l1analysis.connect([(n_infosource, n_selectfiles, [('subj_id', 'subj_id'),
													 ('output_space', 'output_space')]),
					  (n_infosource, n_getOnsets, [('subj_id', 'subj_id')]),
                      (n_infosource, n_getVolNr, [('subj_id', 'subj_id')]),
                      (n_getVolNr, n_extractroi, [('nrVol', 't_size')]),
					  (n_getOnsets, n_modelspec, [('subj_info', 'subject_info')]),
					  (n_infosource, n_level1conest, [('contrasts', 'contrasts')]),
					  (n_selectfiles, n_extractroi, [('func', 'in_file')]),
					  (n_selectfiles, n_gunzip_mask, [('mask_file', 'in_file')]),
					  (n_extractroi, n_smooth, [('roi_file', 'in_files')]),
					  (n_smooth, n_modelspec, [('smoothed_files', 'functional_runs')]),
					  (n_selectfiles, n_modelspec, [('mc_param', 'realignment_parameters')]),
					  (n_modelspec, n_level1design, [('session_info', 'session_info')]),
                      (n_gunzip_mask, n_level1design, [('out_file', 'mask_image')]),
					  (n_level1design, n_level1estimate, [('spm_mat_file', 'spm_mat_file')]),
					  (n_level1estimate, n_level1conest, [('spm_mat_file', 'spm_mat_file'),
														  ('beta_images', 'beta_images'),
														  ('residual_image','residual_image')]),
                      (n_level1estimate, n_datasink, [('RPVimage', results_dir + '.@RPVimage'),
					                                  ('residual_image', results_dir + '.@residual_image'),
					                                  ('mask_image', results_dir + '.@mask_image')]),
					  (n_level1conest, n_datasink, [('spm_mat_file', results_dir + '.@spm_mat'),
													('spmT_images', results_dir + '.@T'),
													('con_images', results_dir + '.@con'),
													('spmF_images', results_dir + '.@F'),
													('ess_images', results_dir + '.@ess'),
													])
					  ])

# -- run workflow and cp graphs
w_l1analysis.write_graph(graph2use='flat', format='pdf', simple_form=False)
#w_l1analysis.run('MultiProc', plugin_args={'n_procs': settings['n_proc'],
#                                           'memory_gb': settings['perc_gb'] * get_system_total_memory_gb()})

#os.system('cp %s %s'% (opj(experiment_dir, working_dir, w_l1analysis.name, 'graph_detailed.png'), opj(experiment_dir, output_dir, results_dir)))
#os.system('cp %s %s'% (opj(experiment_dir, working_dir, w_l1analysis.name, 'graph.png'), opj(experiment_dir, output_dir, results_dir)))
#os.system('rm -rf %s'%opj(experiment_dir, working_dir))                 # del working dir after workflow is finished

if settings['compress_mat_files']:
	print('Compressing SPM.mat files...')
	# compress .mat file to save space
	for subj in settings['subjects']:
		subj_str = 'sub-' + str(subj).zfill(2)
		
		for space in settings['output_space']:
			os.system('gzip %s' % opj(experiment_dir, output_dir, results_dir, subj_str, 'space-' + space, 'SPM.mat'))

print('done.')