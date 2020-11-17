from read_logs import get_classification_onsets as getonsets
from util_functions import createSubjectlist, cleanWorkingOutputDir, plotDesignMatrix, plotTrialInfo, getVolNr
import os
from os.path import join as opj
from os.path import abspath
from sys import platform
import datetime

from nipype import config, logging
from nipype.pipeline.engine import Workflow, Node
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.algorithms.misc import Gunzip
from nipype.algorithms.modelgen import SpecifySPMModel
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.interfaces.fsl import ExtractROI

import nipype.interfaces.spm as spm

# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
	behav_dir = abspath(
		'/Users/jachtzehn/data/fMRI/timePath/behavioural')  # where is the onset and duration information stored?
else:
	experiment_dir = abspath('/mnt/work/achtzehnj/data')
	behav_dir = abspath('/mnt/work/achtzehnj/data/behavioural')  # where is the onset and duration information stored?

working_dir = 'working_nipype'  # working directory for nipype
output_dir = 'output_nipype'  # where should the relevant output be stored?
results_dir = 'tw_beta_estimation'
workflow_name = 'w_l1analysis'  # name of the workflow
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
	print('Output dir %s does not exist, creating it...' % opj(experiment_dir, output_dir))

# ------------ Options ------------
nr_procs = 24  # how many processors should be utilized?
subjects = range(1, 26)  # subjects
ses_ids = range(1, 3)  # session ids to be analysed
run_ids = range(1, 5)  # run iterations, as named in the files
trial_ids = range(1, 66)  # each run consists of 65 trials (including the cross (null) trials)
TR = 2  # TR in seconds
hpf_cutoff = 128  # high-pass filter cutoff in Hz
additionalRegressors = True  # should FD and aCompCor regressors be included?
output_space = ['MNI152NLin2009cAsym']  # T1w and/or MNI152NLin2009cAsym are possible

# nipype config
config_dict = {'execution': {'crashdump_dir': log_dir,
                             'poll_sleep_duration': 10,
                             'crashfile_format': 'txt',
                             'stop_on_first_crash': 'true'},
               'logging': {'util_level': 'INFO',
                           'log_directory': log_dir,
                           'log_to_file': 'true',
                           'log_size': 95400000},
               'monitoring': {'enabled': 'true',
                              'summary_file': log_dir}}

config.update_config(config_dict)
logging.update_logging(config)

# ------------ conditions ------------
condition_names = ['trial', 'time', 'dist', 'lumin', 'dots', 'cross', 'icon', 'comp', 'warning']  # define conditions that should be regarded in the GLM, the first condition is the one for the current trial

# ------------ Preparation -------------
subj_list = createSubjectlist(subjects)  # create a list with subjects in the format of ['sub-01', 'sub-02', ...]
# ------------ Node definition ------------

print('Running analysis for subjects: %s, sessions: %s, runs: %s, trials: %s - %s, output spaces: %s' %
      (str(subj_list), str(ses_ids), str(run_ids), str(trial_ids[0]), str(trial_ids[len(trial_ids) - 1]),
       str(output_space)))

# --- own functions
n_getOnsets = Node(
	Function(input_names=['subj_id', 'base_dir', 'ses_id', 'run_id', 'trial_id', 'addRegressors', 'conditions'],
	         output_names=['subj_info', 'trial_info', 'regr_info'], function=getonsets),
	name='n_getOnsets')

n_getOnsets.inputs.base_dir = experiment_dir
n_getOnsets.inputs.addRegressors = additionalRegressors
n_getOnsets.inputs.conditions = condition_names

n_cleanWorkingOutputDir = Node(
	Function(input_names=['output_space', 'subj_id', 'ses_id', 'run_id', 'trial_id', 'working_dir', 'datasink_out'],
	         function=cleanWorkingOutputDir), name='n_cleanWorkinOutputDir')
n_cleanWorkingOutputDir.inputs.working_dir = opj(experiment_dir, working_dir, workflow_name)

n_plotDesignMatrix = Node(
	Function(input_names=['matFile'], output_names=['fig_filename_full'], function=plotDesignMatrix),
	name='n_plotDesignMatrix')
n_plotTrialInfo = Node(Function(input_names=['trial_info', 'regr_info'], output_names=['fname_trial', 'fname_regr'],
                                function=plotTrialInfo), name='n_plotTrialInfo')

n_getVolNr = Node(
	Function(input_names=['behav_dir', 'subj_id', 'ses_id', 'run_id'], output_names=['nrVol'], function=getVolNr),
	name='n_getVolNr')
n_getVolNr.inputs.behav_dir = behav_dir

# --- utility
n_gunzip_mask = Node(Gunzip(), name='n_gunzip_mask')

n_extractroi = Node(ExtractROI(), name='n_extractroi')
n_extractroi.inputs.output_type = 'NIFTI'
n_extractroi.inputs.t_min = 0

# --- model
# generate SPM-specific model
n_modelspec = Node(SpecifySPMModel(concatenate_runs=False, input_units='secs', output_units='secs', time_repetition=TR,
                                   high_pass_filter_cutoff=hpf_cutoff),
                   name='n_modelspec')
n_modelspec.inputs.parameter_source = 'FSL'  # fmriprep uses FSL's mcflirt algorithm for head motion estimation

# level1Design - Generates an SPM design matrix
n_level1design = Node(spm.Level1Design(bases={'hrf': {'derivs': [0, 0]}}, timing_units='secs', interscan_interval=TR,
                                       model_serial_correlations='FAST',
                                       microtime_resolution=36, microtime_onset=1),
                      name='n_level1design')

# estimate model
n_level1estimate = Node(spm.EstimateModel(estimation_method={'Classical': 1}),
                        name='n_level1estimate')

# --- IO stream

# infosource
n_infosource = Node(IdentityInterface(fields=['subj_id', 'output_space', 'ses_id', 'run_id', 'trial_id']),
                    name="n_infosource")
n_infosource.iterables = [('subj_id', subj_list),
                          ('output_space', output_space),
                          ('ses_id', ses_ids),
                          ('run_id', run_ids),
                          ('trial_id', trial_ids)]

# select files
if output_space[0] == 'MNI152NLin2009cAsym':
	mask_file_name = '{subj_id}_T1w_space-{output_space}_brainmask.nii.gz'
else:
	mask_file_name = '{subj_id}_T1w_brainmask.nii.gz'

templates = {'func': opj('fmriprep', '{subj_id}', 'ses-0{ses_id}', 'func',
                         '{subj_id}_ses-0{ses_id}_task-class_run-0{run_id}_bold_space-{output_space}_preproc.nii.gz'),
             'mc_param': opj('fmriprep', '{subj_id}', 'ses-0{ses_id}', 'func',
                             '{subj_id}_ses-0{ses_id}_task-class_run-0{run_id}_bold_confounds_mc_params_cut.par'),
             'mask_file': opj('fmriprep', '{subj_id}', 'anat', mask_file_name)}

n_selectfiles = Node(SelectFiles(templates,
                                 base_directory=experiment_dir,
                                 sort_filelist=True),
                     name="n_selectfiles")

# datasink
n_datasink = Node(DataSink(base_directory=experiment_dir, container=output_dir),
                  name="n_datasink")

substitutions = []
subjFolders = [('_output_space_%s_run_id_%s_ses_id_%s_subj_id_%s_trial_id_%s' % (o, run, ses, sub, trial),
                '%s/space-%s/ses-%s/run-%s/trial-%s' % (sub, o, ses, run, trial))
               for o in output_space
               for run in run_ids
               for ses in ses_ids
               for sub in subj_list
               for trial in trial_ids]
substitutions.extend(subjFolders)
n_datasink.inputs.substitutions = substitutions

# ------------ Connect workflows ------------
w_l1analysis = Workflow(base_dir=opj(experiment_dir, working_dir), name=workflow_name)

w_l1analysis.connect([(n_infosource, n_selectfiles, [('subj_id', 'subj_id'),
                                                     ('output_space', 'output_space'),
                                                     ('ses_id', 'ses_id'),
                                                     ('run_id', 'run_id')]),
                      (n_infosource, n_getOnsets, [('subj_id', 'subj_id'),
                                                   ('ses_id', 'ses_id'),
                                                   ('run_id', 'run_id'),
                                                   ('trial_id', 'trial_id')]),
                      (n_getOnsets, n_plotTrialInfo, [('trial_info', 'trial_info'),
                                                      ('regr_info', 'regr_info')]),
                      (n_infosource, n_getVolNr, [('subj_id', 'subj_id'),
                                                  ('ses_id', 'ses_id'),
                                                  ('run_id', 'run_id')]),
                      (n_getVolNr, n_extractroi, [('nrVol', 't_size')]),
                      (n_selectfiles, n_extractroi, [('func', 'in_file')]),
                      (n_selectfiles, n_gunzip_mask, [('mask_file', 'in_file')]),
                      (n_getOnsets, n_modelspec, [('subj_info', 'subject_info')]),
                      (n_extractroi, n_modelspec, [('roi_file', 'functional_runs')]),
                      (n_selectfiles, n_modelspec, [('mc_param', 'realignment_parameters')]),
                      (n_gunzip_mask, n_level1design, [('out_file', 'mask_image')]),
                      (n_modelspec, n_level1design, [('session_info', 'session_info')]),
                      (n_level1design, n_level1estimate, [('spm_mat_file', 'spm_mat_file')]),
                      (n_level1estimate, n_plotDesignMatrix, [('spm_mat_file', 'matFile')]),
                      (n_plotDesignMatrix, n_datasink, [('fig_filename_full', results_dir + '.@spm_mat_full')]),
                      (n_plotTrialInfo, n_datasink, [('fname_trial', results_dir + '.@trial_info'),
                                                     ('fname_regr', results_dir + '.@regr_info')]),
                      (n_level1estimate, n_datasink, [('beta_images', results_dir + '.@beta_images'),
                                                      # ('spm_mat_file', results_dir + '.@spm_mat_file'),
                                                      ('RPVimage', results_dir + '.@RPVimage'),
                                                      ('residual_image', results_dir + '.@residual_image'),
                                                      ('mask_image', results_dir + '.@mask_image')]),
                      (n_datasink, n_cleanWorkingOutputDir, [('out_file', 'datasink_out')]),
                      (n_infosource, n_cleanWorkingOutputDir, [('output_space', 'output_space'),
                                                               ('subj_id', 'subj_id'),
                                                               ('ses_id', 'ses_id'),
                                                               ('run_id', 'run_id'),
                                                               ('trial_id', 'trial_id')])
                      ])

# -- run workflow and cp graphs
w_l1analysis.write_graph(graph2use='flat', simple_form=True, format='pdf')
w_l1analysis.write_graph(graph2use='flat', format='pdf')
#w_l1analysis.run('MultiProc', plugin_args={'n_procs': nr_procs})

os.system('cp %s %s' % (opj(experiment_dir, working_dir, w_l1analysis.name, 'graph_detailed.png'),
                        opj(experiment_dir, output_dir, results_dir)))
os.system('cp %s %s' % (
	opj(experiment_dir, working_dir, w_l1analysis.name, 'graph.png'), opj(experiment_dir, output_dir, results_dir)))

# cleanup
#os.system('rm -rf %s' % opj(experiment_dir, working_dir))  # del working dir after workflow is finished
