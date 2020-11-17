import os
from os.path import join as opj
from os.path import abspath
from sys import platform
import datetime
from nipype import config, logging
from nipype.interfaces.utility import Function
from nipype.pipeline.engine import Workflow, Node
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.interfaces.spm import Threshold
from util_functions import conj_stat_maps, ClusterThresh


# ------------ File I/O ------------
if platform == 'darwin':
	experiment_dir = abspath('/Users/jachtzehn/data/fMRI/timePath/')
	behav_dir = abspath(
		'/Users/jachtzehn/data/fMRI/timePath/behavioural')  # where is the onset and duration information stored?
else:
	experiment_dir = abspath('/media/sf_data/fMRI/timePath/')
	behav_dir = abspath(
		'/media/sf_data/fMRI/timePath/behavioural')  # where is the onset and duration information stored?

working_dir = 'working_nipype'  # working directory for nipype
output_dir = 'output_nipype'  # where should the relevant output be stored?
second_level_dir = '2ndlevel'  # parent dir of 1st level results
results_dir = '2ndlevel_conj'  # output of second level results
l2_subdir = opj(experiment_dir, output_dir, second_level_dir)

log_dir_timestamp = datetime.datetime.now().strftime("%y-%m-%d-%H-%M-%S")
log_dir = opj(experiment_dir, 'logs_nipype', results_dir + '_' + log_dir_timestamp)
os.makedirs(log_dir)    # make output dir for current workflow

os.system('rm -rf %s' % opj(experiment_dir, working_dir))                 # del working dir before start
# ------------ Options ------------
space = 'MNI152NLin2009cAsym'
con_list = ['0011', '0015']
con1 = con_list[0]
con2 = con_list[1]
stop_on_first_crash = True
extent_threshold = 0
nr_procs = 4
del_working_dir = True

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


# selectfiles to select maps
n_selectfiles = Node(SelectFiles({'map1': opj(l2_subdir, 'space-' + space, 'contrast-' + con1, 'spmT_0001_thr.nii'),
                                  'map2': opj(l2_subdir, 'space-' + space, 'contrast-' + con2, 'spmT_0001_thr.nii')},
                                 base_directory=experiment_dir,
                                 sort_filelist=True),
                     name="selectfiles")
    
n_conjmaps = Node(Function(input_names=['map1', 'map2'],
                           output_names=['out_file'], function=conj_stat_maps),
                  name='conjmaps')
    
# Threshold - thresholds contrasts
n_level2thresh = Node(ClusterThresh(matlab_path='/Applications/MATLAB_R2016b.app/bin', extent_threshold=60),
                      name='n_level2thresh')

n_datasink = Node(DataSink(base_directory=experiment_dir, container=output_dir),
                  name="n_datasink")

w_l2analysis = Workflow(base_dir=opj(experiment_dir, working_dir), name='w_l2analysis')
w_l2analysis.config['execution'] = {'stop_on_first_crash': stop_on_first_crash}

# Connect up the 2nd-level analysis components
w_l2analysis.connect([(n_selectfiles, n_conjmaps, [('map1', 'map1'), ('map2', 'map2')]),
                      (n_conjmaps, n_level2thresh, [('out_file', 'in_file')]),
                      (n_level2thresh, n_datasink, [('out_file', results_dir + '.@conj_thr')]),
                      (n_conjmaps, n_datasink, [('out_file', results_dir + '.@conj')])])

w_l2analysis.run('MultiProc', plugin_args={'n_procs': nr_procs})

if del_working_dir:
	os.system('rm -rf %s' % opj(experiment_dir, working_dir))                 # del working dir after workflow is finished
