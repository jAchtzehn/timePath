#!/bin/bash

# resample whole brain masks to 3x3x3mm resolution of functional runs
participantNr_start=$1	 	
participantNr_end=$2

data_dir=/Users/jAchtzehn/data/fMRI/timePath
#data_dir=/media/sf_data/fMRI/timePath



for id in $(seq -f "%02g" $participantNr_start $participantNr_end)
do
	
	result_dir=$data_dir/nilearn/sub-$id/space-MNI152NLin2009cAsym/masks
	mkdir -p $result_dir
	
	# resample binary file with AFNI 3dresample to match 3x3x3.3mm resolution of functional images for both hemispheres
	3dresample -master $data_dir/fmriprep/sub-$id/ses-01/func/sub-${id}_ses-01_task-class_run-01_bold_space-MNI152NLin2009cAsym_preproc.nii.gz -prefix $result_dir/'sub-'$id'_wb_mask_binarized.nii.gz' -inset \
		$data_dir/fmriprep/sub-$id/anat/'sub-'$id'_T1w_space-MNI152NLin2009cAsym_brainmask.nii.gz'
	
done
