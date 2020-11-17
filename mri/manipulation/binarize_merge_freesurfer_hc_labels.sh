#!/bin/bash

participantNr_start=$1	 	
participantNr_end=$2

data_dir=/Users/jAchtzehn/data/fMRI/timePath
#data_dir=/media/sf_data/fMRI/timePath



for id in $(seq -f "%02g" $participantNr_start $participantNr_end)
do
	
	result_dir=$data_dir/nilearn/sub-$id/masks
	mkdir -p $result_dir
	
	# binarize label files from freesurfer
	mri_binarize --i $data_dir/freesurfer/sub-$id/mri/lh.hippoSfLabels-T1.v10.FSvoxelSpace.mgz --o $data_dir/freesurfer/sub-$id/mri/lh_binarized.mgz --min 0.0001
	mri_binarize --i $data_dir/freesurfer/sub-$id/mri/rh.hippoSfLabels-T1.v10.FSvoxelSpace.mgz --o $data_dir/freesurfer/sub-$id/mri/rh_binarized.mgz --min 0.0001
	mri_binarize --i $data_dir/freesurfer/sub-$id/mri/rh.hippoSfLabels-T1.v10.FSvoxelSpace.mgz --o $data_dir/freesurfer/sub-$id/mri/hc_binarized.mgz --min 0.0001 --merge $data_dir/freesurfer/sub-$id/mri/lh_binarized.mgz	# binarize and merge

	mri_convert $data_dir/freesurfer/sub-$id/mri/hc_binarized.mgz $data_dir/freesurfer/sub-$id/mri/hc_binarized.nii.gz		# convert to nifti format for both hemispheres
	mri_convert $data_dir/freesurfer/sub-$id/mri/rh_binarized.mgz $data_dir/freesurfer/sub-$id/mri/rh_binarized.nii.gz		# convert to nifti format for rh 
	mri_convert $data_dir/freesurfer/sub-$id/mri/lh_binarized.mgz $data_dir/freesurfer/sub-$id/mri/lh_binarized.nii.gz		# convert to nifti format for lh
	

	# resample binary file with AFNI 3dresample to match 3x3x3.3mm resolution of functional images for both hemispheres
	3dresample -master $data_dir/fmriprep/sub-$id/ses-01/func/sub-${id}_ses-01_task-class_run-01_bold_space-T1w_preproc.nii.gz -prefix $result_dir/'sub-'$id'_hc_mask_binarized.nii.gz' -inset 	$data_dir/freesurfer/sub-$id/mri/hc_binarized.nii
	
	# resample binary file with AFNI 3dresample to match 3x3x3.3mm resolution of functional images for rh
	3dresample -master $data_dir/fmriprep/sub-$id/ses-01/func/sub-${id}_ses-01_task-class_run-01_bold_space-T1w_preproc.nii.gz -prefix $result_dir/'sub-'$id'_hc_rh_mask_binarized.nii.gz' -inset 	$data_dir/freesurfer/sub-$id/mri/rh_binarized.nii
	
	# resample binary file with AFNI 3dresample to match 3x3x3.3mm resolution of functional images for lh
	3dresample -master $data_dir/fmriprep/sub-$id/ses-01/func/sub-${id}_ses-01_task-class_run-01_bold_space-T1w_preproc.nii.gz -prefix $result_dir/'sub-'$id'_hc_lh_mask_binarized.nii.gz' -inset 	$data_dir/freesurfer/sub-$id/mri/lh_binarized.nii
	
	
done
