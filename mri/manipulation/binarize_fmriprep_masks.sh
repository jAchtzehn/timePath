#!/bin/bash


# script that binarizes the anatomical grey matter masks provided by fmriprep and saves them into 
# the nilearn data folder, for both T1w and MNI space

participantNr_start=$1	 	
participantNr_end=$2

data_dir=/Users/jAchtzehn/data/fMRI/timePath
#data_dir=/media/sf_data/fMRI/timePath

for id in $(seq -f "%02g" $participantNr_start $participantNr_end)
do
	
	fmriprep_dir=$data_dir/fmriprep/sub-$id/anat/
	nilearn_dir=$data_dir/nilearn/
	mkdir -p $result_dir
	
	# T1w space
	# binarize label files from freesurfer
	mri_binarize --i $fmriprep_dir/'sub-'$id'_T1w_class-GM_probtissue.nii.gz' --o $nilearn_dir/sub-$id/space-T1w/masks/'sub-'$id'_gm_mask_tmp.nii.gz' --min 0.9 \
		--mask $fmriprep_dir/'sub-'$id'_T1w_brainmask.nii.gz'
	#mri_binarize --i $fmriprep_dir/'sub-'$id'_T1w_space-MNI152NLin2009cAsym_class-GM_probtissue.nii.gz' --o $nilearn_dir/sub-$id/space-MNI152NLin2009cAsym/masks/'sub-'$id'_gm_mask_binarized.nii.gz' --min 0.9 \
	#	--mask $fmriprep_dir/'sub-'$id'_T1w_space-MNI152NLin2009cAsym_brainmask.nii.gz'
	
	# clusterize to get rid of extra voxels
	3dclusterize -inset $nilearn_dir/sub-$id/space-T1w/masks/'sub-'$id'_gm_mask_tmp.nii.gz' -binary -NN 3 -ithr 0 -clust_nvox 100 \
		-pref_map $nilearn_dir/sub-$id/space-T1w/masks/'sub-'$id'_gm_mask_binarized_cluster_tmp.nii.gz'\
			-1sided RIGHT_TAIL 0.5
	
	# resample binary file with AFNI 3dresample to match 3x3x3.3mm resolution of functional images for both hemispheres
	3dresample -master $data_dir/fmriprep/sub-$id/ses-01/func/sub-${id}_ses-01_task-class_run-01_bold_space-T1w_preproc.nii.gz \
		-prefix $nilearn_dir/sub-$id/space-T1w/masks/'sub-'$id'_gm_mask_binarized.nii.gz' \
			-inset $nilearn_dir/sub-$id/space-T1w/masks/'sub-'$id'_gm_mask_binarized_cluster_tmp.nii.gz'
	
	# resample binary file with AFNI 3dresample to match 3x3x3.3mm resolution of functional images for rh
	#3dresample -master $data_dir/fmriprep/sub-$id/ses-01/func/sub-${id}_ses-01_task-class_run-01_bold_space-T1w_preproc.nii.gz -prefix $result_dir/'sub-'$id'_hc_rh_mask_binarized.nii.gz' -inset 	#$data_dir/freesurfer/sub-$id/mri/rh_binarized.nii
	
	# resample binary file with AFNI 3dresample to match 3x3x3.3mm resolution of functional images for lh
	#3dresample -master $data_dir/fmriprep/sub-$id/ses-01/func/sub-${id}_ses-01_task-class_run-01_bold_space-T1w_preproc.nii.gz -prefix $result_dir/'sub-'$id'_hc_lh_mask_binarized.nii.gz' -inset 	#$data_dir/freesurfer/sub-$id/mri/lh_binarized.nii
	
	
done
