#!/bin/bash

participantNr_start=$1	 	
participantNr_end=$2

data_dir=/Users/jAchtzehn/data/fMRI/timePath
#data_dir=/media/sf_data/fMRI/timePath

for id in $(seq -f "%02g" $participantNr_start $participantNr_end)
do
	mkdir -p $data_dir/nilearn/sub-$id/space-MNI152NLin2009cAsym/masks
	
	for mask in hc ifg ips pcg insula ba6
	do
		
		antsApplyTransforms -d 3 -e 0 -v 1\
			-i $data_dir/nilearn/sub-$id/space-T1w/masks/sub-$id'_'$mask'_mask_binarized.nii.gz' \
				-r $data_dir/fmriprep/sub-$id/ses-01/func/sub-$id'_ses-01_task-class_run-01_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz' \
					-t $data_dir/fmriprep/sub-$id/anat/sub-$id'_T1w_target-MNI152NLin2009cAsym_warp.h5' \
						-o $data_dir/nilearn/sub-$id/space-MNI152NLin2009cAsym/masks/sub-$id'_'$mask'_mask.nii.gz'

		mri_binarize --i $data_dir/nilearn/sub-$id/space-MNI152NLin2009cAsym/masks/sub-$id'_'$mask'_mask.nii.gz'\
			--o $data_dir/nilearn/sub-$id/space-MNI152NLin2009cAsym/masks/sub-$id'_'$mask'_mask_binarized.nii.gz' --min 0.0001
		
		rm $data_dir/nilearn/sub-$id/space-MNI152NLin2009cAsym/masks/sub-$id'_'$mask'_mask.nii.gz'
		
		for hs in lh rh
		do
			antsApplyTransforms -d 3 -e 0 -v 1\
				-i $data_dir/nilearn/sub-$id/space-T1w/masks/sub-$id'_'$mask'_'$hs'_mask_binarized.nii.gz' \
					-r $data_dir/fmriprep/sub-$id/ses-01/func/sub-$id'_ses-01_task-class_run-01_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz' \
						-t $data_dir/fmriprep/sub-$id/anat/sub-$id'_T1w_target-MNI152NLin2009cAsym_warp.h5' \
							-o $data_dir/nilearn/sub-$id/space-MNI152NLin2009cAsym/masks/sub-$id'_'$mask'_'$hs'_mask.nii.gz'
	
			mri_binarize --i $data_dir/nilearn/sub-$id/space-MNI152NLin2009cAsym/masks/sub-$id'_'$mask'_'$hs'_mask.nii.gz'\
				--o $data_dir/nilearn/sub-$id/space-MNI152NLin2009cAsym/masks/sub-$id'_'$mask'_'$hs'_mask_binarized.nii.gz' --min 0.0001
			
			rm $data_dir/nilearn/sub-$id/space-MNI152NLin2009cAsym/masks/sub-$id'_'$mask'_'$hs'_mask.nii.gz'
		done
	done
	
done