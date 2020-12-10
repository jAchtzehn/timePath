#!/bin/bash

# resample and convert anatomical ROIs create by the anatomy toolbox for SPM

result_dir=/home/achtzehnj/data/timePath/derivatives/nilearn/group_masks/space-MNI152NLin2009cAsym

# resample binary file with AFNI 3dresample to match 3x3x3.3mm resolution of functional images for both hemispheres
#3dresample -master $result_dir/group_mask_V5_binarized_old.nii -prefix $result_dir/group_mask_hOc4la_lh_binarized.nii.gz -inset \
#	$result_dir/ROI_Visual_hOc4la_L_MNI.nii

#3dresample -master $result_dir/group_mask_V5_binarized_old.nii -prefix $result_dir/group_mask_hOc4la_rh_binarized.nii.gz -inset \
#	$result_dir/ROI_Visual_hOc4la_R_MNI.nii

#3dresample -master $result_dir/group_mask_V5_binarized_old.nii -prefix $result_dir/group_mask_hOc5_lh_binarized.nii.gz -inset \
#	$result_dir/ROI_Visual_hOc5_L_MNI.nii

#3dresample -master $result_dir/group_mask_V5_binarized_old.nii -prefix $result_dir/group_mask_hOc5_rh_binarized.nii.gz -inset \
#	$result_dir/ROI_Visual_hOc5_R_MNI.nii


# combine hOc4la and hOc5 regions to from V5

#fslmaths $result_dir/group_mask_hOc4la_lh_binarized.nii.gz \
#	-add $result_dir/group_mask_hOc5_lh_binarized.nii.gz \
#		$result_dir/group_mask_V5_lh_binarized.nii.gz

#fslmaths $result_dir/group_mask_hOc4la_rh_binarized.nii.gz \
#	-add $result_dir/group_mask_hOc5_rh_binarized.nii.gz \
#		$result_dir/group_mask_V5_rh_binarized.nii.gz

#fslmaths $result_dir/group_mask_V5_lh_binarized.nii.gz \
#	-add $result_dir/group_mask_V5_rh_binarized.nii.gz \
#		$result_dir/group_mask_V5_binarized.nii.gz

#3dresample -master $result_dir/group_mask_V5_binarized_old.nii.gz -prefix $result_dir/group_mask_V5_lh_binarized.nii.gz -inset \
#	$result_dir/group_mask_V5_lh_binarized.nii
	
#3dresample -master $result_dir/group_mask_V5_binarized_old.nii.gz -prefix $result_dir/group_mask_V5_rh_binarized.nii.gz -inset \
#	$result_dir/group_mask_V5_rh_binarized.nii

fslmaths $result_dir/group_mask_V5_lh_binarized.nii.gz \
	-add $result_dir/group_mask_V5_rh_binarized.nii.gz \
		$result_dir/group_mask_V5_binarized.nii.gz
