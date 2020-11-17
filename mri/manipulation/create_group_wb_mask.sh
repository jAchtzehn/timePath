#!/bin/bash

participantNr_start=$1	 	
participantNr_end=$2

data_dir=/Users/jAchtzehn/data/fMRI/timePath
fmriprep_dir=$data_dir/fmriprep/
nilearn_dir=$data_dir/nilearn/
output_dir=$nilearn_dir/group_masks/space-MNI152NLin2009cAsym/
mkdir -p $output_dir

3dMean -prefix $output_dir/group_mask_hc.nii.gz ${nilearn_dir}/sub-*/space-MNI152NLin2009cAsym/masks/sub-*_hc_mask_binarized.nii.gz

mri_binarize --i $output_dir/group_mask_hc.nii.gz --o $output_dir/group_mask_hc_binarized.nii.gz --min 0.5