#!/bin/bash

participantNr_start=$1	 	
participantNr_end=$2

data_dir=/Users/jAchtzehn/data/fMRI/timePath
#data_dir=/media/sf_data/fMRI/timePath



for id in $(seq -f "%02g" $participantNr_start $participantNr_end)
do
	
	result_dir=$data_dir/nilearn/sub-$id/space-T1w/masks
	
	fslmaths $result_dir/'sub-'$id'_hc_mask_binarized.nii.gz' \
		-add $result_dir/'sub-'$id'_ips_mask_binarized.nii.gz' \
			-add $result_dir/'sub-'$id'_ifg_mask_binarized.nii.gz' \
				$result_dir/'sub-'$id'_hc_ips_ifg_mask_binarized.nii.gz'
			
	echo 'Created' $result_dir/'sub-'$id'_hc_ips_ifg_mask_binarized.nii.gz' 
done