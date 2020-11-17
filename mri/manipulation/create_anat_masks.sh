#!/bin/bash

### 
# this script creates label files from the parcellation annotation of freesurfer and converts those to volumes that may be used as 
# anatomical masks
# by changing the input properties of mri_label2vol "--label", other labels can be added/removed
###

participantNr_start=$1	 	
participantNr_end=$2

data_dir=/Users/jAchtzehn/data/fMRI/timePath
#data_dir=/media/sf_data/fMRI/timePath

mask_name='cuneus' # name for the file
label_name=cuneus # label name is in the folder .../freesurfer/labels
fill_th=.1

for id in $(seq -f "%02g" $participantNr_start $participantNr_end)
do
	
	result_dir=/Users/jAchtzehn/data/fMRI/timePath/nilearn/sub-$id/space-T1w/masks
	mkdir -p $result_dir
	
	# create labels from the annotations
	#mri_annotation2label --subject sub-$id --hemi lh --outdir $data_dir/freesurfer/sub-$id/labels
	#mri_annotation2label --subject sub-$id --hemi rh --outdir $data_dir/freesurfer/sub-$id/labels
	
	#mri_annotation2label --subject sub-$id --hemi lh --annotation aparc.a2009s --outdir $data_dir/freesurfer/sub-$id/labels
	#mri_annotation2label --subject sub-$id --hemi rh --annotation aparc.a2009s --outdir $data_dir/freesurfer/sub-$id/labels
	
	label_dir=$data_dir/freesurfer/sub-$id/labels		# directory in which the labels are stored
	
	# first create a transform so we can transform the volume into native space, using the already existing HC segmentation as reference
	tkregister2 --mov $result_dir/'sub-'$id'_hc_mask_binarized.nii.gz' --noedit --s sub-$id --regheader --reg $data_dir/freesurfer/sub-$id/mri/register.dat
	
	# --------------------------------------------
	# now create the volumes for both hemispheres
	mri_label2vol --label $label_dir/lh.$label_name.label --label $label_dir/rh.$label_name.label \
			--temp $result_dir/'sub-'$id'_hc_mask_binarized.nii.gz' \
				--subject sub-$id --fillthresh $fill_th --o $data_dir/freesurfer/sub-$id/mri/$mask_name'_mask.nii.gz' --reg $data_dir/freesurfer/sub-$id/mri/register.dat
	
	# binarize image to create a mask
	mri_binarize --i $data_dir/freesurfer/sub-$id/mri/$mask_name'_mask.nii.gz' --o $result_dir/'sub-'$id'_'$mask_name'_mask_binarized.nii.gz' --min 0.0001
	
	# --------------------------------------------
	# now create the volumes for right hemispheres
	mri_label2vol --label $label_dir/rh.$label_name.label \
		--temp $result_dir/'sub-'$id'_hc_mask_binarized.nii.gz' \
			--subject sub-$id --fillthresh $fill_th --o $data_dir/freesurfer/sub-$id/mri/$mask_name'_rh_mask.nii.gz' --reg $data_dir/freesurfer/sub-$id/mri/register.dat
	
	# binarize image to create a mask
	mri_binarize --i $data_dir/freesurfer/sub-$id/mri/$mask_name'_rh_mask.nii.gz' --o $result_dir/'sub-'$id'_'$mask_name'_rh_mask_binarized.nii.gz' --min 0.0001
	
	# --------------------------------------------
	# now create the volumes for left hemispheres
	mri_label2vol --label $label_dir/lh.$label_name.label \
		--temp $result_dir/'sub-'$id'_hc_mask_binarized.nii.gz' \
			--subject sub-$id --fillthresh $fill_th --o $data_dir/freesurfer/sub-$id/mri/$mask_name'_lh_mask.nii.gz' --reg $data_dir/freesurfer/sub-$id/mri/register.dat
	
	# binarize image to create a mask
	mri_binarize --i $data_dir/freesurfer/sub-$id/mri/$mask_name'_lh_mask.nii.gz' --o $result_dir/'sub-'$id'_'$mask_name'_lh_mask_binarized.nii.gz' --min 0.0001
done
