#!/bin/bash -i

# this script will convert dicom to compress nifti and create according json files following the BIDS standard
# the script assumes that the DICOM files are already uncompressed inside the folder bids_dir/dicom/...

# --- options ----

# what kind of data should be converted?
anatomical=true
functional=true
fmap=true

# which participants should be converted?
participantNr_start=$1	 	
participantNr_end=$2

# define where the data is
#bids_dir=~/data/fMRI/timePath/bids/
#source_dir=~/data/fMRI/timePath/sourcedata/
bids_dir=/media/sf_data/fMRI/timePath/bids
source_dir=/media/sf_data/fMRI/timePath/sourcedata

# --- options ----

# go through the participants and convert data
for id in $(seq -f "%02g" $participantNr_start $participantNr_end)
do
	for ses_id in $(seq 1 2)
	do
		echo " "
		echo "---------- Converting DICOM dataset for subject: sub$id and session: $ses_id ----------"
		echo " "
		
		subjDir_bids=$bids_dir/sub-$id/ses-0$ses_id
		subjDir_source=$source_dir/sub-$id/ses-0$ses_id
		
		# anatomical data
		if [[ "$anatomical" == true ]]; then
			if [[ $ses_id == 1 ]]; then
				# T1 anat
				mkdir -p $subjDir_bids/anat
				dir_t1=($subjDir_source/*t1*)
				echo found T1 folder $dir_t1
				echo ----------- Converting T1 ----------- 
				dcm2niix -b y -z y y -o $subjDir_bids/anat -f "sub-"$id"_ses-0"$ses_id"_T1w" "$dir_t1"

				# T2 anat
				dir_t2=($subjDir_source/*t2*)
				echo found T2 folder $dir_t2
				echo -----------  Converting T2 ----------- 
				dcm2niix -b y -z y y -o $subjDir_bids/anat -f "sub-"$id"_ses-0"$ses_id"_T2w" "$dir_t2"
			fi
		fi

		# func data
		if [[ "$functional" == true ]]; then
			mkdir -p $subjDir_bids/func
			
			dir_func=($subjDir_source/*run*)
			run_id=1
			for dirname in "${dir_func[@]}" ; do		
				echo -----------  Converting functional data from folder $dirname ----------- 
				dcm2niix -b y -z y -o $subjDir_bids/func -f "sub-"$id"_ses-0"$ses_id"_task-class_run-0"$run_id"_bold" "$dirname"
				run_id=$((run_id + 1))	
			done
		fi

		# field maps
		if [[ "$fmap" == true ]]; then
			mkdir -p $subjDir_bids/fmap
			
			dir_fmap=($subjDir_source/*_mapping)
			for dirname in "${dir_fmap[@]}" ; do
				echo -----------  Converting field data from folder $dirname ----------- 
				dcm2niix -b y -z y -o $subjDir_bids/fmap -f "sub-"$id"_m_%s" "$dirname"
			done
		
			# find the file name for the magnitude file name
			magnitude_filename=$(find $subjDir_bids/fmap/ -type f -iname "*_e*.nii.gz" -exec basename {} ';')
			tmp="${magnitude_filename%%e*}"
		
			if [ "$tmp" != "$magnitude_filename" ]; then
				e_idx=${#tmp}
			fi
			magnitude_filename=${magnitude_filename:0:$e_idx-1}
		
			# write correct names for the fieldmaps and create correct .json file
			for fmapFileName in $(find $subjDir_bids/fmap/ -type f -iname "sub*" -exec basename {} ';')
			do	
				# magnitude files
				if [[ $fmapFileName == "$magnitude_filename"* ]]; then
					if [[ $fmapFileName == *".nii.gz" ]]; then
						# magnitude 1
						if [[ $fmapFileName == "$magnitude_filename.nii.gz" ]]; then
							echo "Renamed " $fmapFileName "to: magnitude1.nii.gz"
							mv $subjDir_bids/fmap/$fmapFileName $subjDir_bids'/fmap/sub-'$id'_ses-0'$ses_id'_magnitude1.nii.gz'
						# magnitude 2
						else
							echo "Renamed " $fmapFileName "to: magnitude2.nii.gz"
						mv $subjDir_bids/fmap/$fmapFileName $subjDir_bids'/fmap/sub-'$id'_ses-0'$ses_id'_magnitude2.nii.gz'
						fi
					# if the file is the .json file of the magnitudes, delete	
					else
						rm $subjDir_bids/fmap/$fmapFileName
					fi
				# phase files
				else
					if [[ $fmapFileName == *".nii.gz" ]]; then
						echo "Renamed " $fmapFileName "to: phasediff.nii.gz"
						mv $subjDir_bids/fmap/$fmapFileName $subjDir_bids'/fmap/sub-'$id'_ses-0'$ses_id'_phasediff.nii.gz'
					else
						echo "Renamed " $fmapFileName "to: phasediff.json"
						mv $subjDir_bids/fmap/$fmapFileName $subjDir_bids'/fmap/sub-'$id'_ses-0'$ses_id'_phasediff.json'
					fi
				fi	
			done
		
			# write necessary properties in the phasediff.json file (echotimes and for which files the fieldmap is intended)
			jq '. + { "EchoTime1": 0.00492 }' $subjDir_bids/fmap/sub-$id'_ses-0'$ses_id'_phasediff.json' > $subjDir_bids/fmap/tmp.json && mv $subjDir_bids/fmap/tmp.json $subjDir_bids/fmap/sub-$id'_ses-0'$ses_id'_phasediff.json'
			jq '. + { "EchoTime2": 0.00738 }' $subjDir_bids/fmap/sub-$id'_ses-0'$ses_id'_phasediff.json' > $subjDir_bids/fmap/tmp.json && mv $subjDir_bids/fmap/tmp.json $subjDir_bids/fmap/sub-$id'_ses-0'$ses_id'_phasediff.json'
			jq --arg subjId $id --arg ses_id $ses_id '. + { "IntendedFor": ["ses-0" + $ses_id + "/func/sub-" + $subjId + "_ses-0" + $ses_id + "_task-class_run-01_bold.nii.gz", "ses-0" + $ses_id + "/func/sub-" + $subjId + "_ses-0" + $ses_id + "_task-class_run-02_bold.nii.gz", "ses-0" + $ses_id + "/func/sub-" + $subjId + "_ses-0" + $ses_id + "_task-class_run-03_bold.nii.gz", "ses-0" + $ses_id + "/func/sub-" + $subjId + "_ses-0" + $ses_id + "_task-class_run-04_bold.nii.gz" ] }' $subjDir_bids/fmap/sub-$id'_ses-0'$ses_id'_phasediff.json' > $subjDir_bids/fmap/tmp.json && mv $subjDir_bids/fmap/tmp.json $subjDir_bids/fmap/sub-$id'_ses-0'$ses_id'_phasediff.json'
		fi
	done
	echo " "
		echo "---------- Done converting DICOM dataset for subject: sub$id ----------"
		echo " "
done
