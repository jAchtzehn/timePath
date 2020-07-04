#!/bin/bash

echo "---> MRIQC: $1, FMRIPREP: $2, Participants: $3"

# --- options ----

# what kind of preproc steps should run?
run_mriqc=$1
run_fmriprep=$2

# define where the data is
#data_dir=~/Data/fMRI/pilot_classification/
#data_dir=/Volumes/Serengeti/tmpDir
data_dir=/media/sf_data/fMRI/timePath/

echo "Input data directory: $data_dir"
echo " "
echo " "

# mriqc specific
mriqc_version=0.14.2			# which mriqc version to use
nr_cpus=7						# how many CPUs (real & virtual) should be used?
movement_threshold=1			# threshold in mm for movement estimation	

# fmriprep specific
fmriprep_version=1.1.6			# which mriqc version to use
nr_threads=14					# how many threads should be used?


# which participants should be processed?
part_labels=$3					# participants to run

# --- options ----

# clean data_dir from any .files 
find $data_dir/bids -name ".*" -type f -delete

# ---- mriqc ----
if [[ "$run_mriqc" == true ]]; then
	echo "------ Running MRIQC on participant level ------"
	# 1. run on individual participants
	docker run -it --rm \
	-v $data_dir/bids:/data:ro \
	-v $data_dir/mriqc:/out \
	poldracklab/mriqc:$mriqc_version /data /out participant --participant_label $part_labels \
	--no-sub -vv --nprocs=$nr_cpus --verbose-reports --fd_thres $movement_threshold --mem_gb 32

	echo "------ Running MRIQC on group level ------"
	# 2. run group analysis
	docker run -it --rm \
	-v $data_dir/bids:/data:ro \
	-v $data_dir/mriqc:/out \
	poldracklab/mriqc:$mriqc_version /data /out group \
	--no-sub -vv --nprocs=$nr_cpus --verbose-reports --fd_thres $movement_threshold --mem_gb 32
fi

# ----fmriprep, run with output dir on VM (symbolic linking does not work otherwise) ----
if [[ "$run_fmriprep" == true ]]; then
	echo "------ Running FMRIPREP on participant level, without freesurfer, with fieldmap correction ------"
	docker run -it --rm \
	-v $data_dir/bids:/data:ro \
	-v $data_dir:/out \
	-v $FREESURFER_HOME/license.txt:/opt/freesurfer/license.txt \
	poldracklab/fmriprep:$fmriprep_version /data /out participant --participant_label $part_labels\
	--nthreads $nr_threads --output-space T1w template --work-dir $data_dir/scrtch \
	--fs-no-reconall --mem-mb 32000 --resource-monitor
fi

echo "---> Summary: MRIQC: $1, FMRIPREP: $2, Participants: $3"
echo "---> Summary: Total runtime: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"

