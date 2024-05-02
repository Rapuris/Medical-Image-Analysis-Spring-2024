#!/bin/bash

# Define the FastSurfer home directory
export FASTSURFER_HOME="~/Desktop/fast/FastSurfer"

# Loop through the subjects and run FastSurfer for each one
for i in {1..50}; do
    # Define input and output directories
    input="./MRI_dataset/t1w/subject_${i}_t1w.nii.gz"
    outputDir="./MRI_dataset/t1w_segs/"

    # Run FastSurfer with the specific subject id (sid)
    $FASTSURFER_HOME/run_fastsurfer.sh --t1 $input \
                                       --sd $outputDir \
                                       --sid $i \
                                       --seg_only
done

