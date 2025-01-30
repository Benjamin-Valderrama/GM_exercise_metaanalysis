#!/bin/bash

# activate environment
# eval "$(micromamba shell hook --shell bash)" ; micromamba activate samba

# read input file in (has 'bioproject' and 'run accession' codes)
input=$1
# folder where all the analysed bioprojects will be stored
output=$2


# folder with the nohups for all the projects
mkdir -p ${output}/nohups


# STEP 1 -- Data download:
# counter of analysed projects
n=0

# read the file line by line, skipping the header
tail -n +2 "$input" | while IFS=$'\t' read -r run_accession bioproject; do

    # if it is a new bioproject...
    if [ ! -d $output/$bioproject ]; then

        # keep track of analysed projects
        ((n++))

        # create a directory for the bioproject if it doesn't already exist
        mkdir -p "$output/$bioproject/00.rawdata"
        mkdir -p "$output/$bioproject/logs"
        mkdir -p "$output/$bioproject/outputs"

    fi

    # download data
    echo "PROGRESS -- Downloading raw data of project $n: ${bioproject}" > ${output}/logs/${bioproject}.log
    echo "Downloading : $run_accession" >> $output/$bioproject/logs/download.log

    fastq-dl --accession $run_accession --outdir $output/$bioproject/00.rawdata --silent

    # download one file at a time:

    # Comment: I want the script to behave differently but I don't know how to do it:
    # I want all the run accessions of the same project to be downloaded at once, and
    # all files of the same project have to be downloaded before starting the download
    # of samples in a new project.

    last_pid=$!
    wait $last_pid

done


