#!/bin/bash

# activate micromamba environment. Comment out if running in other machine
eval "$(micromamba shell hook --shell bash)"
micromamba activate samba


# Set path to folder with the scripts used in the pipeline
SCRIPTS_FOLDER="$(dirname "$(realpath "$0")")"
PIPELINE_FOLDER="${SCRIPTS_FOLDER}/pipeline"

function display_usage() {
    echo "Usage: $0 -i accession_codes.tsv -o output_samba/"
    echo "      [-d|--download] [-a | --analyse] [-h|--help]"
    echo ""
    echo "Required arguemnts:"
    echo "  -i, --input          TSV file with two columns: project accessions and sample accessions from ENA."
    echo "  -o, --output         Path to folder where outputs will be saved."
    echo ""
    echo "Optional arguments:"
    echo "  -h, --help           Display this help message."
    echo ""
    exit 1
}

# Display usage if no required arguments were provided
if [[ $# -eq 0 ]]; then
    display_usage
fi


# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        -i|--input)
            input="$2"
            shift
            shift
            ;;
        -o|--output)
            output="$2"
            shift
            shift
            ;;
        -h|--help)
            display_usage
            ;;
        *)
	    echo "Unknown option: $1"
            display_usage
            ;;
    esac
done


# Step 3 -- Data analysis of each project
# counter of analysed projects
n=0

# analyse the downloaded data
cut -f2 $input | tail -n +2 | uniq | while read -r bioproject; do

    ((n++))
    echo "PROGRESS -- Analysing project $n: ${bioproject}" >> ${output}/logs/${bioproject}.log

    # Download ENA metada for the bioproject
    fastq-dl --accession $bioproject --outdir "$output/$bioproject/00.rawdata" --only-download-metadata --silent


    # lunch the analysis of the projects
    # while overall progress of the analysis goes to ${output}/nohups/${bioproject}.out,
    # step-specific logs can be found in ${output}/${bioproject}/nohups/
    bash ${SCRIPTS_FOLDER}/02.analyse_project.sh -s ${output}/${bioproject} --run_dada2 --run_picrust2 --run_modules >> ${output}/logs/${bioproject}.log 2>&1 &

    # save the PID of the process and add that to the log file to keep track of the analysis steps
    last_pid=$!
    wait "$last_pid"

    echo "Project ${accession_number} (number $n) finished\n"

done


# Step 4 -- Project integration
# one .out file is generated for each step of the following script
# bash ${SCRIPTS_FOLDER}/consolidate_projects.sh $output
