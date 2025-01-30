#!/bin/bash

# This script is just to keep track of the steps used. Not that it's actually needed.

messy_input=$1
output_folder=$2

# cut -d, -f3,5 works because our messy_input has the information of interest in those columns.
# tr ',' '\t' transform the csv file to a tsv
# grep "\S" show non-white space characters in the file (used to remove blank lines)
cat ${messy_input} | cut -d, -f3,5 | tr ',' '\t' | grep "\S" > ${output_folder}/clean_data_to_download.tsv

