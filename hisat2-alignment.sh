#!/usr/bin/bash

# Bash script for hisat2 alignment
# Align all .fastq.gz files in the Samples directory 
# Indexed reference genome provided in the Reference directory
# Output is generated in the SAM format in the Output directory
# stdout is saved into a Log directory
# All paths are relative to the base working directory 

#---------------------------------------------

# Define variables:
SAMPLE_DIR="./chrX_data/samples/"
REFERENCE_DIR="./chrX_data/indexes/"
OUTPUT_DIR="./analysis/hisat2-alignment/"
LOG_DIR="./analysis/log/hisat2/"

# Check if output directories exist or create it if it doesn't
[ -d $OUTPUT_DIR ] || mkdir -p $OUTPUT_DIR
[ -d $LOG_DIR ] || mkdir -p $LOG_DIR

# Get the pair of FASTQ files and run hisat2 alignment
# fileone and filetwo will conatain the full path 
for fileone in ${SAMPLE_DIR}*_1.fastq.gz; do
	
	# Replace _1 in fileone with _2 for the seconf FASTQ file
	filetwo=$(echo $fileone | sed 's/_1./_2./g')

	# Get the stem for the sample name using basename minus the _1.fastq.gz extension
	sample_name=$(basename $fileone _1.fastq.gz )
    
	hisat2 --dta -x "${REFERENCE_DIR}/chrX_tran" \
	-1 $fileone -2 $filetwo \
	-S "${OUTPUT_DIR}${sample_name}.sam" \
	--summary-file "${LOG_DIR}${sample_name}-alignment.log"

done

exit 

