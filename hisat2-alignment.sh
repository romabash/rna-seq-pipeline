#!/usr/bin/bash

# Bash script for hisat2 alignment
# Align all .fastq.gz files in the Samples directory 
# Indexed reference genome provided in the Reference directory
# Output is generated in the SAM format in the Output directory
# stdout is saved into a Log directory
# All paths are relative to the base working directory 

#---------------------------------------------

# Define variables:
REFERENCE_DIR="./chrX_data/indexes/"
SAMPLE_DIR="./chrX_data/samples/"
OUTPUT_DIR="./hisat2-alignment/"
LOG_DIR="./log/hisat2/"

# Check if output directories exist or create it if it doesn't
[ -d $OUTPUT_DIR ] || mkdir -p $OUTPUT_DIR
[ -d $LOG_DIR ] || mkdir -p $LOG_DIR

# Get paired files from the sample directory
#for fileone in ${SAMPLE_DIR}*_1.fastq.gz; do
#	filetwo=$(echo $fileone | sed 's/_1./_2./g')
#	sample_name=$(basename $fileone _1.fastq.gz )
#	echo $sample_name
#	echo $fileone
#	echo $filetwo
#	ls -ls $fileone	
#	ls -ls $filetwo
#done


for fileone in ${SAMPLE_DIR}*_1.fastq.gz; do
	
	filetwo=$(echo $fileone | sed 's/_1./_2./g')
	sample_name=$(basename $fileone _1.fastq.gz )
    
	hisat2 --dta -x "${REFERENCE_DIR}/chrX_tran" \
	-1 $fileone -2 $filetwo \
	-S "${OUTPUT_DIR}${sample_name}.sam" \
	--summary-file "${LOG_DIR}${sample_name}-alignment.log"

done

# Working

