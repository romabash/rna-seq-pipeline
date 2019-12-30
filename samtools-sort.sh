#!/usr/bin/bash

# Bash script to sort SAM files and convert to BAM
# Conver and Sort all SAM files in the Samples directory 
# Output is generated in the same directory as a sorted BAM 
# All paths are relative to the base working directory 

#---------------------------------------------

# Define variables:
SAMPLE_DIR="./analysis/hisat2-alignment/"

# Sort and Convert the SAM file from the sample directory, and delete afterwards
for samfile in ${SAMPLE_DIR}*.sam; do

	# Get the stem of the samfile to convert to bam 
	sample_name=$(basename $samfile .sam )

    samtools sort "${SAMPLE_DIR}${sample_name}.sam" -o "${SAMPLE_DIR}${sample_name}.sorted.bam"
	rm "${SAMPLE_DIR}${sample_name}.sam"

done

exit
