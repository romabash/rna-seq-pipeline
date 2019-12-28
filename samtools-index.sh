#!/usr/bin/bash

# Bash script to index BAM files (.bai)
# Index all BAM files in the Samples directory 
# Output is generated in the same directory as a .bai file
# All paths are relative to the base working directory 

#---------------------------------------------

# Define variables:
SAMPLE_DIR="./analysis/hisat2-alignment/"

# Sort and Convert the SAM file from the sample directory, and delete afterwards
for samfile in ${SAMPLE_DIR}*.bam; do

	sample_name=$(basename $samfile )

    samtools index "${SAMPLE_DIR}${sample_name}" "${SAMPLE_DIR}${sample_name}.bai"

done
