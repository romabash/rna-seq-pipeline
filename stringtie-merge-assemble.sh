#!/usr/bin/bash

# Bash script to merge all of the transcrpts and estimate transcript abundance
# Reference Annotation file provided as GTF file 
# Input mergelist file containing all the GTF file names is genearted if it doesn't already exists 
# Output is generated in the same directory as a .gtf file
# All paths are relative to the base working directory 

#---------------------------------------------

# Define variables:
SAMPLE_DIR="./analysis/stringtie-assemble/"
REFERENCE_ANNOTATION="./chrX_data/genes/chrX.gtf"
OUTPUT_DIR="./analysis/stringtie-assemble/"

# Check to see if mergelist.txt exists, if not create it with the names of all the GTF files
if [ ! -e "${SAMPLE_DIR}mergelist.txt" ]; then
    
	for f in ${SAMPLE_DIR}*.gtf; do 
		echo $f >> "${SAMPLE_DIR}mergelist.txt"
	done

fi

#stringtie --merge -p 8 -G chrX_data/genes/chrX.gtf -o stringtie_merged.gtf chrX_data/mergelist.txt
