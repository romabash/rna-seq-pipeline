#!/usr/bin/bash

# Bash script to merge all of the transcrpts
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

# Merge transcript files using stringtie merge
# -p indicates number of threads to use (default: 1)
# -G indicates reference annotation to include in the merging (GTF/GFF3)
# -o indicates output file name for the merged transcripts GTF (default: stdout)
# -l indicates name prefix for output transcripts (default: MSTRG)
# last, provide generated mergelist text filw with names of all GTF files to merge

stringtie --merge -p 1 -G "${REFERENCE_ANNOTATION}" \
-o "${OUTPUT_DIR}stringtie-merged.gtf" \
"${SAMPLE_DIR}mergelist.txt"

















