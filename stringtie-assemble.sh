#!/usr/bin/bash

# Bash script to assemble transcripts using stringtie
# Reference Annotation file provided as GTF file 
# Input BAM files in the Samples directory 
# Output is generated in a separate directory as a .gtf file
# All paths are relative to the base working directory 

#---------------------------------------------

# Define variables:
SAMPLE_DIR="./analysis/hisat2-alignment/"
REFERENCE_ANNOTATION="./chrX_data/genes/chrX.gtf"
OUTPUT_DIR="./analysis/stringtie-assemble/"

# Check if output directory exist or create it if it doesn't
[ -d $OUTPUT_DIR ] || mkdir -p $OUTPUT_DIR

# Loop through each BAM file in Sample directory
for bamfile in ${SAMPLE_DIR}*.bam; do

	sample_name=$( basename $bamfile .sorted.bam)
	
	# -p indicates number of threads to use (default: 1)
	# -G indicates reference annotation to use (GTF/GFF)
	# -o indicates output file name for assembled transcripts GTF (default: stdout)
	# -l indicates name prefix for output transcripts (default: STRG)
	# last, provide input BAM file
    stringtie -p 1 -G "${REFERENCE_ANNOTATION}" -o "${OUTPUT_DIR}${sample_name}.gtf" "${bamfile}"

done

exit
