#!/usr/bin/bash

#------------------------------
# Bash script to estimate transcript abundance:
# Reference Annotation file provided as a stringtie merged GTF file 
# Sample input sorted BAM files
# Output is generated in a separate directory for each sample
# All paths are relative to the base working directory 
#------------------------------

#------------------------------
# Bash options and global variables:
# Abort on nonzero exit status
# Abort on unbound variable
# Don't hide errors within pipes
#------------------------------

set -e
set -u
set -o pipefail

declare -r SAMPLE_DIR="./analysis/hisat2-alignment/"
declare -r REFERENCE_ANNOTATION="./analysis/stringtie-assemble/stringtie-merged.gtf"
declare -r OUTPUT_DIR="./analysis/ballgown/"

#------------------------------
# Function main():
#------------------------------

main() {

  # Check if Annotation file exists
  if [[ ! -e "${REFERENCE_ANNOTATION}" ]]; then
    echo "Annotation file is missing in the provided path"
    exit 1
  fi

  # Loop through each BAM file in Sample directory
  for bamfile in ${SAMPLE_DIR}*.bam; do
    
   local sample_name=$(basename "${bamfile}" .sorted.bam)
   local sample_stem=$(basename "${sample_name}" _chrX)
   
   echo "${sample_stem}"
   
   stringtie -p 1 -e -B -G "${REFERENCE_ANNOTATION}" \
   -o "${OUTPUT_DIR}${sample_stem}/${sample_name}.gtf" \
   "${bamfile}"
    
  done
}

#------------------------------
# Call Main function
#------------------------------
main "${@}"



