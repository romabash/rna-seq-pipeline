#-------------------
# Differential expression analysis of RNA-Seq data using Ballgown
# Set working directory to base directory of this R script
#-------------------

# Load packages
library(ballgown)
library(RColorBrewer)
library(readr)
library(genefilter)
library(dplyr)
library(devtools)

#-------------------
# Read-in Phenotype Data:
# Contains info on Gender and Population of Samples
#-------------------

pheno_data = read.csv("chrX_data/geuvadis_phenodata.csv")

#-------------------
# Read-in the expression data calculated by StringTie.  Provide following arguments:
# directory where the data are stored ( dataDir = “ballgown”),
# pattern that appears in the sample names ( samplePattern = "ERR"),
# phenotypic information loaded in the previous step ( pData )
# Output is a dataframe (Bioconductor S4 object)
#-------------------

bg_chrX = ballgown::ballgown(dataDir="analysis/ballgown", samplePattern="ERR", pData=pheno_data)

#-------------------
# Filter to remove low abundance genes using ballgown::subset function
# Remove all transcripts with a variance across samples less than one
# Arguments:
# a ballgown object (x)
# condition on which to subset (cond, provided as a string)
# genomesubset (TRUE by default)
#-------------------

bg_chrX_filt = subset(bg_chrX, "rowVars(texpr(bg_chrX)) > 1", genomesubset=TRUE)





