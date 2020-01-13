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

###############################
# Data analysis: Part 1
###############################

#-------------------
# Read-in Phenotype Data:
# Contains info on Gender and Population of Samples
#-------------------

pheno_data <- read.csv("chrX_data/geuvadis_phenodata.csv")

#-------------------
# Read-in the expression data calculated by StringTie.  Provide following arguments:
# dataDir: directory where the data are stored (“ballgown”),
# samplePattern: pattern that appears in the sample names ("ERR"),
# pData: phenotypic information loaded in the previous step
# Output is a dataframe (Bioconductor S4 object)
#-------------------

bg_chrX <- ballgown::ballgown(dataDir="analysis/ballgown", samplePattern="ERR", pData=pheno_data)

#-------------------
# Filter to remove low abundance genes using ballgown::subset function
# Remove all transcripts with a variance across samples less than one
# Provide following Arguments:
# x: a ballgown object
# cond: condition on which to subset (provided as a string)
# genomesubset: subset x to a specific part of the genome (TRUE by default)
#-------------------

bg_chrX_filt <- ballgown::subset(x=bg_chrX, cond="rowVars(texpr(bg_chrX)) > 1", genomesubset=TRUE)

#-------------------
# Identify transcripts that show statistically significant differences between groups.
# Account for variation in expression due to other variables.
# Look for transcripts differentially expressed between sexes
# Correct for any differences in expression due to the population variable.
# Use the stattest function from Ballgown:
# Set the getFC=TRUE parameter to look at the confounder-adjusted fold change between the two groups.
# Ballgown’s statistical test is a standard linear model based comparison,
# For small sample sizes (n < 4 per group) it is often better to perform regularization.
# This can be done using the limma package in Bioconductor.
# Other regularized methods such as DESeq and edgeR can be applied to gene or exon counts, 
# but are not appropriate for direct application to FPKM abundance estimates.
# The statistical test uses a cumulative upper quartile normalization
#
# Provide following Arguments:
# x: ballgown object
# feature: the type of genomic feature to be tested for DE ("gene", "transcript", "exon", "intron")
# covariate: string representing the name of the covariate of interest for the DE tests
# adjustvars: optional vector of strings representing the names of potential confounders
# getFC: if TRUE, also return estimated fold changes, adjusted for library size and confounders, between populations
# meas: the expression measurement to use for statistical tests ("cov", "FPKM", "rcount", "ucount", "mrcount", "mcov")
#-------------------

results_transcripts <- ballgown::stattest(gown=bg_chrX_filt, 
                                         feature="transcript", 
                                         covariate="sex", 
                                         adjustvars=c("population"), 
                                         getFC=TRUE, meas="FPKM")

#-------------------
# Identify genes that show statistically significant differences between groups;
# Run the same function used to identify differentially expressed transcripts, with feature=“gene” 
#-------------------

results_genes <- ballgown::stattest(gown=bg_chrX_filt, 
                                   feature="gene",
                                   covariate="sex", 
                                   adjustvars = c("population"), 
                                   getFC=TRUE, meas="FPKM")

#-------------------
# Add gene names and gene IDs to the results_transcripts data frame:
# ballgown::geneNames returns named vector of gene names included in the ballgown object,
# named and ordered by corresponding numeric transcript ID
#-------------------

results_transcripts <- data.frame(geneNames=ballgown::geneNames(x=bg_chrX_filt), 
                                 geneIDs=ballgown::geneIDs(x=bg_chrX_filt), 
                                 results_transcripts)

#-------------------
# Sort the results from smallest p-value to largest:
#-------------------

results_transcripts <- dplyr::arrange(results_transcripts, pval)
results_genes <- dplyr::arrange(results_genes, pval)

#-------------------
# Save result to a CSV file
#-------------------

readr::write_csv(x=results_transcripts, path="analysis/results_transcripts.csv")
readr::write_csv(x=results_genes, path="analysis/results_genes.csv")

#-------------------
# Identify transcripts and genes with a q-value of less than 0.05:
# Chromosome X has 9 transcripts that are differentially expressed between the sexes, 
# At the gene level, chromosome X has 8 differentially expressed genes
#-------------------

de_transcript <- ballgown::subset(x=results_transcripts, results_transcripts$qval < 0.05)
de_genes <- ballgown::subset(x=results_genes, results_genes$qval < 0.05)

###############################
# Data visualization: Part 2
###############################

# Set color palette
palette(rainbow(6))

#-------------------
# Show the distribution of gene abundances (measured as FPKM values) across samples, colored by sex.
# Compare the FPKM measurements for the transcripts. 
# Extract transcript-level expression measurements from ballgown objects using texpr function
# The plots will be easier to visualize if  FPKM data is transformed.
# Use a log 2 transformation that adds one to all FPKM values. The third command actually creates the plot.
#-------------------

fpkm = texpr(x=bg_chrX, meas="FPKM")
fpkm = log2(fpkm+1)

boxplot(fpkm, col=as.numeric(pheno_data$sex), las=2, ylab='log2(FPKM+1)')





