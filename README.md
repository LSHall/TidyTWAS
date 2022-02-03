### TidyTWAS
This repository contains functions for harmonizing, processing, and visualizing cross-trait transcriptome-wide association analysis output.

##TidyTWAS offers three functions:
1) to read in the results from a TWAS performed using FUSION software for up to five traits, to harmonize and process these results, and write out a spreadsheet of cross-trait TWAS Z-scores and P-values for TWAS-significant genes
2) to write out locus files for downstream TWAS finemapping
2) to generate a heatmap of TWAS associations for genes which are TWAS significant in two or more traits

## Prerequisites
#R and the required packages:
install.packages(c('optparse','tidyverse','data.table','xlsx','purrr','plyr','stringr','biomaRt','sqldf'))

#Performing TWAS using FUSION
TidyTWAS is designed to work with the output of https://github.com/gusevlab/fusion_twas/blob/master/FUSION.assoc_test.R
Instructions on how to perform a TWAS are available at http://gusevlab.org/projects/fusion/

#FUSION output to be in one file containing chromosomes 1 - 22

For example
# insert my code

## Parameters
Flag | Description | Default
--twas | Path to genome-wide TWAS results | NA
--output | Path to save output | NA 

## Output files
# Loci text files
# Mirrored (meta-analysed?) manhattan plot
# Cross-trait heatmap of TWAS Z-scores

## Examples


## Help
For questions or comments, please get in touch with me at lynsey.hall@yahoo.co.uk 
