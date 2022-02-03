### TidyTWAS
This repository contains functions for harmonizing, processing, and visualizing cross-trait transcriptome-wide association analysis output.

##TidyTWAS offers three functions:
1) to read in the results from a TWAS performed using FUSION software for up to five traits, to harmonize and process these results, and write out a spreadsheet of cross-trait TWAS Z-scores and P-values for TWAS-significant genes
2) to write out locus files for downstream TWAS finemapping
2) to generate a heatmap of TWAS associations for genes which are TWAS significant in two or more traits

*** this is at a very primordial stage and a work in progress ***
The TidyTWAS.R script will end up being the final script that will be runnable from the command line (similar to FUSION.assoc_test.R)
TidyTWAS.R is taking information from:
1) TWAS.results.processing.R for reading in data from multiple sources and harmonizing
2) deriving_sig_loci_for_conditional_TWAS.R for generating loci for downstream finemapping
3) XDis_heatmap.R for generating the cross-trait heatmap

Mirrored manhattan produces a mirrored manhattan plot of TWAS results for one trait~tissue association, I am figuring out whether this is a feature I want include in TidyTWAS or whether it is unnecessary, as the main function of this software will be for multiple phenotype/tissue results.


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
