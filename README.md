# TidyTWAS
Harmonizing, processing, and visualizing cross-trait TWAS output

This repository contains functions for harmonizing, processing, and visualizing cross-trait TWAS output.

TidyTWAS is designed to work with the output of https://github.com/gusevlab/fusion_twas/blob/master/FUSION.assoc_test.R .

## Prerequisites
R and the required packages:
install.packages(c('optparse','tidyverse','data.table','xlsx','purrr','plyr','stringr','biomaRt','sqldf'))

Performing TWAS using FUSION
Instructions on how to perform a TWAS are available at http://gusevlab.org/projects/fusion/

Per chromosome results files should be combined without duplicating the header.
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

## 

## Help
For questions or comments, please get in touch with me at lynsey.hall@yahoo.co.uk 
