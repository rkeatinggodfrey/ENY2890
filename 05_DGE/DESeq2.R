## Project: Fall 2023 CURE Course 
## Authors: R. Keating Godfrey
## Last updated: 23-10-18

## Resources:
## https://bioconductor.org/packages/release/bioc/html/DESeq2.html
## https://people.duke.edu/~ccc14/duke-hts-2017/Statistics/08032017/DESeq2-Notebook-introduction.html
## https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/DESeqDataSet-class

## Install the packages you need to do the analysis
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)
library(dplyr)
library(ggplot2)

install.packages("htmltools")
library(htmltools)

##############################################################
### (1) Read in counts data and combine into one dataframe ###
##############################################################
## Create an object out of the filepath to the count .tsv files

## change working director to location of counts files
## you can also go to Session > Set working directory > Source file location
## to set directory to the folder where this R script is saved
setwd("{filepath to folder on your computer with your counts.tsv files}")

## create a list object that contains all of the .tsv files in this directory
filedir <-list.files(pattern = ".tsv$")

## Import all files as a list of dataframes
all_files <- lapply(filedir,function(x) {
  read.table(file = x, 
             sep = '\t', 
             header = F,
             row.names=1)
})

## Combine into single dataframe
counts <-bind_cols(all_files)

## create a string of file names
file.names<-c(filedir[1], filedir)

## why does it duplicate the first file name? Remove it!
file.names <- file.names[2:26]

## now get rid of the ".tsv" in the file name
file.names <- gsub(".tsv", "", file.names)

## name the columns of your counts dataframe to match your metadata file
colnames(counts)<-c(file.names)

## change directory back to the one where your metadata is located
## you can also go to Session > Set working directory > Source file location
## to set directory to the folder where this R script is saved
setwd("/Users/rkeatinggodfrey/Documents/2023/03_CURE/DESeq2")

## we used the geneIDs as row names for importing and combining files,
## but now we need them to be a column for our analysis
## so we need to turn the gene id row names into a column (variable)
counts <-tibble::rownames_to_column(counts,"gene_id")

