## Project: Fall 2023 CURE Course 
## Authors: R. Keating Godfrey
## Last updated: 23-11-05

## Resources:
## https://bioconductor.org/packages/release/bioc/html/DESeq2.html
## https://people.duke.edu/~ccc14/duke-hts-2017/Statistics/08032017/DESeq2-Notebook-introduction.html
## https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/DESeqDataSet-class


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

library(DESeq2)
library(dplyr)
library(ggplot2)

install.packages("htmltools")
library(htmltools)



### (1) Read in counts data ###

## OPTION A: read in single file with all counts 
## (compiled using Google Sheets or Excel)

## set working directory to source file location
setwd("~/Documents/2023/03_CURE/DESeq2/DGE_cure")

counts <- read.csv("counts_example.csv", header=T)


## OPTION B: read in counts from folder and combine into one dataframe ###

## change working director to location of counts files
setwd("/Users/rkeatinggodfrey/Documents/2023/03_CURE/DESeq2/DGE_cure/counts")

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

## we used the geneIDs as row names for importing and combining files,
## but now we need them to be a column for our analysis
## so we need to turn the gene id row names into a column (variable)
counts <-tibble::rownames_to_column(counts,"gene_id")



### (2) Import metadata file ###

## change directory back to the one where your metadata is located
setwd("/Users/rkeatinggodfrey/Documents/2023/03_CURE/DESeq2")

## Import metadata (colData) ###
metadata <- read.csv("/Users/rkeatinggodfrey/Documents/2023/03_CURE/DESeq2/DGE_cure/metadata_example.csv")

## **If you imported the directory to create the counts dataframe,
## add file names to metadata**
metadata$sample.name <-file.names

## select only the columns you need
meta <- metadata[,c(1,2,3)]
colnames(meta) <-c("sample","sex","body_part")

## make sure sex and body part are factors
meta$sex <-as.factor(meta$sex)
meta$body_part <-as.factor(meta$body_part)

## check that they are factors
is.factor(meta$sex)
is.factor(meta$body_part)
